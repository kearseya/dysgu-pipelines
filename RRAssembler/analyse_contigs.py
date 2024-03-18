"""
Map contigs using LAST to reconstruct the mapping arrangement.

Generates a Best_alignments file

After this run the remove_common_sites.py script
"""


import pandas as pd
from collections import defaultdict, Counter
from subprocess import run, call, Popen, PIPE

import numpy as np
import sys
import os
import time

import glob
import best_path  # Cython extention
import find_path2  # Python version for debugging

#from ctypes import *
#libbest_path = cdll.LoadLibrary("./best_path.so")

def load_contigs(input_contigs, min_cont_length):
    print("Input", input_contigs)
    contigs = defaultdict(list)
    for item in input_contigs:
        lines = open(item, "r").readlines()
        cont = None
        for l in lines:
            if l[0] == ">":
                #samp = item.split("/")[-2]
                #cont = samp + "_" + l.strip()[1:].replace(" ", ":")
                cont = l[1:].strip()
            else:
                contigs[cont].append(l.strip())
    conts = {k: "".join(v) for k, v in contigs.items() if sum(len(i) for i in v) > min_cont_length}

    return conts


def dodi_contigs(pfix, ref, outdir="out_data", procs=12): #fai, l_match_score, l_mismatch_cost, l_gap_open, l_gap_extend):
    #if os.path.exists(f"{ref}.amb") or os.path.exists(f"{ref}.ann") or os.path.exists(f"{ref}.pac") == False:
    #    print("Indexing reference")
    #    call(f"bwa index {ref}", shell=True)
    print("Running dodi")
    all_contigs_file = os.path.join(outdir, "all_contigs."+pfix+".fq")
    dodi_bam_file = os.path.join(outdir, pfix+".bwa_dodi.bam")
    dodi_psl_file = os.path.join(outdir, pfix+".bwa_dodi.psl")
    file_path = os.path.dirname(os.path.realpath(__file__))
    call(f"{file_path}/./dodi_pipe.sh {all_contigs_file} {ref} {dodi_bam_file}", shell=True)
    #run(f"cat {all_contigs_file} | bwa mem -c 1000 -A2 -B3 -O5 -E2 -T0 -L0 -D 0.25 -r 1.25 -d 200 -k 11 -a -t{procs} {ref} - | dodi --paired False -c 1 -u 21 --ol-cost 2 --max-overlap 50000 - | samtools view -bh - | samtools sort -o {dodi_bam_file}; samtools index {dodi_bam_file}", shell=True)
    call(f"samtools view -h {dodi_bam_file} | ./sam2psl -h | grep -v ^# > {dodi_psl_file}", shell=True)

    # print("Reference database", ref)
    # call("lastal -r {r} -q {q} -a {a} -b {b} -D1000 -P6 -K3 -C3 {ref} all_contigs.{pfix}.fa > last_output.{pfix}.maf".
    #      format(pfix=pfix,
    #             ref=ref,
    #             r=l_match_score,
    #             q=l_mismatch_cost,
    #             a=l_gap_open,
    #             b=l_gap_extend), shell=True)

    # call("maf-convert psl last_output.{p}.maf > last_output.{p}.psl".format(p=pfix), shell=True)

    # # Convert to bam for viewing also
    # commands = [
    #             "maf-convert sam {f}.maf > {f}.sam",
    #             "samtools view -bt {fai} {f}.sam > {f}.bam",
    #             "samtools sort -o {f}.aln.bam {f}.bam",
    #             "samtools index {f}.aln.bam"]

    # keys = {"f": "last_output.{p}".format(p=pfix),
    #         "fai": fai}
    # print(keys)
    # for com in commands:
    #     print(com.format(**keys))
    #     call(com.format(**keys), shell=True)


def collect_mappings(pfix, mapper, split_by_sample=True, contig_name=False,
                     match=1, mis_match=1, gap_open=4, gap_extend=2, min_identity=0.85, outdir="out_data"):

    dodi_psl_file = os.path.join(outdir, pfix+".bwa_dodi.psl")
    #lines = [i.strip().split("\t") for i in open(dodi_psl_file, "r")
    #         if "---" not in i and "Layout" not in i and i != "\n"]
    #if len(lines) == 1:
    #    raise IOError("No contigs found")

    #columns = map(str, "match mis-match rep_match Ns Q_gap_count Q_gap_bases T_gap_count T_gap_bases strand Q_name Q_size Q_start Q_end T_name T_size T_start T_end block_count blockSizes qStarts tStarts".split(" "))

    #['matches', 'misMatches', 'repMatches', 'nCount', 'qNumInsert', 'qBaseInsert', 'tNumInsert', 'tBaseInsert', 'strand', 'qName', 'qSize', 'qStart', 'qEnd', 'tName', 'tSize', 'tStart', 'tEnd', 'blockCount', 'blockSizes', 'qStarts', 'tStarts', 'qBlocks', 'tBlocks', 'RI', 'RNEXT', 'PNEXT', 'TLEN', 'MAPQ', 'AS', 'MAPS', 'FPAIRED', 'FPROPER_PAIR', 'FSECONDARY', 'FQC', 'FDUP']
    #df = pd.DataFrame(data=lines[2:], columns=columns)
    #df.columns = [i.strip() for i in df.columns]
    df = pd.read_csv(dodi_psl_file, sep="\t")

    df.dropna(inplace=True)  #
    df = df.apply(pd.to_numeric, errors='ignore')  # Convert to numeric where possible

    for key in ["blockSizes", "qStarts", "tStarts"]:
        df[key] = [[int(j) for j in i.split(",") if j != ""] for i in list(df[key])]
    df = df.sort_values("qStarts")

    # Drop un-needed clumns
    #df.drop(["Ns", "T_size", "rep_match", "Q_size"], 1, inplace=True)


    #
    df["matches"] = df["matches"] + df["misMatches"]  # Ignore mismatches

    df["Identity"] = (df["matches"] - df["misMatches"]) / df["matches"]
    df = df[df["Identity"] > min_identity]

    # Default score scheme is +1 alignment -1 mismatch, -4 gap open, -2 gap extend for LAST mapper
    df["Last_score"] = ((match * df["matches"]) - (mis_match * df["misMatches"]) -
                        (gap_open * df["tNumInsert"]) - (gap_extend * df["tNumInsert"]))

    # Drop matches of less than 20
    df = df[df["matches"] > 20]

    if split_by_sample:
        # Add a sample column and rename the reads to their original names
        samples, names = zip(*[i.split("_", 1) for i in list(df["qName"])])
        df["sample"] = samples
        df["contig_names"] = df["qName"]  # Keep a reference to old name
        df["Q_name"] = names

        if contig_name:  # Filter out a single contig for debugging
            df = df[df.Q_name == contig_name]
    else:
        df["sample"] = ["None"] * len(df)

    if len(df) == 0:
        raise IOError("No contigs found")
    return df


def slice_range(change_points, n_alns):
    # Generate slice indices based on a list of change points.
    if len(change_points) == 1:
        # Only one contig, return whole list
        yield (0, n_alns)
    else:
        last_seen = 0
        for item in change_points:
            yield (last_seen, item + 1)
            last_seen = item + 1


#if __name__ == "__main__":
def main(regions_dir, ref="/scratch/scwc0010/hg38.fa", outdir="out_data"):
    # PARAMETERS ---------------------------------------------------------
    # INPUT
    input_contigs = glob.glob(f"{regions_dir}/*/contigs.filtered.fasta")  # All
    pfix = "all"

    # input = glob.glob("random_regions_to_assemble/contigs/*.fasta")
    # pfix = "random"

    # DODI
    run_dodi = True
    #ref = "~/Documents/Data/db/hg19_last_ref/hg19last"
    #fai = "/Users/kezcleal/Documents/Data/db/ucsc.hg19.fasta.fai" #"/Users/kezcleal/Documents/Data/db/hg38/hg38.fa.fai"
    l_match_score = 1       # Last default 1
    l_mismatch_cost = 4     # 1
    l_gap_open = 6          # 7
    l_gap_extend = 1        # 1

    min_cont_length = 150

    # Other
    target = False  #"NODE_42_length_315_cov_0.240196"          # select target contig, enter string
    inter_cost = 15  # 20
    intra_cost = 14  # 10
    ins_cost = 1            # Score goes down quicker with an insertion
    hom_cost = 3            # A match gives + 2 due to both seqs being counted
    min_identity = 0.86

    keep_single_alignments = True  # Contigs with no SVs have only one alignment. Keep these?

    # --------------------------------------------------------------------

    np.random.seed(12345)

    print("="*40)
    print(pfix)
    print("="*40)

    split = False
    contigs = load_contigs(input_contigs, min_cont_length)
    
    all_contigs_file = os.path.join(outdir, "all_contigs."+pfix+".fq")
    print("Contig names:", list(contigs.keys())[0:3], "etc")
    if run_dodi:
        print("Mapping {} contigs with bwa mem ({} kb)".format(len(contigs),
                                                            sum([len(i) for i in contigs.values()]) / 1000.))
        with open(all_contigs_file, "w") as out:
            for k, v in contigs.items():
                out.write("@{}\n{}\n+\n{}\n".format(k, v, "I"*len(v)))

        dodi_contigs(pfix, ref)# fai, l_match_score, l_mismatch_cost, l_gap_open, l_gap_extend)

    #
    all_mapped = collect_mappings(pfix, "bwa", contig_name=target,
                                  match=l_match_score, mis_match=l_mismatch_cost,
                                  gap_open=l_gap_open, gap_extend=l_gap_extend,
                                  min_identity=min_identity)  # Similar parameters as bwa mem

    print("Processing alignments:", len(all_mapped))

    #
    # Add a chromosome integer key column
    unique, keys = np.unique(all_mapped["tName"], return_inverse=True)
    all_mapped["chrom_key"] = keys

    #
    all_mapped = all_mapped.sort_values(["sample", "Q_name", "qStart"])
    v = np.array(all_mapped.Q_name)
    #sub = all_mapped.as_matrix(columns=["chrom_key", "qStart", "qEnd", "Last_score"]).astype(float)
    sub = all_mapped[["chrom_key", "qStart", "qEnd", "Last_score"]].values.astype(np.dtype("f"))

    #
    change_points = np.where(v[:-1] != v[1:])[0]
    # Need to add a reference to the last contig, after the last change point
    if len(change_points) == 0:
        change_points = [0]
    else:
        change_points = list(change_points) + [change_points[-1] + 1]

    contig_lengths = [len(contigs[all_mapped["contig_names"].iloc[index]]) for index in change_points]

    #
    count = 0
    all_t = []
    keep_indexes = []  # Grow the list of good alignment indexes
    t = time.time()

    for l_idx, array_idx in enumerate(slice_range(change_points, len(all_mapped))):

        length = contig_lengths[l_idx]
        data = sub[slice(*array_idx)]

        try:
            indexes, max_s = best_path.optimal_path(data, length,
                                                    inter_cost=inter_cost,
                                                    ins_cost=ins_cost,  # Score goes down quicker with an insertion
                                                    hom_cost=hom_cost,  # A match is counted twice for each overlap
                                                    intra_cost=intra_cost)
        except:
            print(change_points)
            print(array_idx)
            print(data)
            print(length)
            raise ValueError
        # Debug using python script:
        # indexes, max_s = find_path2.optimal_path(data, length, inter_cost=20, ins_cost=1, hom_cost=2.5, intra_cost=10)
        # Check if oython and cython versions match
        # assert ((indexes, max_s) == find_path2.optimal_path(data, length,
        #                                                     inter_cost=inter_cost,
        #                                                     ins_cost=ins_cost,
        #                                                     hom_cost=hom_cost,
        #                                                     intra_cost=intra_cost))
        # print max_s == max_s
        # print data
        # print find_path2.optimal_path(data, length,
        #                         inter_cost=inter_cost,
        #                         ins_cost=ins_cost,
        #                         hom_cost=hom_cost,
        #                         intra_cost=intra_cost)

        # Filter single alignments or alignments with low score
        if (len(indexes) > 1 or keep_single_alignments) and max_s > 0:

            # Add the array index start point to the returned indexes
            keep_indexes += [i + array_idx[0] for i in indexes]

    print("Time to find paths:", time.time() - t)

    if not target:
        print("Saving results")
        dfs = all_mapped.iloc[keep_indexes]

        # Drop contigs with low overall identity, these may be repetitive, poorly described contigs
        # Drop contigs if they have only one alignment to 'normal' chromosome (chr1, chr2 etc, not chr_un..)
        t = []
        keepers = set([])
        for cont, d in dfs.groupby("contig_names"):
            idm = d.Identity.mean()
            t.append(idm)
            passed = False
            if idm > 0.95:
                normal_chrom_counts = Counter([len(i.replace("chr", "")) < 3 for i in list(d["tName"])])
                if len(d) == 1:
                    continue  # SKip single alignments

                if len(d) <= 2:
                    if 1 in normal_chrom_counts:
                        passed = True
                else:
                    if 1 in normal_chrom_counts and normal_chrom_counts[1] > 1:
                        passed = True
            if passed:
                keepers.add(cont)

        print("Fraction of contigs with high enough identity:", len([i for i in t if i > 0.95]) / float(len(t)))
        dfs = dfs[dfs["contig_names"].isin(keepers)]

        print("Kept alignments", len(dfs), "kept contigs", len(keepers))

        # import matplotlib.pyplot as plt
        # plt.figure()
        # plt.hist(t, bins=25)
        # plt.show()
        print("Saving as: ", "Best_alignments.{}.csv".format(pfix))
        dfs.to_csv("Best_alignments.{}.csv".format(pfix))
    # print dfs

    else:
        d = all_mapped.iloc[keep_indexes] #[["T_name", "T_start", "T_end"]]
        print(d)
