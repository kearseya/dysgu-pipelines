"""
This is run after the remove_common_sites script and takes the "All_mapped_contigs.unique.csv" file as input.
After this run the plot_sequence.py or to add the chain label annotation run chain_finder scripts,
the finder_contigs script will add the chained annotation.
"""

import glob
import numpy as np
from collections import defaultdict
import pandas as pd
import pickle
import sys
import random
from skbio.alignment import local_pairwise_align_ssw, StripedSmithWaterman
from skbio.sequence import DNA
from operator import mul
from math import pow
from multiprocessing import Pool
import time


def load_contigs(input_contigs, min_cont_length=150):

    contigs = defaultdict(list)

    lines = open(input_contigs, "r").readlines()
    cont = None
    for l in lines:
        if l[0] == ">":
            cont = l.strip()[1:]
        else:
            contigs[cont].append(l.strip())

    return {k: "".join(v) for k, v in contigs.items() if sum(len(i) for i in v) > min_cont_length}


def load_maf(f):

    maf = open(f, "r")
    records = []
    # fields = ["score", "E", "chrom", "ref_pos", "ref_end", "ref_seq", "cont_name", "cont_pos", "cont_end",
    #           "cont_strand", "cont_seq"]
    r = {"Maf_indexes": []}  # Fill this out and then update for each new record
    count = 0
    maf_header = ""
    for index, l in enumerate(maf):
        if l[0] == "#":
            maf_header += l
            continue

        a = [i for i in l.strip().split(" ") if i != '']
        r["Maf_indexes"].append(index)  # Keep track of the line number for filtering later

        if len(a) == 0:
            records.append(r)
            r = {"Maf_indexes": []}
            count = 0
            continue

        if a[0] == "a":
            r["score"] = int(a[1].split("=")[1])  # Keep as strings for now
            r["E"] = a[3].split("=")[1]

        if a[0] == "s" and count == 0:
            # The reference sequence
            r["chrom"] = a[1]
            r["ref_pos"] = a[2]
            r["ref_end"] = int(a[2]) + int(a[3])
            r["ref_seq"] = a[6]
            count += 1

        elif a[0] == "s" and count == 1:
            sample, name = a[1].split("_", 1)
            r["cont_name"] = a[1]
            r["sample"] = sample
            r["cont_pos"] = a[2]
            r["cont_end"] = int(a[2]) + int(a[3])
            r["cont_strand"] = a[4]
            r["cont_seq"] = a[6]
            r["cont_length"] = a[5]

    data = pd.DataFrame.from_records(records).to_csv("out_data/alignments.all.csv")

    # Write a header to use
    with open("out_data/maf_header.txt", "w") as h:
        h.write(maf_header)


def match_alignments(maf, psl):

    # Filter out other chromosomes
    chroms = set(psl.T_name)
    maf = maf[maf["chrom"].isin(chroms)]

    start = set(psl.T_start)
    end = set(psl.T_end)
    Qstart = set(psl.Q_start)

    # The maf alignments can be the reverse complement of the psl which is confusing. Therefore the start and ends
    # may be flipped around in their meaning.
    # Filter maf alignments which dont have a matching start or end coordinate
    # If the maf alignment is on the reverse strand the contig position is counted from the end!
    abs_pos = []
    for i, r in maf.iterrows():
        if r["cont_strand"] == "+":
            abs_pos.append(r["cont_pos"])
        else:
            abs_pos.append(r["cont_length"] - r["cont_end"])

    maf["Q_start"] = abs_pos

    maf = maf[maf["Q_start"].isin(Qstart)]

    maf = maf[(maf["ref_pos"].isin(start) & maf["ref_end"].isin(end))
            | (maf["ref_pos"].isin(end) & maf["ref_end"].isin(start))]

    # Keep a subset of the maf file
    maf = maf.sort_values("Q_start")[["score", "E", "ref_pos", "ref_end", "ref_seq", "cont_seq", "Maf_indexes"]].reset_index()

    assert all(i == j for i, j in zip(list(maf.ref_pos), list(psl.T_start)))

    return maf


def rev_comp(s):
    d = {"A": "T", "C": "G", "T": "A", "G": "C", "N": "N"}
    return "".join([d[i] for i in s][::-1])


def match_segments(a, b, start_first, start_second, cont, m):
    # +2 match -3 for mismatch
    aln = local_pairwise_align_ssw(DNA(a), DNA(b), gap_open_penalty=10, gap_extend_penalty=10)

    # Now everything is on the forward strand, normal indexing can be used instead of counting from the end
    ap, bp = aln[2]
    score = aln[1]

    startA, endA = int(ap[0] + start_first), int(ap[1] + start_first + 1)
    startB, endB = int(bp[0] + start_second), int(bp[1] + start_second + 1)

    sequenceA = cont[startA: endA]  # Use the alignment rather than the contig. keeps gaps.
    sequenceB = cont[startB: endB]

    if str(sequenceA) != str(aln[0][0]).replace("-", "") or str(sequenceB) != str(aln[0][1]).replace("-", ""):
        raise AssertionError(m, sequenceA, sequenceB, aln)

    sequenceA = str(aln[0][0])  # Use the alignment rather than the contig. keep gaps.
    sequenceB = str(aln[0][1])

    # Pval cal derived from
    # https://stats.stackexchange.com/questions/26988/probability-of-finding-a-particular-sequence-of-base-pairs
    prob_match = 0.25
    prob_mismatch = 0.75

    seq1 = aln[0][0]
    seq2 = aln[0][1]

    # Calculate the probability of a sequence match
    prob = reduce(mul, (prob_match if i == j else prob_mismatch for i, j in zip(seq1, seq2)))
    # Put the probability of Not-a-match to the number of possible seqeunce positions
    seq_pos = max((len(a), len(b)))  # Choose the londer of the two sequences to calculate the probability

    match_len = len(aln[0][0])
    p = 1 - pow(1 - prob, seq_pos - match_len + 1)

    # Use random sequences to calculate a probability distribution, then use this to calculate the probability of
    # obtaining an equvalent or higher alignment score

    return {"score_ssw": score, "A_start": startA, "A_end": endA, "B_start": startB, "B_end": endB, "pval": p,
                        "sequenceA": sequenceA, "sequenceB": sequenceB, "LengthA": len(a), "LengthB": len(b)}


def add_simularity(m, cont):
    """
    Adjacent alignments in the contig are compared for simularity.
    Sections of microhomology are always included in the comparison.
    For insertions, the insertion is only added to one of the adjacent sequences. Insertions from which fall before the
    first alignment or after the second alignment are ignored. The decision of which alignment the insertion is joined
    to is governed by the highest simularity score.
    :param m: The contig dataframe
    :param cont: The contig sequence
    :return: The contig dataframe with appended simularity information
    """
    # Add simulairy in the forward and reverse directions, there might be a clear winner and direction to events.
    ds = [i[1] for i in m.iterrows()]
    records = []

    for i in range(len(ds)-1):
        first = ds[i]
        second = ds[i+1]

        # Isolate contig sequence
        # Calculate insertion length
        ins_length = ds[i+1]["Q_start"] - ds[i]["Q_end"]

        comparisons = []  # Make comparisons between adjacent sequences, choose the best comparison from here

        a = cont[ds[i]["Q_start"]: ds[i]["Q_end"]]
        b = cont[ds[i + 1]["Q_start"]: ds[i + 1]["Q_end"]]

        start_first = ds[i]["Q_start"]
        start_second = ds[i + 1]["Q_start"]

        if ins_length > 0:
            a_ins = cont[ds[i]["Q_start"]: ds[i]["Q_end"] + ins_length]  # Insertion inlcuded with segment A
            b_ins = cont[ds[i+1]["Q_start"] - ins_length: ds[i+1]["Q_end"]]  # Insertion included with segment B

            comparisons.append( [(a_ins, b, start_first, start_second, cont, m),
                                (a, b_ins, start_first, start_second - ins_length, cont, m)] )
        else:
            comparisons.append([(a, b, start_first, start_second, cont, m)])

        #
        for comp in comparisons:

            if len(comp) == 1:
                r = match_segments(*comp[0])  # No insertions
                r["choice"] = "no_ins"

            else:  # Insertion, choose which side the insertion lies
                r1 = match_segments(*comp[0])
                r2 = match_segments(*comp[1])

                # Choose the one with the minimal pval or the first one if identical
                if r2["pval"] < r1["pval"]:
                    r = r2
                    r["choice"] = "ins_add_right"
                else:
                    r = r1
                    r["choice"] = "ins_add_left"

            records.append(r)

    if len(records) == 0:
        return []

    records.append({k: None for k in records[-1]})  # Add an empty record to keep indexes same length
    return pd.concat([m, pd.DataFrame(records)], axis=1)


def get_psl(f, pfix):
    samps = pd.read_csv(f)
    d = [i for i in samps.columns if "pct" in str(i) or str(i) in set(["A", "C", "G", "T"]) or "rho" in str(i)]
    samps = samps[[i for i in samps.columns if i not in d]]

    print(map(str, samps.columns))

    # This is for use to plot the resulting SVs as curves on a chromosome, as for the delly calls.
    # See the plot_filtered_unique_from_vcf.py script
    # Sample: svtype: chromosome: [(start, end) ... ]

    dvcf = defaultdict(lambda: defaultdict(lambda: defaultdict(list)))#{t: defaultdict(list) for t in ["DEL", "INV", "TRA", "DUP"]})}

    # for key in ["DB" + str(i) for i in range(37, 59)] + ["DB" + str(i) for i in range(67, 117)]:
    #     dvcf[key] = {t: defaultdict(list) for t in ["DEL", "INV", "TRA", "DUP"]}

    # Load contigs from contigs file
    contigs = load_contigs("out_data/all_contigs.all.fa")

    # Load alignments from maf file
    load_maf("last_output.all.maf")
    maf = pd.read_csv("out_data/alignments.all.csv")
    maf = {k: df for k, df in maf.groupby(["sample", "cont_name"])}

    dfs_to_concat = []
    records = []
    for key, df in samps.groupby("contig_names"):

        # if key != "NODE_1_length_3106_cov_11.1813":  # For debuging
        #     continue
        print(key)

        alns_maf = maf[(df["sample"].iloc[0], key)]
        matched = match_alignments(alns_maf, df)

        # Join the psl and maf alignments
        df = pd.concat([df.reset_index(), matched], axis=1)

        # Add measure of simularity between adjacent segments
        df = add_simularity(df, contigs[key]) #"{}_{}".format(df["sample"].iloc[0], key)])

        if len(df) == 0:
            # Only a single alignment in the contig, no SVs
            continue
        dfs_to_concat.append(df)

        #
        if len(df) > 2:
            complex = True
        else:
            complex = False
        rows = [d.to_dict() for i, d in df.iterrows()]
        for idx in range(len(rows)-1):

            seg1 = rows[idx]
            seg2 = rows[idx+1]

            microhomology = 0
            insertion = 0
            if seg1["Q_end"] > seg2["Q_start"]:
                microhomology = seg1["Q_end"] - seg2["Q_start"]

            elif seg1["Q_end"] < seg2["Q_start"]:
                insertion = seg2["Q_start"] - seg1["Q_end"]

            chr1, chr2 = seg1["T_name"], seg2["T_name"]
            pos1, pos2 = seg1["T_end"], seg2["T_start"]

            strand1, strand2 = seg1["strand"], seg2["strand"]

            # Need to check join type is correct
            if chr1 != chr2:
                svtype = "TRA"
                if strand1 == strand2:
                    if strand1 == "+":
                        if pos1 < pos2:
                            jointype = "3to5"
                        else:
                            jointype = "5to3"
                    else:
                        if pos1 < pos2:

                            jointype = "5to3"
                        else:
                            jointype = "3to5"
                else:
                    if strand1 == "+":
                        jointype = "3to3"
                    else:
                        jointype = "5to5"

            else:
                if strand1 == strand2:
                    if strand1 == "+":
                        if pos1 < pos2:
                            jointype = "3to5"
                        else:
                            svtype = "DUP"
                            jointype = "5to3"
                    else:
                        if pos1 < pos2:
                            svtype = "DUP"
                            jointype = "5to3"
                        else:
                            svtype = "DEL"
                            jointype = "3to5"

                else:
                    svtype = "INV"
                    if strand1 == "+":
                        jointype = "3to3"
                    else:
                        jointype = "5to5"

            # Use the sequence information to find the insertion or microhomology sequences
            contig = contigs["{}_{}".format(seg1["sample"], seg1["Q_name"])]

            micro_h_seq = ""
            try:
                if microhomology > 0:
                    micro_h_seq = contig[int(seg1["Q_end"]) - int(microhomology): int(seg1["Q_end"])]
            except:
                print(seg1, seg2)
                print(contig)
                print("hom")
                continue
            insertion_seq = ""
            try:
                if insertion > 0:
                    insertion_seq = contig[int(seg1["Q_end"]): int(seg2["Q_start"])+1]
            except:
                print(seg1)
                print(seg2)
                print(contig)
                print("insertion")
                continue

            records.append({"sample": seg1["sample"],
                            "read": seg1["Q_name"],
                            "chrom1": chr1, "chrom2": chr2,
                            "pos1": pos1, "pos2": pos2,
                            "strand1": strand1, "strand2": strand2,
                            "jointype": jointype, "svtype": svtype,
                            "microhomology": microhomology, "insertion": insertion,
                            # Add sequence information
                            "micro_h_seq": micro_h_seq, "insertion_seq": insertion_seq,
                            "seg1_E": seg1["E"], "seg2_E": seg2["E"], "Complex": complex})

            samp = seg1["sample"]
            if chr1 == chr2:
                dvcf[samp][svtype][chr1].append((pos1, pos2))
            else:
                dvcf[samp][svtype][chr1].append((pos1, pos1 + 1))  # Translocation will appear as lines
                dvcf[samp][svtype][chr2].append((pos2, pos2 + 1))


    # Save the SV data
    pd.DataFrame.from_records(records).to_csv("out_data/alignments.all.SVdata.csv")

    # Save the dvcf for plotting elsewhere
    #pickle.dump(dvcf, open("dvcf.pkl", "wb"))

    drop = set([u"index", u"Unnamed: 0", u"to_split", u"key", u"level_0"])
    df = pd.concat(dfs_to_concat).reset_index()
    df = df[[i for i in df.columns if i not in drop]]

    # Save the updated psl alignments
    df.to_csv("out_data/alignments.all.updated_psl.csv")
    return df


def worker(lengths):
    l1 = lengths[0]
    l2 = lengths[1]
    alphabet = ["A", "C", "T", "G"]
    samples = 1000

    scores = [lengths]

    rand1 = "".join([random.choice(alphabet) for i in range(l1)])

    for i in range(samples):
        rand2 = "".join([random.choice(alphabet) for j in range(l2)])
        aln = StripedSmithWaterman(rand1, gap_open_penalty=10, gap_extend_penalty=10, score_only=True)(
            rand2)
        scores.append(aln.optimal_alignment_score)
    return scores


def make_aln_score_table(lengths):
    # Make a list of alignment scores for random sequences of length1 and length2. Use multiprocessing for speed.
    table = {}
    pool = Pool(processes=8)
    args = [map(int, i) for i in lengths]
    res = pool.map(worker, args)
    for item in res:
        table[tuple(item[0])] = item[1:]

    pickle.dump(table, open("prob_table.pkl", "wb"))


def add_alignment_probability():
    df = pd.read_csv("out_data/alignments.all.updated_psl.csv")

    df1 = df.dropna()  # Must remove na values
    seq_lengths = set([tuple(sorted(i)) for i in zip(df1.LengthA, df1.LengthB)])

    print("Making table for n comparisons: ", len(seq_lengths))
    make_aln_score_table(seq_lengths)
    prob_table = pickle.load(open("prob_table.pkl", "rb"))

    probs = []
    for idx, r in df.iterrows():

        try:
            key = tuple(sorted([int(r.LengthA), int(r.LengthB)]))
        except:
            # The last alignment in a contig will have lengths (nan, nan) due to pairwise comparison leaving out the
            # last alignment.
            probs.append(None)
            continue

        dist = prob_table[key]
        score = r.score_ssw
        prob = 1 - (len([i for i in dist if score >= i]) / float(len(dist)))
        probs.append(prob)

    assert (len(probs) == len(df))
    df["Pairwise_aln_prob"] = probs
    print("Fraction of significant pairwise alignments:", len([i for i in probs if i!=None and i < 0.01]) / float(len(probs)))
    df.to_csv("out_data/alignments.all.updated_psl.csv")


if __name__ == "__main__":
    get_psl("All_mapped_contigs.unique.csv", "all")

    # Now calculate an alignment score probability for all the A and B adjacent comparisons
    add_alignment_probability()

