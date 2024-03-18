"""
Given a target region recursively get the discordant reads which link to this region, then assemble.
"""

import pysam
from collections import defaultdict, deque, Counter
import sys
import glob
import intervaltree
from subprocess import call, Popen, PIPE
import pandas as pd
import random
import os
import itertools
from collections import Counter
import time
from rich.progress import Progress, TextColumn, BarColumn, TaskProgressColumn, TimeRemainingColumn, TimeElapsedColumn
import recursive_find
from pathlib import Path
random.seed(1234)


def load_svs(data, pad, pfix=""):

    data["chr1_start"] = data["start1"] - pad
    data["chr1_end"] = data["end1"] + pad
    data["chr2_start"] = data["start2"] - pad
    data["chr2_end"] = data["end2"] + pad

    return data


def run_kmer(singles, pairs):
    print("running kmer to downsample data")
    # call("cutadapt -q 15,15 -o {s}.trim.fastq {s}".format(s=singles), shell=True)
    call("normalize-by-median.py -k 21 -C 20 -x 5e8 {s} -o {s}.keep.fastq".format(s=singles), shell=True)


def load_vcfs(chain_file):
    all_breaks = defaultdict(list)
    breaks_list = defaultdict(lambda: defaultdict(list))  # Sample: Chromosome: list of breakpoints
    count = 0
    vcfs = glob.glob(chain_file)
    if len(vcfs) > 0 :
        print("vcfs: ", vcfs)
        # type_thresholds = {"INS": 0.2, "DEL": 0.3, "INV": 0.15, "DUP": 0.15, "TRA": 0.35}
        for v in vcfs:
            #print(v)
            f = pysam.VariantFile(v)
            samp = os.path.basename(v).replace(".vcf", "")
            for line in f.fetch():
                #if list(line.filter.keys())[0] != "lowProb" or float(line.samples[s]["PROB"]) > 0.2 or float(line.samples[s]["SU"]) > 4: 
                # if float(line.samples[s]["PROB"]) >= type_thresholds[line.info["SVTYPE"]]:
                if (line.info["SVTYPE"] != "TRA" and line.stop-line.pos >= 50000) or line.info["SVTYPE"] == "TRA":
                    chr1 = line.contig
                    pos1 = line.pos
                    end1 = line.stop
                    chr2 = line.info["CHR2"]
                    if line.info["SVTYPE"] == "TRA":
                        pos2 = line.info["CHR2_POS"]
                        end2 = line.info["CHR2_POS"]
                    else:
                        pos2 = line.pos
                        end2 = line.stop
                    samp = v.split("/")[-1].split(".")[0]
                    breaks_list[samp][chr1].append(pos1)
                    breaks_list[samp][chr2].append(pos2)
                    all_breaks[samp].append([chr1, pos1, end1, chr2, pos2, end2])
                    count += 1
    return all_breaks


def load_beds(chain_file):
    all_breaks = {}
    breaks_list = defaultdict(lambda: defaultdict(list))  # Sample: Chromosome: list of breakpoints
    count = 0
    beds = glob.glob(chain_file)
    if len(beds) > 0:
        print("beds: ", beds)
        for v in beds:
            df = pd.read_csv(v, sep="\t", index_col=None)
            for idx, r in df.iterrows():
                chr1 = r["chrom1"]
                pos1 = r["start1"]
                end1 = r["end1"]
                chr2 = r["chrom2"]
                pos2 = r["start2"]
                end2 = r["end2"]
                samp = r["Sample"]
                breaks_list[samp][chr1].append(pos1)
                breaks_list[samp][chr2].append(pos2)
                all_breaks[samp].append([chr1, pos1, end1, chr2, pos2, end2])
                count += 1
    return all_breaks



#if __name__ == "__main__":
def main(bams_glob, chain_file, ref, pad, pfix="chained_and_unchained", OUT="regions_all", outfix="all"):
    # Paths to bam files
    bam_paths = {i.split("/")[-1].split(".")[0]: i for i in glob.glob(bams_glob)} # "/mnt/breast/*.*am"
    print(bam_paths)

    # Chained and unchained # "../chain_link/chain_out_50/found_chains_p0.4_t20/all_svs.unique.chains.csv"
    if chain_file.endswith(".csv"):
        data = pd.read_csv(chain_file, index_col=None)
    elif chain_file.endswith(".vcf"):
        all_breaks = load_vcfs(chain_file)
        for s in all_breaks.keys():
            #print(pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "chrom2", "start2"]))
            df = pd.concat([df, pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2"])])
        data = df
        print(df)
    elif chain_file.endswith(".bed"):
        all_breaks = load_vcfs(chain_file)
        for s in all_breaks.keys():
            #print(pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "chrom2", "start2"]))
            df = pd.concat([df, pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2"])])
        data = df
        print(df)
    elif os.path.isdir(chain_file):
        all_vcfs = load_vcfs(f"{chain_file}/*.vcf")
        all_beds = load_beds(f"{chain_file}/*.bed")
        all_breaks = {**all_vcfs, **all_beds}
        df = pd.DataFrame(columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2"])
        for s in all_breaks.keys():
            #print(pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "chrom2", "start2"]))
            df = pd.concat([df, pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2"])])
        data = df
        print(df)


    # pfix = "chained_and_unchained"
    # OUT = "regions_all"
    # outfix = "all"

    print("=" * 40)
    print(pfix, OUT)
    print("=" * 40)
    data = load_svs(data, pad=650, pfix=pfix)  # Pad = read length + insert length
    samples = [Path(i).stem for i in bam_paths]
    data = data[data["Sample"].isin(samples)]

    # Group data by sample and recursively find reads for each breaksite
    grps = data.groupby("Sample")

    for sample, df in grps:

        print("-" * 40)
        print(sample)

        if not os.path.exists("./{o}/{s}".format(o=OUT, s=sample)):
            os.mkdir("./{o}/{s}".format(o=OUT, s=sample))
        else:
            continue # maybe add prompt that can skip or overwrite?

        # Collect reads from the .bam file
        if True:
            sample_bam = bam_paths[sample]
            bam = pysam.AlignmentFile(sample_bam, "rb")

            out_singles = pysam.AlignmentFile("./{o}/{s}/{s}_singles.bam".format(o=OUT, s=sample), "wb", template=bam)
            out_pairs = pysam.AlignmentFile("./{o}/{s}/{s}_paired.bam".format(o=OUT, s=sample), "wb", template=bam)

            with Progress(TextColumn("[progress.description]{task.description}"), BarColumn(), TaskProgressColumn(), TimeRemainingColumn(), TimeElapsedColumn()) as progress: 
                task = progress.add_task("[red]Recursive serach...", total=len(df)*2)
                discordants = []
                for idx, r in df.iterrows():
                    d1 = recursive_find.recurive_find(bam, r.chrom1, r.chr1_start, r.chr1_end, pad, max_depth=3)
                    progress.update(task, advance=1)
                    d2 = recursive_find.recurive_find(bam, r.chrom2, r.chr2_start, r.chr2_end, pad, max_depth=3)
                    discordants += recursive_find.filter_duplicates(d1+d2)
                    progress.update(task, advance=1)

            unique = recursive_find.filter_duplicates(discordants)
            print("Reads to assemble", len(unique))

            # Sort reads by Q name and split into singles and paired reads for best assembly
            singles, pairs = recursive_find.sort_disc(unique)

            # Write a fastQ file
            with out_singles, out_pairs:
                for read in unique:   # WRITE OUT ALL READS, NOT JUST SINGLES
                    out_singles.write(read)
                for read in pairs:
                    out_pairs.write(read)

            call("samtools sort -o ./{o}/{s}/{s}_singles.srt.bam ./{o}/{s}/{s}_singles.bam ; \
            samtools index ./{o}/{s}/{s}_singles.srt.bam".format(o=OUT, s=sample), shell=True)
            call("samtools sort -o ./{o}/{s}/{s}_paired.srt.bam ./{o}/{s}/{s}_paired.bam ; \
            samtools index ./{o}/{s}/{s}_paired.srt.bam".format(o=OUT, s=sample), shell=True)

            # Use bedtools to convert to fastq for singles and pairs
            sr = "./{o}/{s}/{s}_singles.fastq".format(o=OUT, s=sample)
            call("bedtools bamtofastq -i {} -fq {}".format("./{o}/{s}/{s}_singles.bam".format(o=OUT, s=sample),
                                                       sr), shell=True)

            pr = "./{o}/{s}/{s}_paired.fastq".format(o=OUT, s=sample)
            call("bedtools bamtofastq -i {} -fq {}".format("./{o}/{s}/{s}_paired.bam".format(o=OUT, s=sample),
                                                           pr), shell=True)

            # Use kmer to downsample reads for assembly
            run_kmer(sr, pr)

        # Assemble
        if True:
            # Assemble using SPAdes
            # Using too long a kmer can sometimes fail, if this is the case, try with a shorter set of kmers:
            kmers = [21, 33, 55, 77, 111][::-1]
            for k in [",".join(map(str, kmers[:j])) for j in range(len(kmers), 0, -1)]:
                print("Assembling with", k)
                s = "./{o}/{s}/{s}_singles.fastq".format(o=OUT, s=sample)
                pe = "./{o}/{s}/{s}_paired.fastq".format(o=OUT,s=sample)
                base = "./{o}/{s}/{s}".format(o=OUT, s=sample)
                                                            # 1          2      3      4   5
                out = Popen("bash ./assemble_and_map.sh {s}.keep.fastq {base} {pfix} {k} {pe} ${ref}".format(
                    s=s,
                    pe=pe,
                    pfix=pfix,
                    k=k,
                    base=base),
                    shell=True, stderr=PIPE, stdout=PIPE)
                c = out.communicate()
                #print("c1")
                #print(str(c[1]))
                #print("c0")
                #print(str(c[0]))
                if "======= SPAdes pipeline finished abnormally and WITH WARNINGS!" in str(c[0]):
                    print(k, "======= SPAdes pipeline finished abnormally and WITH WARNINGS!")
                else:
                    break

            else:
                raise StopIteration("ASSEMBLY FAILED: {}".format(sample))

            # Save contigs as well
            call("cp ./assembly/{p}/contigs.fasta ./{o}/{s}/".format(p=pfix, o=OUT, s=sample), shell=True)

        # Run bwa mem and samtools for gw viewing
        if False:
            call("bwa mem -a ~/Documents/Data/db/hg38/hg38.fa \
            regions_{ofix}/{samp}/contigs.fasta > regions_{ofix}/{samp}/{samp}.sam".format(ofix=outfix, samp=sample),
                 shell=True)

            call("samtools view -bF4 regions_{ofix}/{samp}/{samp}.sam \
            | samtools sort - regions_{ofix}/{samp}/{samp}".format(ofix=outfix, samp=sample, shell=True))

            call("samtools index regions_{ofix}/{samp}/{samp}.bam".format(ofix=outfix, samp=sample), shell=True)

