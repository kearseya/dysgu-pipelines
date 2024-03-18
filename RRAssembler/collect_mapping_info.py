import pysam 
import pandas as pd
import sys
from collections import defaultdict
import os
import glob
from intervaltree import Interval, IntervalTree
#import pkg_resources

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



def get_query_pos_from_cigartuples(r):
    # Infer the position on the query sequence of the alignment using cigar string
    start = 0
    query_length = r.infer_read_length()  # Note, this also counts hard-clips
    end = query_length
    if r.cigartuples[0][0] == 4 or r.cigartuples[0][0] == 5:
        start += r.cigartuples[0][1]
    if r.cigartuples[-1][0] == 4 or r.cigartuples[-1][0] == 5:
        end -= r.cigartuples[-1][1]
    return start, end, query_length

def flip_not_reverse(s):
    rev = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join([rev[i] for i in [*s]])

def reverse_comp(s):
    rev = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join([rev[i] for i in [*s][::-1]])


def mapping_info(f, outf, chain_file, pad):
    # load vcf files for breakpoints
    vcf_breaks = {}
    if chain_file.endswith(".csv"):
        df = pd.read_csv(chain_file, index_col=None)
    elif chain_file.endswith(".vcf"):
        all_breaks = load_vcfs(chain_file)
        for s in all_breaks.keys():
            #print(pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "chrom2", "start2"]))
            df = pd.concat([df, pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2"])])
            vcf_breaks[s] = {}
        print(df)
    elif chain_file.endswith(".bed"):
        all_breaks = load_vcfs(chain_file)
        for s in all_breaks.keys():
            #print(pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "chrom2", "start2"]))
            df = pd.concat([df, pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2"])])
            vcf_breaks[s] = {}
        print(df)
    elif os.path.isdir(chain_file):
        all_vcfs = load_vcfs(f"{chain_file}/*.vcf")
        all_beds = load_beds(f"{chain_file}/*.bed")
        all_breaks = {**all_vcfs, **all_beds}
        df = pd.DataFrame(columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2"])
        for s in all_breaks.keys():
            #print(pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "chrom2", "start2"]))
            df = pd.concat([df, pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2"])])
            vcf_breaks[s] = {}
        print(df)

    for st in df.groupby("Sample"):
        s = st[1]
        sample = st[0]
        chroms = set(s["chrom1"].tolist())
        for c in chroms:
            chrom_tuples = list(s[s["chrom1"] == c][["start1", "end1"]].itertuples(index=False, name=None)) + list(s[s["chrom2"] == c][["start2", "end2"]].itertuples(index=False, name=None))
            #print(chrom_tuples)
            #print(IntervalTree(Interval(b-100, e+100) for b, e in chrom_tuples[0:10]))
            vcf_breaks[sample][c] = IntervalTree(Interval(begin-pad, end+pad) for begin, end in chrom_tuples)
    

    af = pysam.AlignmentFile(f, 'r')
    d = defaultdict(list)
    for a in af.fetch(until_eof=True):
        if not a.flag & 4:
            d[a.qname].append(a)

    res = []
    no = 0
    yes = 0
    for qname, v in d.items():
        sample = qname.split("_")[0]
        flag = [(index, i) for index, i in enumerate(v) if not i.flag & 2304]
        if len(flag) > 1:  # todo check bug in dodi, not currently setting primary alignment flag properly
            flag = [flag[flag.index(max(flag, key=lambda x: x[1].get_tag('AS')))]]

        if len(flag) != 1:
            print('Error in ', f, 'flag problem', len(flag), [i.flag for i in v])
            quit()
        pri_index, pri_read = flag[0]
        primary_reverse = bool(pri_read.flag & 16)
        if primary_reverse == False:
            seq = pri_read.get_forward_sequence()
        else:
            seq = flip_not_reverse(pri_read.query_sequence)
        n_aligns = len(v)
        any_seq = False
        for index, a in enumerate(v):
            qstart, qend, qlen = get_query_pos_from_cigartuples(a)
            align_reverse = bool(a.flag & 16)
            if primary_reverse != align_reverse:
                start_temp = qlen - qend
                qend = start_temp + qend - qstart
                qstart = start_temp
            pri = index == pri_index
            if not pri:
                no += 1
            else:
                yes += 1
                any_seq = len(seq) if seq else 0
            
            current_chrom = af.get_reference_name(a.rname)
            vcf_val = False
            if vcf_breaks[sample].get(current_chrom) != None:
                vcf_val = vcf_breaks[sample][current_chrom].overlaps(a.reference_start+1, a.reference_end)
            print(vcf_val)

            res.append(
                {'sample': sample,
                 'qname': a.qname,
                 'n_alignments': n_aligns,
                 'chrom': current_chrom,
                 'rstart': a.reference_start + 1,
                 'rend': a.reference_end,
                 'strand': '-' if align_reverse else '+',
                 'qstart': qstart,
                 'qend': qend,
                 'qlen': qlen,
                 'aln_size': qend - qstart,
                 'mapq': a.mapq,
                 'alignment_score': a.get_tag('AS'),
                 'seq': seq if pri else '',
                 'cigar': a.cigarstring,
                 'vcf': vcf_val
                }
            )

        if not any_seq:
            print('missing', qname, [(len(vv.seq), vv.infer_query_length()) if vv.seq else vv.infer_query_length() for vv in v])
            quit()

    df = pd.DataFrame.from_records(res).sort_values(['qname', 'qstart'])

    bad_anchors = []
    # flag reads with small anchoring alignments
    for grp, d in df.groupby('qname'):
        aln_s = list(d['aln_size'])
        if aln_s[0] < 50 or aln_s[-1] < 50:
            bad_anchors += [1] * len(d)
        else:
            bad_anchors += [0] * len(d)
    df['short_anchor<50bp'] = bad_anchors

    df = df.sort_values(['n_alignments', 'qname', 'qstart'], ascending=[False, True, True])

    df = df[['chrom', 'rstart', 'rend', 'qname', 'n_alignments', 'aln_size', 'qstart', 'qend', 'strand', 'mapq', 'qlen', 'alignment_score', 'short_anchor<50bp', 'seq', 'cigar', 'vcf']]
    df.to_csv(outf, index=False, sep="\t")


# if __name__ == '__main__':
#     import argparse
#     parse = argparse.ArgumentParser()
#     parse.add_argument('--bam', help='bam file to assess')
#     parse.add_argument('--out', help='out put bed file')
#     args = parse.parse_args()
#     mapping_info(args.bam, args.out)
#     print('Done')
