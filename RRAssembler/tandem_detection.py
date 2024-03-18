import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from pyfaidx import Fasta
from dysgu_repeats import compute_rep
from skbio.alignment import StripedSmithWaterman
from plot_contigs import plot_individual, set_colors

def base_perc(seq):
    bases = {"A": 0, "T": 0, "C": 0, "G": 0}
    seq = seq.upper()
    for b in bases.keys():
        bases[b] = seq.count(b)
    for b in bases.keys():
        bases[b] = bases[b]/len(seq)
    bases = {k: v for k, v in sorted(bases.items(), key=lambda item: item[1], reverse=True)}
    bpv = sum(list(bases.values())[:2])
    return bases, bpv

def simple_check(repeat, seq):
    start = seq[:4]
    end = seq[-4:]
    #print("start: ", start, " end: ", end)
    if start == repeat*2 or start == repeat[::-1]*2:
        #print("simple start")
        return True
    elif end == repeat*2 or end == repeat[::-1]*2:
        #print("simple end")
        return True
    else:
        #print("failed")
        return False

def print_details(df, seq=None, gap_seq=None, mapq=None, bp=None, bpv=None):
    print("--------------------")
    if seq != None:
        print(list(df["qname"])[0], len(seq))
        print(df)
        print("seq: ", seq)
        print("raw rep: ", compute_rep(seq))
    if gap_seq != None:
        print("gap: ", gap_seq)
        print("len: ", len(gap_seq))
        print("rep: ", compute_rep(gap_seq))
    if mapq != None:
        print("map: ", mapq)
    if bp != None:
        print("rbp: ", bp)
    if bpv != None:
        print("bpv: ", bpv)


def tandem_finder(df, qname, cd):
    ddf = df.copy()
    dropped = False
    ### FLANKING
    diff = {}
    ## if can have tandem sandwhiched
    if len(df) > 2:
        ## for chromosomes, check if close sections
        for d in df.groupby("chrom"):
            if len(d[1]) > 1:
                seq = list(df["seq"].dropna())[0]
                se = list(d[1][["rstart", "rend"]].itertuples())
                diff[d[0]] = []
                ## calc gap between chrom sections
                for pos, i in enumerate(se[:-1]):
                    diff[d[0]].append((se[pos+1][2]-i[1], se[pos+1][0]-i[0], i[0], se[pos+1][0]))
                for x in diff[d[0]]:
                    ## if gap in adjacent chrom sections
                    if -1000 < x[0] < 1000 and x[1] > 1:
                        gap_seq = seq[df.loc[x[2], "qend"]:df.loc[x[3], "qstart"]]
                        flank_seq = [seq[df.loc[x[2], "qstart"]:df.loc[x[2], "qend"]], seq[df.loc[x[3], "qstart"]:df.loc[x[3], "qend"]]] 
                        bp, bpv = base_perc(gap_seq)
                        if compute_rep(gap_seq) > 0.2:
                            # ddf = ddf.drop([i in range(x[2]+1, x[3])], axis=0)
                            sus_index = [i for i in range(x[2]+1, x[3])]
                            ## remove probable tandem repeat
                            for i in sus_index:
                                if df.loc[i, "mapq"] < 30:
                                    df = df.drop([i], axis=0)
                                    #print_details(df=df, seq=seq, gap_seq=gap_seq, mapq=list(df["mapq"]), bp=bp, bpv=bpv)
                            ## if contig now two close alignments (repeat in middle) remove
                            if len(df) == 2:
                                #ddf = ddf.drop([i[0] for i in se], axis=0)
                                return None

    ### ADJACENT
    if len(df) > 1 and dropped == False:
        mapq = list(df["mapq"])
        ## suspect low mapping quality of being tandem repeat
        if 0 in mapq:
            seq = list(df["seq"].dropna())[0]
            ## remove contigs that are mainly repetative
            if compute_rep(seq) > 0.2:
                return None
            ## remove contigs with short mapping
            # if max(df["qlength"]) < 50:
            #     return None
            ## get coordinates to extract section strings
            qse = list(df[["qstart", "qend"]].itertuples())
            seqs = []
            for t in qse:
                seqs.append((seq[t[1]:t[2]], df.loc[t[0]]["mapq"], t[0]))
            nseqs = len(seqs)
            ## for each section check repeatedness of low mapq
            for p, (s, m, i) in enumerate(seqs):
                if m < 10:
                    bp, bpv = base_perc(s)
                    if compute_rep(s) > 0.2 or bpv > 0.9:
                        repeat = list(bp.keys())[0]+list(bp.keys())[1]
                        #print("\n========================\nfound mapq < 10")
                        #print("repeat: ", repeat)
                        #print_details(df=df, seq=seq, gap_seq=s, mapq=mapq, bp=bp, bpv=bpv)
                        #plot_individual((qname, df), cd, show_plot=True)
                        ## if section first 
                        sb, ssb = False, False
                        if p == 0:
                            for x in seqs[1:]:  # check all next seqs
                                if x[1] > 0:   # if mapq decent
                                    cs = x[0]   # compare string
                                    sb = simple_check(repeat, cs)
                                    if sb:
                                        ssb = True
                        ## if section middle
                        elif 0 < p < nseqs:
                            for x in seqs[0:p-1]:
                                if x[1] > 0:
                                    cs = x[0]
                                    sb = simple_check(repeat, cs)
                            for x in seqs[p+1:nseqs]:
                                if x[1] > 0:
                                    cs = x[0]
                                    sb = simple_check(repeat, cs)
                                    if sb:
                                        ssb = True
                        ## if last
                        else:
                            for x in seqs[:-1]:
                                if x[1] > 30:
                                    cs = x[0]
                                    sb = simple_check(repeat, cs)
                                    if sb:
                                        ssb = True
                        if ssb:
                            #print_details(df=df, seq=seq, gap_seq=s, mapq=mapq, bp=bp, bpv=bpv)
                            df = df.drop([i], axis=0)
                            #print_details(df=df, seq=seq, gap_seq=s, mapq=list(df["mapq"]), bp=bp, bpv=bpv)
                            #plot_individual((qname, df), cd, show_plot=True)
    return df




def parse_bed(fn):
    bed = pd.read_csv(fn, sep="\t")
    bed["qlength"] = bed["qend"] - bed["qstart"]
    bed[["sample", "1", "2", "3", "4", "cov"]] = bed['qname'].str.split('_', expand=True)
    bed = bed.drop(["1", "2", "3", "4"], axis=1)
    cd = set_colors()
    #cd["rep"] = "white"
    bed["vcf"] = bed["vcf"].astype(bool) # specify to remove warning
    print(bed)
    flank_segs, flank_contigs, raw_rep, small_map, adj_seg, adj_contigs = 0, 0, 0, 0, 0, 0
    filtered = pd.DataFrame(columns=list(bed.columns))
    filtered["vcf"] = filtered["vcf"].astype(bool) # specify to remove warning
    for df in bed.groupby("qname"):
        filt = tandem_finder(df[1], df[0], cd)
        if isinstance(filt, pd.DataFrame):
            filtered = pd.concat([filtered, filt]).reset_index(drop=True)

    print(filtered)
    filtered.to_csv("out_data/filtered_contigs.bed", sep="\t", index=False)

parse_bed("out_data/contigs.bed")

