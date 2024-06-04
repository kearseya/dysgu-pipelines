import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.collections import BrokenBarHCollection
import matplotlib.cm as cm
import matplotlib.transforms as mtransforms
import matplotlib.ticker as mticker
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import random
import itertools
import glob
import sys
from collections import defaultdict
import networkx as nx
import seaborn as sns
from scipy import stats
import pywt
from statsmodels.robust import mad
from scipy.stats.mstats import winsorize
from scipy.signal import butter, filtfilt, find_peaks
from pathlib import Path

#import pwlf
import os
import numbers
import glob
import pickle
from rich.progress import Progress
from scipy.interpolate import RBFInterpolator
from pypdf import PdfMerger
from itertools import groupby
# import rpy2.robjects as ro
# from rpy2.robjects.packages import importr
# from rpy2.robjects import pandas2ri
import click

from kneed import KneeLocator
from sklearn.datasets import make_blobs
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

@click.group(name='tools', invoke_without_command=True)
@click.option("-i", "--indir", default="raw_cov/hawk", help="input dir with raw coverages")
@click.option("-r", "--ref", default=None)
@click.option("-w", "--winsize", default=10000, help="window size to generate files for")
@click.option("-b", "--bg", nargs=4, default=("sample_pairs.csv", "tumor_db", "tumor_stela", "normal_db"), help="pairs csv with stela [file, tumour col, tumour stela, normal col]")
@click.option("-s", "--subsample", default=250, help="n samples to plot")
@click.option("-t", "--threshold", default=3.81, help="threshold to split groups")
@click.option("--split", default=0.5, help="percentage split of groups in subsample to plot")
@click.option("--subsets", default=1, help="n subsets of size subsample")
@click.option("-c", "--categorical", is_flag=True, default=False, help="Use if split categorical")
@click.option("-t", "--tag", default="dataset", help="tag to identify dataset")
@click.option("-g", "--gamma", default=1000, help="pcf gamma (lower = more relaxed)")
@click.option("-k", "--kmin", default=5, help="pcf kmin (smallest segment size = kmin*winsize)")
@click.option("-p", "--plot-inter", default=False, is_flag=True, help="plot output from intermediary steps")
@click.option("--name-crop", default=3, type=int, help="y-tick name crop int", hidden=True)
@click.option("--prefix", default="None", type=str, help="y-tick prefix remove", hidden=True)
@click.option("--sigfig", default=2, type=int, help="y-tick significant figures", hidden=True)
@click.option("--yaxis-scale", default=False, is_flag=True, help="replace y-tick text with color scale")
@click.option("--ignore-score", default=False, is_flag=True, help="remove score boxplot from perif")
@click.option("--ignore-box", default=False, is_flag=True, help="dont plot score box in perif")
@click.option("--skip-genome", default=False, is_flag=True, help="skip plotting individual genome coverages")
@click.pass_context
def cli_tools(ctx, indir, ref, winsize, bg, subsample, threshold, split, subsets, categorical, tag, gamma, kmin, plot_inter, name_crop, prefix, sigfig, yaxis_scale, ignore_score, ignore_box, skip_genome):
    """Tool related commands"""
    if ctx.invoked_subcommand is None:
        if ref == None:
            click.echo(ctx.get_help())
            ctx.exit()
            pass
        else:
            click.echo("Running pipeline")
            ctx.ensure_object(dict)
            ctx.obj["winsize"] = winsize
            newwin = winsize
            winunit = "bp"
            if winsize//1000 >= 1:
                newwin = winsize//1000
                winunit = "kb"
                if winsize//1000000 >= 1:
                    newwin = winsize//1000000
                    winunit = "mb"
            # analysis
            ctx.obj["newwin"] = newwin
            ctx.obj["winunit"] = winunit
            ctx.obj["raw_cov_dir"] = indir
            ctx.obj["ref"] = ref
            ctx.obj["gamma"] = gamma
            ctx.obj["kmin"] = kmin
            # indput data
            ctx.obj["bgfile"] = bg[0]
            ctx.obj["bgtc"] = bg[1]
            ctx.obj["bgts"] = bg[2]
            ctx.obj["bgnc"] = bg[3]
            # plot params
            ctx.obj["sub"] = subsample
            ctx.obj["thresh"] = threshold
            ctx.obj["split"] = split
            ctx.obj["ssets"] = subsets
            ctx.obj["cat"] = categorical
            ctx.obj["tag"] = tag
            ctx.obj["plot_inter"] = plot_inter
            ctx.obj["namecrop"] = name_crop
            ctx.obj["prefix"] = prefix
            ctx.obj["sigfig"] = sigfig
            ctx.obj["yaxis_scale"] = yaxis_scale
            ctx.obj["ignore_score"] = ignore_score
            ctx.obj["ignore_box"] = ignore_box
            if skip_genome == False:
                ctx.obj["skip_genome"] = None
            else:
                ctx.obj["skip_genome"] = True
            ctx.invoke(normalise)
            ctx.invoke(compare)
            ctx.invoke(plot)
            ctx.invoke(gainloss)
            ctx.invoke(mergeout)
            #compare(ctx, ref)
    else:
        ctx.ensure_object(dict)
        newwin = winsize
        winunit = "bp"
        if winsize//1000 >= 1:
            newwin = winsize//1000
            winunit = "kb"
            if winsize//1000000 >= 1:
                newwin = winsize//1000000
                winunit = "mb"
        # analysis
        ctx.obj["winsize"] = newwin
        ctx.obj["winunit"] = winunit
        ctx.obj["raw_cov_dir"] = indir
        ctx.obj["ref"] = ref
        ctx.obj["gamma"] = gamma
        ctx.obj["kmin"] = kmin
        # input data
        ctx.obj["bgfile"] = bg[0]
        ctx.obj["bgtc"] = bg[1]
        ctx.obj["bgts"] = bg[2]
        ctx.obj["bgnc"] = bg[3]
        # plot params
        ctx.obj["sub"] = subsample
        ctx.obj["thresh"] = threshold
        ctx.obj["split"] = split
        ctx.obj["ssets"] = subsets
        ctx.obj["cat"] = categorical
        ctx.obj["tag"] = tag
        ctx.obj["plot_inter"] = plot_inter
        ctx.obj["namecrop"] = name_crop
        ctx.obj["prefix"] = prefix
        ctx.obj["sigfig"] = sigfig
        ctx.obj["yaxis_scale"] = yaxis_scale
        ctx.obj["ignore_score"] = ignore_score
        ctx.obj["ignore_box"] = ignore_box
        if skip_genome == False:
            ctx.obj["skip_genome"] = None
        else:
            ctx.obj["skip_genome"] = True
    pass


def get_backgrounds(f="sample_pairs.csv", tc="tumor_db", ts="tumor_stela", nc="normal_db"):
    df = pd.read_csv(f)
    df = df.rename({tc: "tumour", nc: "normal", ts: "stela"}, axis=1)
    df["tumour"] = df["tumour"].astype(str)
    df["normal"] = df["normal"].astype(str)
    df = df.dropna()
    tl = df.set_index("tumour")["stela"].to_dict()
    return dict(zip(df["tumour"], df["normal"])), tl



@cli_tools.command(name="preprocess")
@click.option("-r", "--ref", default=None)
@click.option("-w", "--winsize", default=10000, help="window size to generate files for")
@click.option("-c", "--call", default=False, is_flag=True, help="run bash commands")
@click.option("-g", "--gc", default=False, is_flag=True, help="create gc_pct.bed")
@click.option("-m", "--medmap", default=None, help="File to convert to median mappability file")
@click.pass_context
def preprocess(ctx, ref, winsize, call, gc, medmap):
    if ref == None:
        print("please provide reference type, if not present, add lengths below")
        quit()
    elif ref.lower() in {"hg38", "grch38"}:
        chrom_lengths_dict={'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                        'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422, 
                        'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 
                        'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167, 
                        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415} 
    elif ref.lower() in {"hg19", "grch37"}:
        chrom_lengths_dict={'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260, 
                        'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747, 
                        'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 
                        'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520, 
                        'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566}
  
    newwin = winsize
    winunit = "bp"
    if winsize//1000 >= 1:
        newwin = winsize//1000
        winunit = "kb"
        if winsize//1000000 >= 1:
            newwin = winsize//1000000
            winunit = "mb"

    if gc:
        if os.path.exists(f"{ref}_prep.bed") == False:
            with open(f"{ref}_prep.bed", "w") as f:
                for c in chrom_lengths_dict:
                    for i in np.arange(0, chrom_lengths_dict[c], winsize)[:-1]:
                        f.write(f"{c}\t{i}\t{i+winsize}\n")
                    f.write(f"{c}\t{np.arange(0, chrom_lengths_dict[c], winsize)[-1]}\t{chrom_lengths_dict[c]}\n")
            print(f"Created: {ref}.{newwin}{winunit}_windows.gc_pct.bed")

        if call:
            print(f"Running command: bedtools nuc -fi {ref}.fa -bed {ref}_prep.bed | cut -f1,2,3,5 > {ref}.{newwin}{winunit}_windows.gc_pct.bed")
            os.system(f"bedtools nuc -fi {ref}.fa -bed {ref}_prep.bed | cut -f1,2,3,5 > {ref}.{newwin}{winunit}_windows.gc_pct.bed")
        else:
            print(f"run to create gc perc bed:\n\nbedtools nuc -fi {ref}.fa -bed {ref}_prep.bed | cut -f1,2,3,5 > {ref}.{newwin}{winunit}_windows.gc_pct.bed\n\n")

    if medmap != None:
        print("Reading mapping file")
        m = pd.read_csv(medmap, sep="\t")
        chroms = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
        total_regions = 0
        chrom_regions = []
        for chrom in chroms:
            num_regions = len(np.arange(0, chrom_lengths_dict[chrom], winsize))
            total_regions += num_regions
            chrom_regions.append(num_regions)
        chrom_cumsum = np.cumsum(chrom_regions)
        bin_tag = {'chr1': 0}
        for idx, c in enumerate(chroms[1:]):
            bin_tag[c] = chrom_cumsum[idx]
        print(bin_tag)
        avg_df = pd.read_csv(f"{ref}_prep.bed", sep="\t", header=None)
        avg_df = avg_df.rename({0: "chrom", 1: "start", 2: "end"}, axis=1)
        
        m["bin"] = m["chrom"].map(bin_tag)
        avg_df["bin"] = avg_df["chrom"].map(bin_tag)
        m["bin"] = m["bin"] + m["start"]//winsize
        avg_df["bin"] = avg_df["bin"] + avg_df["start"]//winsize

        medians = m.groupby("bin")["map"].median()
        medians = medians.to_frame()
        medians = medians.reset_index()
        medians = medians.rename({"map": "median"}, axis=1)
        print("medians", medians)
        means = m.groupby("bin")["map"].mean()
        means = means.to_frame()
        means = means.reset_index()
        means = means.rename({"map": "mean"}, axis=1)
        print("means", means)
        counts = m.groupby("bin")["map"].size()
        counts = counts.to_frame()
        counts = counts.reset_index()
        counts = counts.rename({"map": "n"}, axis=1)
        print("counts", counts)

        avg_df = avg_df.merge(medians, on="bin", how="left")
        avg_df = avg_df.merge(means, on="bin", how="left")
        avg_df = avg_df.merge(counts, on="bin", how="left")

        avg_df[["median", "mean"]] = avg_df[["median", "mean"]].fillna(-1)
        avg_df["n"] = avg_df["n"].fillna(0)
        avg_df = avg_df.astype({"n": int})
        print(avg_df)

        # with Progress() as progress:
        #     l = []
        #     task_med = progress.add_task("Calculating average mappability:", total=total_regions)
        #     for chrom in [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]:
        #         cm = m[m["chrom"] == chrom]
        #         for i in np.arange(0, chrom_lengths_dict[chrom], winsize)[:-1]:
        #             tmp = cm[(cm["start"].between(i, i+winsize) & cm["end"].between(i, i+winsize))]
        #             if len(tmp) == 0:
        #                 l.append([chrom, i, i+winsize, -1, -1, 0])
        #             else:
        #                 l.append([chrom, i, i+winsize, np.median(tmp["map"]), np.mean(tmp["map"]), len(tmp)])
        #             progress.update(task_med, advance=1)
        #         last_window = np.arange(0, chrom_lengths_dict[chrom], winsize)[-1]
        #         tmp = cm[cm["start"].between(last_window, chrom_lengths_dict[chrom])]
        #         if len(tmp) == 0:
        #             l.append([chrom, last_window, chrom_lengths_dict[chrom], -1, -1, 0])
        #         else:
        #             l.append([chrom, last_window, chrom_lengths_dict[chrom], np.median(tmp["map"]), np.mean(tmp["map"]), len(tmp)])
        #         progress.update(task_med, advance=1)
        # avg_df = pd.DataFrame(l, columns=["chrom", "start", "end", "map", "mean", "n"])

        if ref.lower() in {"hg19", "grch37"}:
            avg_df.to_csv(f"hg19.{newwin}{winunit}_windows.median_mappability.bed", sep="\t")
        elif ref.lower() in {"hg38", "grch38"}:
            avg_df.to_csv(f"hg38.{newwin}{winunit}_windows.median_mappability.bed", sep="\t")



@cli_tools.command(name="normalise")
@click.pass_context
def normalise(ctx):
    """
    Generates the x.cov.bed files (1st)

    Next run the compare_depth_gc_map_normalized.py
    """
    sns.set_style("whitegrid")
    ref = ctx.obj["ref"] 
    normal_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                   'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    if ref.lower() in {"hg38", "grch38"}:
        df = pd.read_csv(f"hg38.{ctx.obj['newwin']}{ctx.obj['winunit']}_windows.median_mappability.bed", sep="\t")#, header=None)
        gc_df = pd.read_csv(f"hg38.{ctx.obj['newwin']}{ctx.obj['winunit']}_windows.gc_pct.bed", sep="\t")
    elif ref.lower() in {"hg19", "grch37"}:
        df = pd.read_csv(f"hg19.{ctx.obj['newwin']}{ctx.obj['winunit']}_windows.median_mappability.bed", sep="\t")#, header=None)
        gc_df = pd.read_csv(f"hg19.{ctx.obj['newwin']}{ctx.obj['winunit']}_windows.gc_pct.bed", sep="\t")
    #df.columns = ["chrom", "start", "end", "map"]
    gc_df = gc_df[gc_df["#chrom"].isin(normal_list)]
    df["GC"] = gc_df["pct_gc"]
    print(gc_df)
    print(df)
    assert np.array_equal(df["start"], gc_df["start"]), "GC and mapping have differnt lengths"

    # Round the map value
    df["map"] = (np.round(df["map"], 3) * 100).astype(int)
    df["GC"] = (df["GC"] * 100).astype(int)

    print(list(df["map"][:20]))
    print(list(df["GC"][:20]))
    print(df.head())

    map_gc_dict = {k: v for k, v in zip(zip(df["chrom"], df["start"], df["end"]), zip(df["map"], df["GC"]))}
    print(list(map_gc_dict.items())[:10])

    normal_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                   'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    # normal_list = [i.replace("chr", "") for i in normal_list]
    if ref.lower() in {"hg38", "grch38"}:
        chrom_lengths = (248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 
                         138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 
                         83257441, 80373285, 58617616, 64444167, 46709983, 50818468, 156040895, 57227415)
    elif ref.lower() in {"hg19", "grch37"}:
        chrom_lengths = (249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022,
                         141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753,
                         81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566)
    chrom_sizes = dict(zip(normal_list, chrom_lengths))


    # Load in read counts for samples
    sample_reads = glob.glob(f"{ctx.obj['raw_cov_dir']}/*_cov.bed")
    print(len(sample_reads), sample_reads)

    def flex_gap(df, chrom, window_size):
        pd.options.mode.chained_assignment = None
        df["denom"] = df["start"] // window_size
        df["gap"] = df.start - df.end.shift()
        df["gap"] = df["gap"].astype(bool).astype(int)
        df["gap"] = df.gap.cumsum()
        df["group"] = df["denom"] + df["gap"]
        new = df.groupby(df["group"]).agg({"start": "min", "end": "max", "coverage": "mean"})
        new["size"] = new["end"] - new["start"]
        new = new[new["size"] == window_size]
        new["chromosome"] = chrom
        new = new[["chromosome", "start", "end", "coverage"]]
        return new

    def load(s, map_gc_dict):
        data = pd.read_csv(s, sep="\t", low_memory=False)
        data.columns = ["chromosome", "start", "end", "coverage"]
        data["chromosome"] = data["chromosome"].apply(lambda x: "chr"+x if "chr" not in x else x)
        # if len(data[data["chromosome"].str.contains("chr")==False]) > 0:
        #     data["chromosome"] = "chr" + data["chromosome"].astype(str)

        new = pd.DataFrame()
        max_df = pd.DataFrame(map_gc_dict.keys())
        max_df = max_df.rename({0: "chromosome", 1: "start", 2: "end"}, axis=1)
        max_df = max_df[max_df["chromosome"].isin(normal_list)]
        max_df = max_df.loc[max_df.groupby("chromosome")["end"].idxmax()]
        for chrom in normal_list:
            length = max_df[max_df["chromosome"] == chrom]["end"].values[0]
            tmp = flex_gap(data[data["chromosome"] == chrom], chrom, ctx.obj["winsize"])
            try:
                for _ in range(10):
                    if tmp.iloc[-1]["end"] > length:
                        tmp = tmp.iloc[:-1 , :]
                    else:
                        break
            except:
                print("trimming failed")
                continue
            new = pd.concat([new, tmp])
        data = new

        mapb, gcpct = zip(*[map_gc_dict[k] for k in zip(data["chromosome"], data["start"], data["end"])])

        data["GC"] = gcpct
        data["map"] = mapb
        return data

    def norm(s, df, sample):
        samp = sample.split("/")[-1].split(".")[0].replace("_cov", "")
        print(samp)
        norm_array = defaultdict(list)
        for c, g, m in zip(list(s["coverage"]), list(df["GC"]), list(df["map"])):

            if c > 500:
                continue

            if 30 < g < 60:
                if 60 < m < 101:
                    norm_array[(g, m)].append(c)
            # if 0 < g < 100:
            #     if 0 < m < 101:
            #         norm_array[(g, m)].append(c)

        s = s[(s["GC"].between(30, 60)) & (s["map"].between(70, 101))]
        #s = s[(s["GC"].between(0, 101)) & (s["map"].between(0, 101))]
        # Genome median
        median = np.median(s["coverage"])
        print("Median", median)
        # The median value at this GC and Mappabilty value
        norm_array = {k: np.median(np.array(v)) for k, v in norm_array.items()}

        x, y, z = [], [], []
        arr = np.zeros((101, 101))
        for (i, j), v in norm_array.items():
            if v < 1000:
                sub = v - median
                arr[i, j] = sub
                x.append(j), y.append(i), z.append(sub)

        ## legacy so might get removed
        # rbf = Rbf(x, y, z)
        # z = rbf(x, y)
        ## calculate epsilon as previous for same kernel
        xi = np.stack([x, y])
        ximax = np.amax(xi, axis=1)
        ximin = np.amin(xi, axis=1)
        edges = ximax - ximin
        edges = edges[np.nonzero(edges)]
        epsilon = np.power(np.prod(edges)/len(x), 1.0/len(edges))
        
        rbf = RBFInterpolator(np.column_stack([x, y]), np.array(z), smoothing=5, kernel="multiquadric", epsilon=epsilon)
        z = rbf(np.column_stack([x, y]))

        if ctx.obj["plot_inter"] == True:
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            ax.set_zlabel("Difference from median cov")
            p3d = ax.scatter(x, y, z, s=30, c=z, cmap=cm.coolwarm, vmin=-1, vmax=1)
            plt.ylabel("GC %")
            plt.xlabel("Mappability %")
            plt.title(f"{samp} Mappability vs GC")
            ax.set_zlim(-4, 4)
            plt.savefig(os.path.join("output", ctx.obj["tag"], "map_gc_median_diff", f"{samp}_map_gc_median_diff.thresholded.pdf"))
            plt.close()
            #plt.grid(False)
            #quit()
        # Do the normalization
        norm_value = {(int(xx), int(yy)): zz for xx, yy, zz in zip(x, y, z)}

        normed_counts = []
        bad_indexes = []
        count = 0
        for g, m, v in zip(list(s["GC"]), list(s["map"]), list(s["coverage"])):
            if (m, g) in norm_value:
                normed_counts.append(v - norm_value[(m, g)])
            else:
                normed_counts.append(0)
                bad_indexes.append(count)
            count += 1

        # Now do a wavelet transform to get final step
        d = np.array(normed_counts)

        noisy_coefs = pywt.wavedec(d, 'haar', level=2, mode='per')
        sigma = mad(noisy_coefs[-1])
        uthresh = sigma * np.sqrt(2 * np.log(len(d)))

        denoised = noisy_coefs[:]
        denoised[1:] = (pywt.threshold(i, value=uthresh) for i in denoised[1:])
        sig = pywt.waverec(denoised, 'haar', mode='per')

        if len(sig) > len(s["GC"]):
            sig = sig[:-1]

        bad_indexes = set(bad_indexes)  # Remove bad indexes from final result
        s["normed"] = [sig[i] if i not in bad_indexes else np.NAN for i in range(len(sig))]  #normed_counts

        print("Result on stdev after normalizing before/after:", s["coverage"].std(), s["normed"].std())

        if ctx.obj["plot_inter"] == True:
            all_chroms = {c: df for c, df in s.groupby("chromosome")}
            for chrom in ["chr{}".format(i) for i in range(1, 23)] + ["chrX"]:

                d = all_chroms[chrom]

                # Plot the result of normalizing
                fig, ax = plt.subplots(figsize=(7,2))
                plt.subplots_adjust(bottom=0.25)
                plt.title("{}: {}".format(samp, chrom))

                plt.scatter(d["start"], d["normed"], s=2, color="black", rasterized=True)
                plt.ylim(0, 60)
                plt.xlabel("Chromosome position (Mb)")
                plt.ylabel("Coverage")

                plt.xlim(0, d.end.max())
                #ax.yaxis.grid(True, lw=1)
                ax.set_yticks(np.arange(0, 60, 10))
                #ax.xaxis.grid(False)

                plt.xticks(plt.xticks()[0], [i / 1e6 for i in plt.xticks()[0]])
                sampledir = os.path.join("output", ctx.obj["tag"], "norm_graphs_no_bg_subtract", samp)
                outdir = os.path.join("output", ctx.obj["tag"], "norm_graphs_no_bg_subtract", samp)
                if not os.path.exists(sampledir):
                    os.mkdir(sampledir)
                if not os.path.exists(outdir):
                    os.mkdir(outdir)
                plt.savefig(os.path.join(outdir, f"{samp}_{chrom}.after_norm.pdf"), dpi=200)
                plt.close()

            # Plot a whole genome profile
            fig, ax = plt.subplots(figsize=(9, 2))
            plt.subplots_adjust(bottom=0.25)
            plt.title("{}".format(samp[2:]))
            plt.ylim(0, 60)

            plt.ylabel("Coverage")
            plt.xlabel("Chromosome")
            plt.xticks(plt.xticks()[0], "")
            plt.grid(False)
            coor = itertools.cycle(["black", "chocolate"])
            max_pos = 0
            for idx, chrom in enumerate(["chr{}".format(i) for i in range(1, 23)] + ["chrX"]):
                d = all_chroms[chrom]

                if idx > 0:
                    index = sum(chrom_lengths[:idx])
                else:
                    index = 0
                color = next(coor)
                plt.plot([index, index], [0, 100], "--", color="black", lw=0.5, alpha=0.3)
                max_pos = (d.start + index).max()
                plt.scatter(d["start"] + index, d["normed"], s=0.5, color=color, rasterized=True)

                if color == "black":
                    plt.text(index + 0.5*chrom_lengths[idx], 50, chrom.strip("chr"), horizontalalignment='center', fontsize=10)
                else:
                    plt.text(index + 0.5 * chrom_lengths[idx], 4, chrom.strip("chr"), horizontalalignment='center',
                             fontsize=10)
            plt.xlim(0, max_pos)
            plt.savefig(os.path.join(outdir, f"{samp}.All.after_norm.pdf"), dpi=200)
            plt.close()

        # sys.exit()
        # Save normalized coverage values
        s = s[["chromosome", "start", "end", "normed"]]    
        s.to_csv(os.path.join("cnpinter", ctx.obj["tag"], "normalized_gc_map_wav", f"{samp}.cov.bed"), sep="\t", index=False)

    if not os.path.exists(os.path.join("cnpinter")):
        os.mkdir(os.path.join("cnpinter"))
    if not os.path.exists(os.path.join("cnpinter", ctx.obj["tag"])):
        os.mkdir(os.path.join("cnpinter", ctx.obj["tag"]))
    if not os.path.exists(os.path.join("output", ctx.obj["tag"])):
        os.mkdir(os.path.join("output", ctx.obj["tag"]))
    if not os.path.exists(os.path.join("cnpinter", ctx.obj["tag"], "normalized_gc_map_wav")):
        os.mkdir(os.path.join("cnpinter", ctx.obj["tag"], "normalized_gc_map_wav"))

    norm_no_bg_sub_path = os.path.join("output", ctx.obj["tag"], "norm_graphs_no_bg_subtract")
    map_gc_med_path = os.path.join("output", ctx.obj["tag"], "map_gc_median_diff")
    if ctx.obj["plot_inter"] == True:
        if not os.path.exists(norm_no_bg_sub_path):
            os.mkdir(norm_no_bg_sub_path)
        if not os.path.exists(map_gc_med_path):
            os.mkdir(map_gc_med_path)
    
    # Now normalize read counts for each sample
    for s in sample_reads:
        samp = Path(s).stem.replace("_cov", "")
        if os.path.isfile(os.path.join("cnpinter", ctx.obj["tag"], "normalized_gc_map_wav", f"{samp}.cov.bed")) == False:
            data = load(s, map_gc_dict)
            norm(data, df, s)
        else:
            print(os.path.join("cnpinter", ctx.obj["tag"], "normalized_gc_map_wav", f"{samp}.cov.bed"), " exists")



@cli_tools.command(name="compare")
@click.pass_context
def compare(ctx):
    """
    Plots normalised read counts (2nd)

    Plots the normalized read counts. Counts have been normalized to GC and mappability using the normalize script.
    Makes the input for the copyNumber package
    See ./output/graphs/{}/{}_normalized.pdf for graphs

    Runs Rscript run_copyNumber.R
    """
    if not os.path.exists(os.path.join("cnpinter", ctx.obj["tag"], "copyNumber_io")):
        os.mkdir(os.path.join("cnpinter", ctx.obj["tag"], "copyNumber_io"))

    plt.ioff()

    plots = True

    normal_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
    'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    ref = ctx.obj["ref"]
    if ref.lower() in {"hg38", "grch38"}:
        chrom_lengths_dict={'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                        'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422,
                        'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189,
                        'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167,
                        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415}
    elif ref.lower() in {"hg19", "grch37"}:
        chrom_lengths_dict={'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260,
                        'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747,
                        'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392,
                        'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520,
                        'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566}


    chrom_lengths = [chrom_lengths_dict[i] for i in normal_list]

    colors = ['#FB141B', '#FB461D', '#FC8524', '#FEC72E', '#EFEF17', '#B8FD34', '#7DFC31', '#4EFC30', '#41FC36', '#41FD64', '#40FD9F',
              '#3EFEDF', '#32DFFE', '#1A9EFB', '#0063FC', '#0034FC', '#232DFC', '#6A2DFC', '#AD30FC', '#EF32FD', '#FA29C6', '#FB1E86',
              '#FB1748',  '#FB141B']

    colors = ['#CA1F1F', '#051F8A'] * 12

    # Background is tumor: normal
    backgrounds, tumour_tl = get_backgrounds(ctx.obj["bgfile"], ctx.obj["bgtc"], ctx.obj["bgts"], ctx.obj["bgnc"])
    print(backgrounds)
    print(tumour_tl)

    binsize=ctx.obj["winsize"]
    read_length = 125.0

    f = glob.glob(os.path.join("cnpinter", ctx.obj["tag"], "normalized_gc_map_wav")+"/*.cov.bed")
    print(len(f), f)

    def get_coverage(d):
        df = pd.read_csv(d, sep="\t")
        median_cov = np.nanmedian(np.array(df.normed))
        df["norm2median"] = df["normed"] / median_cov * 2.
        return df, median_cov  # Skip any NaNs

    def g_medians(l):
        chrm = [item for sublist in l for item in sublist]
        return np.median(chrm)

    def mung_data():
        dfs = {}
        dds = []  # Hold x and y data for plotting

        signals = {}
        read_counts = {}
        x_all = {}
        y_all = {}

        for d in f:
            samp = d.split("/")[-1].split(".")[0].replace("_cov", "")
            #samp = samp.replace("_cov", "")

            df, median_cov = get_coverage(d)
            print("'{}': {},".format(samp, median_cov))
            grouped = {chrm: v for chrm, v in df.groupby("chromosome")}
            if all([i in grouped for i in normal_list]) == False:
                print(f"chromosome missing in {samp}")
                continue

            sigs = []
            r_c = []
            xc = []
            yc = []
            for c in normal_list:
                m = grouped[c]

                y = m['norm2median']
                x = m['start']

                sigs.append(y)
                r_c.append([x, y])
                yc.append(y)

            signals[samp] = sigs
            dfs[samp] = df
            read_counts[samp] = r_c

            x_all[samp] = xc
            y_all[samp] = yc

        database = [dfs, dds, signals, read_counts, x_all, y_all]
        #pickle.dump(database, open("database.pkl", "wb"), protocol=2)  # Protocol 2 is much faster for loading
        return database

    dfs, dds, signals, read_counts, x_all, y_all = mung_data()
    #print("dfs", dfs.keys())
    ##print("dds", dds[0])
    #print("signals", signals.keys())
    #print("read_counts", read_counts.keys())
    #print("x_all", x_all.keys())
    #print("y_all", y_all.keys())
    #print("Loading DB")
    #dfs, dds, signals, read_counts, x_all, y_all = pickle.load(open("database.pkl", "rb"))
    ##print("DFS:\n", dfs)#, "\nDDS:\n", dds, "\nSIGNALS:\n", signals, "\nREAD COUNTS\n:", read_counts)
    #print("DB loaded")
    #print(tumour_tl)
    #print(signals)

    def make_copynumber_input():
        df_cn = pd.DataFrame()
        # Save a file for the copynumber R package
        # Order by telomere length, shortest to longest
        #print(backgrounds)
        samples_order = sorted(tumour_tl.items(), key=lambda x: x[1])
        #print(samples_order)

        count = 0
        for samp, tl in samples_order:
            # catch missing from backgrounds
            try:
                parental = backgrounds[samp]
                if parental not in signals or samp not in x_all:
                    print("WARNING not in x_all")
                    continue
            except:
                print("Parental missing")
                continue

            p_df = dfs[parental]
            #print(p_df)
            c_df = dfs[samp]
            #print(c_df)

            # Add generic bed info
            if count == 0:
                df_cn["chromosome"] = c_df["chromosome"]
                df_cn["pos"] = c_df["start"] + 5000  # The middle point
                count += 1

            # Sample subtract the parental background:
            df_cn[samp] = c_df["norm2median"] - p_df["norm2median"]  # np.log(c_df["norm2median"] / p_df["norm2median"])

        final_order = samples_order

        # Save a copy of the order in this plot
        with open(os.path.join("cnpinter", ctx.obj["tag"], "copyNumber_io", "sample_order.txt"), "w") as fo:
            fo.write("\t".join([str(v[0]) for v in final_order]))

        # Drop NaN values
        print(df_cn)
        df_cn = df_cn.dropna()
        print(df_cn)

        df_cn["chromosome"] = [i.strip("chr") for i in list(df_cn["chromosome"])]
        df_cn = df_cn[df_cn["chromosome"] != "Y"]  # Skip Y chromosome

        # Sort chromosomes
        chroms = {k: k.strip("chr") for k in set(df_cn["chromosome"])}
        chroms["X"] = '23'
        chroms = {k: int(v) for k, v in chroms.items()}

        df_cn["sortkey"] = [chroms[i] for i in list(df_cn["chromosome"])]
        df_cn = df_cn.sort_values(by=["sortkey", "pos"])
        del df_cn["sortkey"]

        # df_cn = df_cn[["chromosome", "pos", "DB37", "DB41"]]  # Keep a small selection for now for testing
        ## Remove DB from column names
        ## df_cn.columns = [i.strip("DB") for i in df_cn.columns]
        df_cn.to_csv(os.path.join("cnpinter", ctx.obj["tag"], "copyNumber_io", "copynumbers.all.csv"), index=False)

    def subtract_bg():
        for samp, parental in backgrounds.items():
            #print(samp, signals, samp, x_all)
            print(samp)
            if parental not in signals or samp not in x_all:
                print("WARNING parental or sample missing")
                continue
            p_sig = signals[parental]
            #p_df = dfs[parental]
            p_y = y_all[parental]

            #x = x_all[samp]
            y = y_all[samp]
            #sig = signals[samp]

            #y = [a - b for a, b in zip(y, p_y)]

            colors1 = ['#CA1F1F', '#051F8A'] * 12

            l = -2
            h = 3

            fig = plt.figure(figsize=(12,3))
            plt.title("{}, {}".format(samp, tumour_tl[samp]), fontsize=14)
            plt.xlim(0, sum(chrom_lengths))

            plt.ylabel("Relative copy number change", fontsize=12)
            plt.xlabel("Chromosome number", fontsize=12)

            plt.tick_params(top=False, bottom=False, left=True, right=False,
                            labelleft=True, labelbottom=True)

            ax = fig.add_subplot(111)

            # ax.spines['left'].set_color('black'), ax.spines['right'].set_color('#ECECEC')
            # ax.spines['top'].set_color('#ECECEC'), ax.spines['bottom'].set_color('black')
            # ax.axes.set_xticks([])
            # for tic in ax.yaxis.get_major_ticks():  # Remove ticks but keep labels
            #     tic.tick1On = tic.tick2On = False
            # ax.yaxis.grid(linestyle='-', linewidth=0.5, alpha=0.2)

            for i, key in enumerate(normal_list):
                x = read_counts[samp][i][0]
                y = read_counts[samp][i][1]
                ysig = signals[samp][i]

                mm = p_sig[i] #mm_chroms[i]

                if len(mm) != len(ysig):  # Sometimes has an extra value
                    mm = mm[:-1]

                p = sum([chrom_lengths[k] for k in range(i)])

                endx = [chrom_lengths[i]+p]*2
                endy = [l, h]

                ysig = ysig - mm
                if x.shape != ysig.shape:
                    ysig = ysig[:-1]
                ax.text(endx[0] - chrom_lengths[i]/2, -1.5, normal_list[i][3:], fontsize=8, horizontalalignment='center')

                plt.scatter(x + p, ysig, s=1, color=colors1[i], rasterized=True) # facecolor=colors[i], edgecolor="black")
                plt.plot(endx, endy, '--', dashes=(1,2), color="black", linewidth=1)

            if ctx.obj["plot_inter"] == True:
                directory = os.path.join("output", ctx.obj["tag"], "graphs", samp)
                if not os.path.exists(directory):
                    os.makedirs(directory)

                plt.savefig(os.path.join(directory, f"{samp}_subtract_bg.pdf"), dpi=400)
            plt.close()

    if not os.path.exists(os.path.join("output", ctx.obj["tag"], "graphs")):
        os.makedirs(os.path.join("output", ctx.obj["tag"], "graphs"))

    # Subtract background
    make_copynumber_input()
    subtract_bg()


def pcf(wd, assembly, samples, tag, rplot, gamma, kmin):
    if rplot == True:
        if not os.path.exists(os.path.join(wd, "output", tag, "R")):
            os.makedirs(os.path.join(wd, "output", tag, "R"))
        if not os.path.exists(os.path.join(wd, "output", tag, "graphs")):
            os.makedirs(os.path.join(wd, "output", tag, "graphs"))
        for s in samples:
            sampledir = os.path.join(wd, "output", tag, "R", "graphs", s)
            if not os.path.exists(sampledir):
                os.makedirs(sampledir)
    os.system(f"Rscript run_copynumber.R {wd} {assembly} {tag} {rplot}")


def get_start_end_mean(df, c):
    df = df[df["chrom"] == c]
    return list(df[["mean", "start.pos", "end.pos"]].itertuples(index=False, name=None))


def lowpass_score(df, score_thresh=0.5):
    T = 5.0         # sample period (time)
    fs = 30.0       # sample rate (Hz)
    cutoff = 0.5    # desired cutoff frequency of the filter
    nyq = 0.5 * fs  # Nyquist Frequency
    order = 2       # sin wave can be approx represented as quadratic
    n = int(T * fs) # total number of samples

    def butter_lowpass_filter(data, cutoff, fs, order):
        normal_cutoff = cutoff / nyq
        # Get the filter coefficients
        b, a = butter(order, normal_cutoff, btype='low', analog=False)
        y = filtfilt(b, a, data, axis=0)
        return y

    samples = list(df.columns)[2:]
    scores = {}
    for ypos, sample in enumerate(samples):
        scores[sample] = {}
        data = df[["chromosome", sample]]
        data = data.rename({"chromosome": "chrom"}, axis=1)
        chroms = [str(i) for i in range(1, 23)] + ["X"]
        for idx, chrom in enumerate(chroms):
            filt_data = np.array(data[data["chrom"] == chrom][sample].tolist())
            ## low pass
            y = butter_lowpass_filter(filt_data, cutoff, fs, order)
            p, properties = find_peaks(y)
            t, tproperties = find_peaks(-y)

            c = []
            d = []
            if len(p) > len(t):
                c = [i for x in zip(p, t) for i in x] + [p[-1]]
            else:
                c = [i for x in zip(t, p) for i in x] + [t[-1]]
          
            c = list(sorted(c))
            for didx, i in enumerate(c[:-1]):
                d.append(abs(y[c[didx]] - y[c[didx+1]]))
            
            scores[sample][chrom] = int(round(sum([i for i in d if i > score_thresh]), 0))
    
    return scores


def jitterer(df):
    width = 0.8
    for x in df.groupby("sample"):
        data = x[1]["scores"].tolist()
        ypos = x[1]["ypos"].values[0]
        counts, edges = np.histogram(data, bins=int(max(data)-min(data)))
        centres = (edges[:-1] + edges[1:]) / 2.
        yvals = centres.repeat(counts)
        max_offset = width / counts.max()
        offsets = np.hstack([(np.arange(cc) - 0.1 * (cc - 1)) for cc in counts])
        x[1]["ypos_j"] = ypos + (offsets * max_offset)
        df = df.merge(x[1], how="left", on="scores")
    return df


def simple_beeswarm2(y, nbins=None, width=0.1):
    """
    Returns x coordinates for the points in ``y``, so that plotting ``x`` and
    ``y`` results in a bee swarm plot.
    # from https://stackoverflow.com/questions/36153410/how-to-create-a-swarm-plot-with-matplotlib
    """
    y = np.asarray(y)
    if nbins is None:
        # nbins = len(y) // 6
        nbins = np.ceil(len(y) / 6).astype(int)
    # Get upper bounds of bins
    x = np.zeros(len(y))
    nn, ybins = np.histogram(y, bins=nbins)
    nmax = nn.max()
    #Divide indices into bins
    ibs = []#np.nonzero((y>=ybins[0])*(y<=ybins[1]))[0]]
    for ymin, ymax in zip(ybins[:-1], ybins[1:]):
        i = np.nonzero((y>ymin)*(y<=ymax))[0]
        ibs.append(i)
    # Assign x indices
    dx = width / (nmax // 2)
    for i in ibs:
        yy = y[i]
        if len(i) > 1:
            j = len(i) % 2
            i = i[np.argsort(yy)]
            a = i[j::2]
            b = i[j+1::2]
            x[a] = (0.5 + j / 3 + np.arange(len(b))) * dx
            x[b] = (0.5 + j / 3 + np.arange(len(b))) * -dx
    return x

# def rand_jitter(arr):
#     stdev = .01 * (max(arr) - min(arr))
#     return arr + np.random.randn(len(arr)) * stdev


def conv_mathtext(value, bold=True):
    form = mticker.ScalarFormatter(useOffset=False, useMathText=True)
    print(value)
    if value != 0:
        s = form.format_data(float("{:.1E}".format(float(value)))).replace(" ", "")
    else:
        s = "0.0"
    if bold == True:
        return f"$\mathbf{{{s}}}$"
    else:
        return f"${s}$"

def to_sigfig(v, n):
    return float('{:g}'.format(float('{:.{p}g}'.format(v, p=n))))

def bool_to_label(l, tlab="S", flab="L"):
    if l:
        return str(tlab)
    else:
        return str(flab)


def add_periferal_stats(ax, scores, samples, tumour_tl, cat, threshold, ss=1, boty=-2, offy=-1, bottom=1, ignore_score=False, ignore_box=False):
    scores = pd.DataFrame.from_records(scores)
    scores = scores[samples]
    swarm_df = pd.DataFrame({"scores": [], "sample": [], "ypos": []})
    total_scores = []
    for ypos, sample in enumerate(samples):
        swarm_df = pd.concat([swarm_df, pd.DataFrame({"scores": scores[sample].tolist(), "sample": [sample]*23, "ypos": [(ypos*ss)+1+offy]*23})], ignore_index=True)
        total_scores.append(np.sum(scores[sample]))
    print("total scores")
    print(total_scores)
    #swarm_df = pd.concat([swarm_df, pd.DataFrame({"scores": [-10], "sample": ["spacer"]})], ignore_index=True)
    print("scores df")
    print(scores)
    print("swarm df")
    print(swarm_df)
    telomere_lengths = [tumour_tl[i] for i in samples]
    
    # sns.boxplot(ax=ax["score"], data=swarm_df, x="scores", y="sample",
    #         meanprops={'marker' : 'D', 'markeredgecolor' : 'black', 'markersize' : 6, 'label':'mean'},
    #         medianprops={'visible': False}, whiskerprops={'visible': False},
    #         showmeans=True, showfliers=False, showbox=True, showcaps=False, boxprops=dict(alpha=.5), zorder=20, hue="sample", palette="pastel", legend=False) #, transform=mtransforms.Affine2D().translate(0, 1+offy))
    # sns.stripplot(ax=ax["score"], data=swarm_df, x="scores", y="sample", zorder=100, hue="sample", palette="pastel", legend=False, transform=mtransforms.Affine2D().translate(0, -(1+offy)))
   
    score_colors = (1. - 0.0) * plt.get_cmap("tab10")(np.linspace(0., 1., len(samples))) + 0.0 * np.ones((len(samples), 4)) # 0.0 change to add whiteness
    if ignore_score == False:
        for ci, x in enumerate(swarm_df.groupby("sample")):
            # swarm plot
            y = simple_beeswarm2(x[1]["scores"].tolist())
            offset = x[1]["ypos"].values[0]
            ypos = [i+offset for i in y]
            x[1]["ypos_j"] = [i+offset for i in y]
            ax["score"].scatter(x=x[1]["scores"], y=ypos, color=score_colors[ci], s=2)
            # boxplot
            ax["score"].boxplot(x[1]["scores"], positions=[offset], vert=False, patch_artist=True, showmeans=True, widths=(ss-(ss/8)), showfliers=False, boxprops=dict(facecolor=score_colors[ci], alpha=0.5, linewidth=0.5), whiskerprops=dict(linewidth=0.5, linestyle="--"), capprops=dict(linewidth=0.5), medianprops=dict(color="red", linewidth=0.5), meanprops=dict(marker='D', markeredgecolor="black", markerfacecolor="green", markersize=5))

        # ax["score"].scatter(x=swarm_df["scores"], y=swarm_df["ypos_j"])
        # ax["score"].boxplot(swarm_df["scores"], labels=swarm_df["sample"], vert=False, positions=[i+offy for i in range(1, len(samples)+1)])
        
        ax["score"].set_xlim(-0.5, max(swarm_df["scores"])+5)
        #ax["score"].set_ylim(((boty+offy)/ss), len(samples)+offy)
        ax["score"].set_ylim(boty, len(samples)*ss)
        ax["score"].tick_params(left=False, labelleft=False)
        ax["score"].set_xscale("symlog", base=10)
        ax["score"].set(ylabel=None)
        #ax["score"].legend().remove()
        ax["score"].grid(which="major", axis="x")
        ax["score"].spines[["left", "bottom", "top", "right"]].set_edgecolor("black")
        ax["score"].tick_params(axis="x", which="major", bottom=True, labelsize=5, length=3, direction="out")
        ax["score"].tick_params(axis="x", which="minor", bottom=True, length=1.5, direction="out")
        #ax["score"].minorticks_off()

    recursive_partition = {"total": [], "chrom": []}
    for p in range(1, len(total_scores)):
        recursive_partition["total"].append(stats.mannwhitneyu(total_scores[:p], total_scores[p:], alternative="greater", axis=0).pvalue)
        recursive_partition["chrom"].append(stats.mannwhitneyu(np.array(list(itertools.chain(*scores[scores.columns[:p]].values.tolist()))), np.array(list(itertools.chain(*scores[scores.columns[p:]].values.tolist()))), alternative="greater", axis=0).pvalue)
    print(recursive_partition)

    ax2 = ax["rp"].twiny()
    ax2.plot(recursive_partition["chrom"], [(i*ss)+offy+1 for i in range(0, len(samples)-1)], color="C1")
    ax["rp"].plot(recursive_partition["total"], [(i*ss)+offy+1 for i in range(0, len(samples)-1)], color="C0")
    ax["rp"].tick_params(left=False, labelleft=False)
    ax2.tick_params(left=False, labelleft=False)
    #ax["rp"].plot(recursive_partition["chrom"], [(i*ss)-offy for i in range(len(samples)-1)], color="C1")
    ax["rp"].set_ylim(boty, len(samples)*ss) #ax["rp"].set_ylim(boty, (len(samples)*ss))
    ax["rp"].set_xlim(min([float(min(recursive_partition["total"])-(min(recursive_partition["total"])/8)), 0.001]), 1)
    ax2.set_xlim(min([float(min(recursive_partition["chrom"])-(min(recursive_partition["chrom"])/8)), 0.001]), 1)
    #ax["rp"].set_xlim(min(recursive_partition["total"]+recursive_partition["chrom"]), 1)
    ax["rp"].set_xticks([0.005,0.010,0.050,0.100,0.500,1],[0.005,0.010,0.050,0.100,0.500,1])
    ax["rp"].set_xscale("log", base=100)
    ax2.set_xscale("log", base=100)
    ax["rp"].set_xticks(ticks=[1, 0.05, 0.01, 0.001], labels=[1, 0.05, 0.01, 0.001], minor=False)
    ax2.set_xticks([])
    ax["rp"].set_xlabel("pvalue")
    ax["rp"].tick_params(axis="x", labelsize=5, labelrotation=45, bottom=True, top=False, left=False, right=False)
    ax["rp"].minorticks_off()
    ax2.tick_params(axis="x", labelsize=5, length=0, labelrotation=45, bottom=False, top=False, left=False, right=False)
    
    ax["rp"].vlines(0.001, -2, (len(samples)*ss)+ss, linestyle="--", color="grey", linewidth=0.5)
    ax["rp"].vlines(0.005, -2, (len(samples)*ss)+ss, linestyle="--", color="grey", linewidth=0.5)
    ax["rp"].vlines(0.01, -2, (len(samples)*ss)+ss, linestyle="--", color="grey", linewidth=0.5)
    ax["rp"].vlines(0.05, -2, (len(samples)*ss)+ss, linestyle="--", color="grey", linewidth=0.5)

    ## multiple peaks
    # rps, _ = find_peaks(-np.array(recursive_partition["chrom"]))
    # rps = [i for i in rps if recursive_partition["chrom"][i] <= min([0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05], key=lambda x:abs(x-min(recursive_partition["chrom"])))]
    # ax["rp"].hlines([i*ss for i in rps], 0, 1, linestyle="--", linewidth=0.5, color="red")
    # for pos, i in enumerate(recursive_partition["chrom"]):
    #     if i <= min([0.0001, 0.0005, 0.001, 0.005, 0.01, 0.05], key=lambda x:abs(x-min(recursive_partition["chrom"]))) and pos in rps:
    #         ax["rp"].text(0.001, pos*ss, f"{'%.2E' % i}")
    ## minimum peak
    if offy == -1:
        offy = offy*ss
    cmin = np.argmin(recursive_partition["chrom"])
    tmin = np.argmin(recursive_partition["total"])
    if cmin != tmin:
        print(f"optimal chrom split: {telomere_lengths[cmin]}")
        print(f"optimal genom split: {telomere_lengths[tmin]}")
        ax["rp"].hlines([((cmin+1)*ss)+offy], 0, 1, linestyle="--", linewidth=0.5, color="C1")
        ax["rp"].hlines([((tmin+1)*ss)+offy], 0, 1, linestyle="--", linewidth=0.5, color="C0")
        ax["rp"].text(0.91, ((cmin+1)*ss)+offy, conv_mathtext(recursive_partition['chrom'][cmin]), va="bottom", ha="right", fontweight="bold", fontsize=5, color="C1")
        ax["rp"].text(0.91, ((tmin+1)*ss)+offy, conv_mathtext(recursive_partition['total'][tmin]), va="bottom", ha="right", fontweight="bold", fontsize=5, color="C0")
        ax["rp"].text(0.5, (((len(samples)-(len(samples)*0.1))*ss)-(ss/2))/(len(samples)*ss), f"$chrom$ {round(telomere_lengths[cmin], 2)}\n$genome$ {round(telomere_lengths[tmin], 2)}", ha="center", va="top", fontsize=4, fontweight="bold", clip_on=True, transform=ax["rp"].transAxes, zorder=100) 

    else:
        print(f"optimal split: {telomere_lengths[cmin]}")
        ax["rp"].hlines([((tmin+1)*ss)+offy], 0, 1, linestyle="--", linewidth=0.5, color="C0")
        ax["rp"].text(0.91, ((cmin+1)*ss)+offy, conv_mathtext(recursive_partition['chrom'][cmin]), va="top", ha="right", fontsize=5, fontweight="bold", color="C1")
        ax["rp"].text(0.91, ((cmin+1)*ss)+offy, conv_mathtext(recursive_partition['total'][tmin]), va="bottom", ha="right", fontsize=5, fontweight="bold", color="C0")
        ax["rp"].text(0.5, (((len(samples)-1)*ss)-(ss/2))/(len(samples)*ss), f"$opt$ {round(telomere_lengths[cmin], 2)}", ha="center", va="top", fontsize=4, fontweight="bold", clip_on=True, transform=ax["rp"].transAxes, zorder=100) 

    ax["rp"].set(ylabel=None)
    ax["rp"].tick_params(left=False, labelleft=False)
    ax["rp"].grid(False)
    ax2.grid(False)
    ax["rp"].spines[["left", "bottom", "top", "right"]].set_edgecolor("black")
    ax2.spines[["left", "bottom", "top", "right"]].set_edgecolor("black")
    ax["rp"].tick_params(axis="x", which="major", length=3, direction="out")
    #ax["rp"].tick_params(axis="x", which="minor", length=0.5, direction="out")
    
    if cat == False:
        ax["rp"].text(0.5, ((len(samples)*ss)-(ss/2))/(len(samples)*ss), f"$r_s$ {conv_mathtext(stats.spearmanr(telomere_lengths, total_scores).correlation)}\n$p$ {conv_mathtext(stats.spearmanr(telomere_lengths, total_scores).pvalue)}", ha="center", va="top", fontsize=4, fontweight="bold", clip_on=True, transform=ax["rp"].transAxes, zorder=100) 

    if bottom == 1 and ignore_box == False:
        # split = max([cmin, tmin])+1 # split by rp
        score_cols = list(scores.columns)
        split = len([i for i in score_cols if tumour_tl[i] <= threshold])
        chrsplit = pd.DataFrame(columns=["scores", "group"])
        chrsplit = pd.concat([chrsplit, pd.DataFrame.from_dict({"scores": np.array(list(itertools.chain(*scores[scores.columns[:split]].values.tolist()))), "group": "short"})])
        chrsplit = pd.concat([chrsplit, pd.DataFrame.from_dict({"scores": np.array(list(itertools.chain(*scores[scores.columns[split:]].values.tolist()))), "group": "long"})])

        gensplit = pd.DataFrame(columns=["scores", "group"])
        gensplit = pd.concat([gensplit, pd.DataFrame.from_dict({"scores": total_scores[:split], "group": "short"})])
        gensplit = pd.concat([gensplit, pd.DataFrame.from_dict({"scores": total_scores[split:], "group": "long"})])

        print(chrsplit)
        print(gensplit)

        sns.boxplot(ax=ax["chrsplit"], data=chrsplit, x="scores", y="group", hue="group",
                meanprops={'marker' : 'D', 'markeredgecolor' : 'black', 'markersize' : 5, 'label':'mean'},
                #medianprops={'visible': False}, whiskerprops={'visible': False},
                showmeans=True, showfliers=False, showbox=True, showcaps=False, boxprops=dict(alpha=.5), zorder=10)
        sns.stripplot(ax=ax["chrsplit"], data=chrsplit, x="scores", y="group", hue="group", zorder=0)
        sns.boxplot(ax=ax["gensplit"], data=gensplit, x="scores", y="group", hue="group",
                meanprops={'marker' : 'D', 'markeredgecolor' : 'black', 'markersize' : 5, 'label':'mean'},
                #medianprops={'visible': False}, whiskerprops={'visible': False},
                showmeans=True, showfliers=False, showbox=True, showcaps=False, boxprops=dict(alpha=.5), zorder=10)
        sns.stripplot(ax=ax["gensplit"], data=gensplit, x="scores", y="group", hue="group", zorder=0)
        ax["chrsplit"].set_xscale("symlog", base=10)
        ax["gensplit"].set_xscale("symlog", base=10)
        ax["chrsplit"].set_xlim(-0.5, max(chrsplit["scores"])+1)
        ax["gensplit"].set_xlim(min(gensplit["scores"])-0.5, max(gensplit["scores"])+1)
        ax["chrsplit"].set(ylabel=None, xlabel=None, yticks=[])
        ax["gensplit"].set(ylabel=None, xlabel=None, yticks=[])
        chrval = stats.mannwhitneyu(chrsplit[chrsplit['group'] == 'short']['scores'].tolist(), chrsplit[chrsplit['group'] == 'long']['scores'].tolist(), alternative='greater')
        genval = stats.mannwhitneyu(gensplit[gensplit['group'] == 'short']['scores'].tolist(), gensplit[gensplit['group'] == 'long']['scores'].tolist(), alternative='greater')

        ax["chrsplit"].legend(labels=[], labelcolor=["C0", "C1"], loc=2, frameon=False, title=f"chromosome\n$p =$ {conv_mathtext(chrval.pvalue)}", facecolor="none", fontsize=5)
        ax["gensplit"].legend(labels=[], labelcolor=["C0", "C1"], loc=2, frameon=False, title=f"genome\n$p =$ {conv_mathtext(genval.pvalue)}", facecolor="none", fontsize=5)
        ax["chrsplit"].spines[["left", "bottom", "top", "right"]].set_edgecolor("black")
        ax["chrsplit"].tick_params(axis="x", which="major", bottom=True, length=3, direction="out")
        ax["chrsplit"].tick_params(axis="x", which="minor", bottom=True, length=1.5, direction="out")
        ax["gensplit"].spines[["left", "bottom", "top", "right"]].set_edgecolor("black")
        ax["gensplit"].tick_params(axis="x", which="major", bottom=True, length=3, direction="out")
        ax["gensplit"].tick_params(axis="x", which="minor", bottom=True, length=1.5, direction="out")
        #ax["nseggap"].axis("off")
        ax["dummy"].axis("off")
        

    return ax



@cli_tools.command()
@click.option("-p", "--rplot", default=False, is_flag=True, help="Produce copynumber library plots")
@click.option("-s", "--skip", default=False, is_flag=True, help="Skip winsorizing and pcf stages")
@click.option("-t", "--threshold", default=None, help="Fusogenic threshold (heatmap line)")
@click.option("-c", "--chrom-plot", default=False, is_flag=True, help="Plot individual chromosomes")
@click.option("-g", "--skip-genome", default=False, is_flag=True, help="Skip plotting genome CN")
@click.pass_context
def plot(ctx, rplot, skip, threshold, chrom_plot, skip_genome):
    """
    Plot normalised coverage of genomes (3rd)

    Calls an R script that uses a modified version of the copynumber library to allow for more
    reference genomes to be used. Cytoband outside of the library can be used.
    """
    # overwrite if pre commend flag used, allows skipping in whole pipe
    if ctx.obj["skip_genome"] != None:
        skip_genome = True

    mpl.rc('axes.formatter', use_mathtext=True)
    ref = ctx.obj["ref"]
    pd.options.mode.chained_assignment = None

    df = pd.read_csv(os.path.join("cnpinter", ctx.obj["tag"], "copyNumber_io", "copynumbers.all.csv"), low_memory=False)
    cols = list(df.columns)
    samples = cols[2:]
    
    if skip == False:
        print("Running Rscipt...\n")
        if rplot == True:
            pcf(os.getcwd(), ref, samples, ctx.obj["tag"], "yes", ctx.obj["gamma"], ctx.obj["kmin"])
        if rplot == False:
            pcf(os.getcwd(), ref, samples, ctx.obj["tag"], "no", ctx.obj["gamma"], ctx.obj["kmin"])
        print("\nFinished Rscript")
    
    df = pd.read_csv(os.path.join("cnpinter", ctx.obj["tag"], "copyNumber_io", "winsorized.all.csv"), low_memory=False, index_col=0)
    df = df.rename({"chrom": "chromosome"}, axis=1)
    print(df)
    cols = list(df.columns)
    samples = cols[2:]

    segs = pd.read_csv(os.path.join("cnpinter", ctx.obj["tag"], "copyNumber_io", "segmented.copynumber.csv"))    
    segs["sampleID"] = segs["sampleID"].astype(str)
    segs["chrom"] = segs["chrom"].astype(str)
 
    cytoband = True
    if cytoband == True:
        color_lookup = {
                  'gneg': "#ffffff",
                'gpos25': "#999999",
                'gpos50': "#666666",
                'gpos75': "#333333",
               'gpos100': "#000000",
                  'acen': "#cc6666",
                  'gvar': "#cccccc",
                 'stalk': "#e6e6e6"
               }
        if ref.lower() in {"hg19", "grch37"}:
            cyto = pd.read_csv("hg19_cytoBand.txt", sep="\t", header=None)
        elif ref.lower() in {"hg38", "grch38"}:
            cyto = pd.read_csv("hg38_cytoBand.txt", sep="\t", header=None)
        cyto = cyto.rename({0: "c", 1: "s", 2: "e", 3: "a", 4: "t"}, axis=1)
        cyto["w"] = cyto["e"] - cyto["s"]
        cyto["t"] = cyto["t"].replace(color_lookup)

    # for i in samples:
    #     df[i] = winsorize(df[i])
    #     print(i, " ", min(df[i]), " ", max(df[i]))

    chromosomes = sorted([int(i) for i in list(set(df["chromosome"])) if i not in ("X", "Y")])
    chromosomes = [str(i) for i in chromosomes]
    if "X" in set(df["chromosome"]):
        chromosomes.append("X")
    if "Y" in set(df["chromosome"]):
        chromosomes.append("Y")
    #print(chromosomes)

    #chrom_plot = False
    if skip_genome == False:
        whole_genome_plot = True
    else:
        whole_genome_plot = False

    heatmap_plot = True

    squiggle_plot = True
    log_seg = True # False
    squiggle_scale = True # False
    squiggle_perif = True
    squiggle_score_type = "pcf"

    low_pass_plot = True
    low_pass_perif = True
    
    score_box_plot = True

    if not os.path.exists(os.path.join("output", ctx.obj["tag"], "graphs")):
        os.makedirs(os.path.join("output", ctx.obj["tag"], "graphs"))

    for s in samples:
        directory = os.path.join("output", ctx.obj["tag"], "graphs", s)
        if not os.path.exists(directory):
            os.makedirs(directory)
   
    with Progress() as progress:
        ## individual chromosomes
        if chrom_plot == True:
            print("Plotting chromosomes...")
            task_chrom = progress.add_task("Chromosomes: ", total=(len(samples)*len(chromosomes)))
            plt.rcParams["figure.figsize"] = (20,3)
            for s in samples:
                for c in chromosomes:
                    tmp = df[df["chromosome"] == c][["chromosome", "pos", s]]  
                    xmax = max(df[df["chromosome"] == c]["pos"])
                    
                    plt.title(f"{s} chr{c}")
                    plt.grid(False)
                    plt.xlim(0, xmax)
                    plt.ylim(-2.1, 3.1)
                    if cytoband == True:
                        xwidths = list(zip(cyto[cyto["c"] == f"chr{c}"]["s"], cyto[cyto["c"] == f"chr{c}"]["w"]))
                        plt.broken_barh(xwidths, (-2.2, 0.2), facecolors=cyto[cyto["c"] == f"chr{c}"]["t"])
                    
                    plt.hlines(0, 0, xmax, linestyles="dashed", color="black")
                    plt.scatter(tmp["pos"], tmp[s], color="darkgrey", marker="s", linewidths=0, s=5)
                    
                    nover = len(tmp[tmp[s] > 3])
                    nunder = len(tmp[tmp[s] < -2])
                    plt.scatter(tmp[tmp[s] > 3]["pos"], [3]*nover, color="darkgrey", marker="*", linewidths=0, s=10)
                    plt.scatter(tmp[tmp[s] < -2]["pos"], [-2]*nunder, color="darkgrey", marker="*", linewidths=0, s=10)

                    #sampleID	chrom	arm	start.pos	end.pos	n.probes	mean
                    pcf_segs = get_start_end_mean(segs[segs["sampleID"] == s], c)
                    for i in pcf_segs:
                        plt.hlines(i[0], i[1], i[2], color="red", linewidths=3)

                    plt.ylabel("Relative Copy Number")
                    plt.tick_params(top=False, bottom=False, left=True, right=False,
                        labelleft=True, labelbottom=False)
                    plt.tight_layout()
                    plt.savefig(os.path.join("output", ctx.obj["tag"], "graphs", s, f"{s}_chr{c}.pdf"))
                    plt.close()
                    progress.update(task_chrom, advance=1)

        if whole_genome_plot or heatmap_plot:
            task_cumsum = progress.add_task("Cumsum: ", total=len(chromosomes))
            prev = 0
            midpoints = {}
            for c in chromosomes:
                xmax = max(df[df["chromosome"] == c]["pos"])
                midpoints[c] = xmax/2 + prev
                df.loc[df["chromosome"] == c, "pos"] += prev
                segs.loc[segs["chrom"] == c, "start.pos"] += prev
                segs.loc[segs["chrom"] == c, "end.pos"] += prev
                if cytoband == True:
                    cyto.loc[cyto["c"] == f"chr{c}", "s"] += prev
                    cyto.loc[cyto["c"] == f"chr{c}", "e"] += prev
                prev += xmax
                progress.update(task_cumsum, advance=1)


        ## stitch all samples into one
        if whole_genome_plot == True:
            print("Plotting genomes...")
            task_genomes = progress.add_task("Genomes: ", total=len(samples)*2)
            plt.rcParams["figure.figsize"] = (25,5)
            
            ## normalised (with segs)
            for s in samples:
                plt.title(f"{s}")
                plt.grid(False)
                for c in chromosomes:
                    tmp = df[df["chromosome"] == c][["chromosome", "pos", s]]
                    xmax = max(df[df["chromosome"] == c]["pos"])
                    plt.vlines(xmax, -2.1, 3.1, linestyles="dashed", color="black", alpha=0.5)

                    if cytoband == True:
                        xwidths = list(zip(cyto[cyto["c"] == f"chr{c}"]["s"], cyto[cyto["c"] == f"chr{c}"]["w"]))
                        plt.broken_barh(xwidths, (-2.2, 0.2), facecolors=cyto[cyto["c"] == f"chr{c}"]["t"])

                    plt.scatter(tmp["pos"], tmp[s], color="darkgrey", marker="s", linewidths=0, s=1)
                    
                    nover = len(tmp[tmp[s] > 3])
                    nunder = len(tmp[tmp[s] < -2])
                    plt.scatter(tmp[tmp[s] > 3]["pos"], [3]*nover, color="darkgrey", marker="*", linewidths=0, s=4)
                    plt.scatter(tmp[tmp[s] < -2]["pos"], [-2]*nunder, color="darkgrey", marker="*", linewidths=0, s=4)

                    #sampleID	chrom	arm	start.pos	end.pos	n.probes	mean
                    pcf_segs = get_start_end_mean(segs[segs["sampleID"] == s], c)
                    for i in pcf_segs:
                        plt.hlines(i[0], i[1], i[2], color="red", linewidths=3)

                plt.ylim(-2.1, 3.1)
                plt.xlim(0, xmax)
                plt.hlines(0, 0, xmax, linestyle="dashed", color="black", alpha=0.5)
                plt.xticks(ticks=list(midpoints.values()), labels=list(midpoints.keys()))
                plt.xlabel("chromosome")
                plt.ylabel("Relative Copy Number")
                plt.tick_params(top=False, bottom=False, left=True, right=False,
                                labelleft=True, labelbottom=True)
                plt.tight_layout()
                plt.savefig(os.path.join("output", ctx.obj["tag"], "graphs", s, f"{s}_segments.pdf"))
                plt.close()
                progress.update(task_genomes, advance=1)

                ## unclipped
                alternate = 1
                plt.title(f"{s}")
                plt.grid(False)
                for c in chromosomes:
                    tmp = df[df["chromosome"] == c][["chromosome", "pos", s]]
                    xmax = max(df[df["chromosome"] == c]["pos"])
                    #plt.ylim(-2.1, 3.1)
                    plt.vlines(xmax, -2.1, 3.1, linestyles="dashed", color="black", alpha=0.5)
 
                    if alternate % 2 == 0:
                        plt.scatter(tmp["pos"], tmp[s], color="blue", marker="s", linewidths=0, s=1)
                    else:
                        plt.scatter(tmp["pos"], tmp[s], color="orange", marker="s", linewidths=0, s=1)
                    
                    #sampleID	chrom	arm	start.pos	end.pos	n.probes	mean
                    pcf_segs = get_start_end_mean(segs[segs["sampleID"] == s], c)
                    for i in pcf_segs:
                        plt.hlines(i[0], i[1], i[2], color="red", linewidths=3)
                    alternate += 1

                plt.xlim(0, xmax)
                xmin, xmax, ymin, ymax = plt.axis()
                cytoheight = (ymax-ymin)*0.02
                plt.ylim(ymin-cytoheight, ymax)
                if cytoband == True:
                    for c in chromosomes:
                        xwidths = list(zip(cyto[cyto["c"] == f"chr{c}"]["s"], cyto[cyto["c"] == f"chr{c}"]["w"]))
                        plt.broken_barh(xwidths, (ymin-cytoheight, cytoheight), facecolors=cyto[cyto["c"] == f"chr{c}"]["t"])


                plt.hlines(0, 0, xmax, linestyle="dashed", color="black", alpha=0.5)
                plt.xticks(ticks=list(midpoints.values()), labels=list(midpoints.keys()))
                plt.xlabel("chromosome")
                plt.ylabel("Relative Copy Number")
                plt.tick_params(top=False, bottom=False, left=True, right=False,
                                labelleft=True, labelbottom=True)
                plt.tight_layout() 
                plt.savefig(os.path.join("output", ctx.obj["tag"], "graphs", s, f"{s}_win.pdf"))
                plt.close()
                progress.update(task_genomes, advance=1)

        if heatmap_plot or squiggle_plot or low_pass_plot:
            tot_segs = segs.copy()
            backgrounds, tumour_tl = get_backgrounds(ctx.obj["bgfile"], ctx.obj["bgtc"], ctx.obj["bgts"], ctx.obj["bgnc"])
            nsamples = len(set(list(segs["sampleID"]))) 
            if ctx.obj["cat"] == False:
                shorts = [k for k, v in tumour_tl.items() if v <= ctx.obj["thresh"]]
                longs = [k for k, v in tumour_tl.items() if v > ctx.obj["thresh"]]
            else:
                shorts = [k for k, v in tumour_tl.items() if v == True]
                longs = [k for k, v in tumour_tl.items() if v == False]
            all_samples = list(set(segs["sampleID"]))
            all_samples = list({s: l for s, l in sorted(tumour_tl.items(), key=lambda item: item[1]) if s in all_samples}.keys())
            shorts = [i for i in shorts if i in all_samples]
            longs = [i for i in longs if i in all_samples]
            short_set = []
            long_set = []
            ytick_labels_set = []
            nsubsets = 0
            for i in range(ctx.obj["ssets"]):
                if nsamples > ctx.obj["sub"]:
                    if ctx.obj["cat"] == False:
                        short_samples = [shorts[i] for i in sorted(random.sample(range(len(shorts)), min(len(shorts), int(ctx.obj["sub"]*ctx.obj["split"]))))]
                        long_samples = [longs[i] for i in sorted(random.sample(range(len(longs)), min(len(longs), int(ctx.obj["sub"]*(1-ctx.obj["split"])))))]
                    else:
                        short_samples = [shorts[i] for i in sorted(random.sample(range(len(shorts)), min(len(shorts), ctx.obj["sub"]*ctx.obj["split"])))]
                        long_samples = [longs[i] for i in sorted(random.sample(range(len(longs)), min(len(longs), ctx.obj["sub"]*(1-ctx.obj["split"]))))]
                    short_set.append(short_samples)
                    long_set.append(long_samples)
                    shorts = list(set(shorts).difference(short_samples))
                    longs = list(set(longs).difference(long_samples))
                    nsamples = nsamples - len(short_samples) - len(long_samples)
                elif 0 < nsamples < ctx.obj["sub"]:
                    short_set.append(shorts)
                    long_set.append(longs)
                    nsamples = 0

                subset_samples = dict((k, tumour_tl[k]) for k in short_set[i] + long_set[i])
                subset_samples = list(dict(sorted(subset_samples.items(), key=lambda k: k[1])).keys())
                if ctx.obj["cat"] == False:
                    ytick_labels = [f"{s.replace(ctx.obj['prefix'], '')[:ctx.obj['namecrop']]} - {to_sigfig(tumour_tl[s], ctx.obj['sigfig'])}" for s in subset_samples]
                else:
                    ytick_labels = [f"{s.replace(ctx.obj['prefix'], '')[:ctx.obj['namecrop']]} - {bool_to_label(tumour_tl[s])}" for s in subset_samples]
                ytick_labels_set.append(ytick_labels)

                nsubsets += 1

                if nsamples == 0:
                    break

        for z in range(nsubsets):
            ytick_labels = ytick_labels_set[z]
            samples = short_set[z] + long_set[z]
            samples = [f"{s}" for s in samples]
            samples = list({s: l for s, l in sorted(tumour_tl.items(), key=lambda item: item[1]) if s in samples}.keys())
            print("samples: ", samples)
            segs = tot_segs[tot_segs["sampleID"].isin(samples)]
            print("segs: ", segs)
            ## heatmap all samples
            if heatmap_plot == True:
                print("Plotting heatmap...")
                plt.rcParams["figure.figsize"] = (20, 20) #plt.rcParamsDefault["figure.figsize"]
                segs["width"] = segs["end.pos"] - segs["start.pos"]
                segs["mean_clip"] = segs["mean"].clip(-1.5, 1.5)

                norm = mpl.colors.Normalize(vmin=-1.5, vmax=1.5)
                cmap = mpl.colormaps.get_cmap("RdBu_r")  # PiYG

                print(tumour_tl)
                # samples = [f"{s}" for s in samples]
                # samples = list({s: l for s, l in sorted(tumour_tl.items(), key=lambda item: item[1]) if s in samples}.keys())
                # ytick_labels = [f"{s[:5]} - {'%.3f' % tumour_tl[s]}" for s in samples]
                print(samples)
                print(ytick_labels)

                for y, s in enumerate(samples):
                    xwidths = list(segs[segs["sampleID"] == s][["start.pos", "width"]].itertuples(index=False, name=None))
                    plt.broken_barh(xwidths, (y, 1), facecolors=cmap(norm(segs[segs["sampleID"] == s]["mean_clip"])))
                
                for idx, c in enumerate(chromosomes):
                    xmax = max(df[df["chromosome"] == c]["pos"])
                    if idx != len(chromosomes)-1:
                        plt.vlines(xmax, 0, len(samples), linestyles="dashed", color="black", alpha=0.5)

                if threshold != None:
                    thresh_line = len([tumour_tl[i] for i in samples if float(tumour_tl[i]) <= float(threshold)])
                    print(thresh_line)
                    plt.axhline(y=thresh_line, color="black", linestyle="dashed")

                plt.grid(False)
                plt.xlim(0, xmax) 
                plt.ylim(0, len(samples))
                plt.xticks(ticks=list(midpoints.values()), labels=list(midpoints.keys()))
                plt.yticks(ticks=[i-0.5 for i in range(1, len(samples)+1)], labels=ytick_labels)
                plt.tick_params(top=False, bottom=False, left=True, right=False,
                                    labelleft=True, labelbottom=True)
                plt.xlabel("chromosome")
                plt.ylabel("sample")
                plt.tight_layout()
                plt.savefig(os.path.join("output", ctx.obj["tag"], f"heatmap_{z}.pdf"))
                plt.close()


            ## squiggle plot
            if squiggle_plot == True:
                ref = ctx.obj["ref"]
                if ref.lower() in {"hg38", "grch38"}:
                    chrom_lengths_dict={'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                                    'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422, 
                                    'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 
                                    'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167, 
                                    'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415} 
                elif ref.lower() in {"hg19", "grch37"}:
                    chrom_lengths_dict={'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260, 
                                    'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747, 
                                    'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 
                                    'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520, 
                                    'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566}
                print("Plotting squiggle...")
                # task_squiggle = progress.add_task("Squiggles: ", total=len(samples))
                #plt.rcParams["figure.figsize"] = (25,5)
                # if cytoband == True:
                #     xwidths = list(zip(cyto[cyto["c"] == f"chr{c}"]["s"], cyto[cyto["c"] == f"chr{c}"]["w"]))
                #     plt.broken_barh(xwidths, (-2.2, 0.2), facecolors=cyto[cyto["c"] == f"chr{c}"]["t"])
                # backgrounds, tumour_tl = get_backgrounds(ctx.obj["bgfile"], ctx.obj["bgtc"], ctx.obj["bgts"], ctx.obj["bgnc"])
                # ytick_labels = [f"{s[:5]} - {'%.3f' % tumour_tl[s]}" for s in samples]
                #plt.grid(False)
                uc = 3
                lc = 2
                tc = lc+uc
                hc = uc/tc
                dc = lc/tc
                nsegs = []
                greyout = 0.05

                boxB = ctx.obj["ignore_box"]
                scoreB = ctx.obj["ignore_score"]
                scaleB = ctx.obj["yaxis_scale"]

                fig = plt.figure(figsize=(8.27,11.69), constrained_layout=True)
                # set axes depending on ignore flags
                if squiggle_perif == False:
                    ax = fig.subplot_mosaic([["squig"]])
                else:
                    if ctx.obj["ssets"] == 1 and boxB == False:
                        if scoreB == False and scaleB == False:
                            ax = fig.subplot_mosaic([["squig", "squig", "score", "rp"], ["chrsplit", "gensplit", "dummy", "dummy"]], gridspec_kw={"width_ratios": [10, 10, 2, 2], "height_ratios": [20, 1]})
                        elif scoreB == True and scaleB == False:
                            ax = fig.subplot_mosaic([["squig", "squig", "rp"], ["chrsplit", "gensplit", "dummy"]], gridspec_kw={"width_ratios": [10, 10, 2], "height_ratios": [20, 1]})
                        elif scoreB == False and scaleB == True:
                            ax = fig.subplot_mosaic([["scale", "squig", "squig", "score", "rp"], ["dummy_scale", "chrsplit", "gensplit", "dummy", "dummy"]], gridspec_kw={"width_ratios": [1, 10, 10, 2, 2], "height_ratios": [20, 1]})
                            ax["dummy_scale"].axis("off")
                        elif scoreB == True and scaleB == True:
                            ax = fig.subplot_mosaic([["scale", "squig", "squig", "rp"], ["dummy_scale", "chrsplit", "gensplit", "dummy"]], gridspec_kw={"width_ratios": [1, 10, 10, 2], "height_ratios": [20, 1]})
                            ax["dummy_scale"].axis("off")
                    elif ctx.obj["ssets"] > 1 or boxB == True:
                        if scoreB == False and scaleB == False:
                            ax = fig.subplot_mosaic([["squig", "squig", "score", "rp"]], gridspec_kw={"width_ratios": [10, 10, 2, 2]})
                        elif scoreB == True and scaleB == False:
                            ax = fig.subplot_mosaic([["squig", "squig", "rp"]], gridspec_kw={"width_ratios": [10, 10, 2]})
                        elif scoreB == False and scaleB == True:
                            ax = fig.subplot_mosaic([["scale", "squig", "squig", "score", "rp"]], gridspec_kw={"width_ratios": [1, 10, 10, 2, 2]})
                        elif scoreB == True and scaleB == True:
                            ax = fig.subplot_mosaic([["scale", "squig", "squig", "rp"]], gridspec_kw={"width_ratios": [1, 10, 10, 2]})

                total_genome_size = sum(chrom_lengths_dict.values())

                if log_seg == True:
                    log_x = {}
                    for ypos, s in enumerate(samples):
                        log_x[s] = {}
                        nsegs.append(len(segs[(segs["sampleID"] == s) & (abs(segs["mean"]) > greyout)]))
                        for idx, c in enumerate(chromosomes):
                            log_x[s][c] = []
                            #sampleID	chrom	arm	start.pos	end.pos	n.probes	mean
                            pcf_segs = get_start_end_mean(segs[segs["sampleID"] == s], c)
                            for i in pcf_segs:
                                if i[2]-i[1] == 0:
                                    log_x[s][c].append(1)
                                else:
                                    log_x[s][c].append(np.log10(i[2]-i[1]))
                else:
                    log_x = {}
                    for ypos, s in enumerate(samples):
                        log_x[s] = {}
                        nsegs.append(len(segs[(segs["sampleID"] == s) & (abs(segs["mean"]) > greyout)]))
                        for idx, c in enumerate(chromosomes):
                            log_x[s][c] = []
                            #sampleID	chrom	arm	start.pos	end.pos	n.probes	mean
                            pcf_segs = get_start_end_mean(segs[segs["sampleID"] == s], c)
                            for i in pcf_segs:
                                if i[2]-i[1] == 0:
                                    log_x[s][c].append(1)
                                else:
                                    log_x[s][c].append(i[2]-i[1])

                # norm = mpl.colors.TwoSlopeNorm(vmin=-2, vmax=10, vcenter=0)
                # cmap = cm.get_cmap("RdBu")
                unorm = mpl.colors.Normalize(vmin=0.0, vmax=10.0)
                #ucmap = mpl.colormaps.get_cmap("autumn")  # PiYG
                ucmap = plt.get_cmap("autumn")
                lnorm = mpl.colors.Normalize(vmin=0.0, vmax=2.0)
                #lcmap = mpl.colormaps.get_cmap("winter")  # PiYG
                lcmap = plt.get_cmap("winter")
                scores = {}
                for ypos, s in enumerate(samples):
                    nsegs.append(len(segs[(segs["sampleID"] == s) & (abs(segs["mean"]) > greyout)]))
                    scores[s] = {}
                    xtick_locs = []
                    if squiggle_scale == False:
                        relative_idx = 0
                        relative_size = 0
                    for idx, c in enumerate(chromosomes):
                        xvals = []
                        yvals = []
                        means = []
                        widths = []
                        #sampleID	chrom	arm	start.pos	end.pos	n.probes	mean
                        pcf_segs = get_start_end_mean(segs[segs["sampleID"] == s], c)
                        for i in pcf_segs:
                            means.append(i[0])
                            if log_seg == False:
                                xvals.append((i[1]/chrom_lengths_dict[f"chr{c}"])+idx)
                                #xvals.append((i[2]/chrom_lengths_dict[f"chr{c}"])+idx)
                                widths.append((i[2]-i[1])/chrom_lengths_dict[f"chr{c}"])
                            if i[0] > 0:
                                if i[0] > uc: # clip
                                    yvals.append(hc)
                                    #yvals.append(ypos+1)
                                else:
                                    yvals.append((i[0]/tc))
                                    #yvals.append((i[0]/tc)+dc+ypos)
                            elif i[0] < 0:
                                if i[0] < -lc: # clip
                                    yvals.append(-dc)
                                    #yvals.append(ypos)
                                else:
                                    yvals.append((i[0]/tc))
                                    #yvals.append((i[0]/tc)+dc+ypos)
                            else:
                                yvals.append(0)
                                #yvals.append(ypos+dc)
                        
                        #if log_seg == True:
                        ## plot line
                        log_len = sum(log_x[s][c])
                        ## bar plot
                        if squiggle_scale == True:
                            xvals = [idx] + [(i/log_len)+idx for i in list(np.cumsum(log_x[s][c][:-1]))]
                            widths = [i/log_len for i in log_x[s][c]]
                        else:
                            relative_size = chrom_lengths_dict[f"chr{c}"]/total_genome_size
                            xvals = [relative_idx] + [((i/log_len)*relative_size)+relative_idx for i in list(np.cumsum(log_x[s][c][:-1]))]
                            widths = [(i/log_len)*relative_size for i in log_x[s][c]]
                            idx = relative_idx
                            xtick_locs.append(relative_idx+(relative_size/2))
                            relative_idx += relative_size
                            # (0.0249, 25000, 92075000)
                            # (0.0275, 94185000, 242095000)

                        ## seperating chromosomes lines
                        ax["squig"].vlines(idx, 0, len(samples), color="black", alpha=0.4, linewidths=0.01, linestyle='--')

                        #ax["squig"].plot(xvals, yvals, '.-', color="black", linewidth=0.2, markersize=1)
                        mid = [dc+ypos]*len(yvals)
                        upp = [dc+(greyout/tc)+ypos]*len(yvals)
                        lww = [dc-(greyout)+ypos]*len(yvals)
                        xvals = np.array(xvals)
                        yvals = np.array(yvals)
                        mid = np.array(mid)
                        #ax["squig"].plot(xvals, mid, '--', alpha=0.3)
                   
                        #print(len(xvals), len(yvals), len(mid))
                        if len(pcf_segs) > 0:
                            #ax["squig"].bar(xvals, yvals, width=widths, color=cmap(norm([-1 for i in means])), alpha=0.9, align="edge")
                            ax["squig"].bar([i for i, x in zip(xvals, yvals) if x > 0], [i for i in yvals if i > 0], width=[i for i, x in zip(widths, yvals) if x > 0], bottom=dc+ypos, color=ucmap(unorm([i for i in [i for i in means if i > 0]])), alpha=1, linewidth=0, align="edge", zorder=50)
                            ax["squig"].bar([i for i, x in zip(xvals, yvals) if x < 0], [i for i in yvals if i < 0], width=[i for i, x in zip(widths, yvals) if x < 0], bottom=dc+ypos, color=lcmap(lnorm([-i for i in [i for i in means if i < 0]])), alpha=1, linewidth=0, align="edge", zorder=50)

                        # ax["squig"].fill_between(xvals, yvals, mid, where=(yvals <= lww), color='blue', alpha=0.5, linestyle="None")
                        # ax["squig"].fill_between(xvals, yvals, mid, where=(yvals >= upp), color='red', alpha=0.5, linestyle="None")
                        # # ax["squig"].fill_between(xvals, yvals, mid, where=(yvals == ypos+1), color=cmap(norm(rmeans)), alpha=0.5, linestyle="None")
                        # ax["squig"].fill_between(xvals, yvals, mid, where=(yvals < upp), color='black', alpha=0.1, linestyle="None")
                        # ax["squig"].fill_between(xvals, yvals, mid, where=(yvals > lww), color='black', alpha=0.1, linestyle="None")
                        means = [abs(i) for i in means]
                        scores[s][c] = sum(means)
                        if squiggle_scale == True:
                            if 1 <= sum(means)/10 < 3:
                                #ax["squig"].text(idx+0.5, ypos-0.05, "*", fontsize="medium", color="black", fontweight="bold", ha="center")
                                ax["squig"].add_patch(mpl.patches.Rectangle((idx, ypos), 1, 1, color="darkgrey", alpha=0.25, zorder=1, linewidth=0, edgecolor=None))
                            elif 3 <= sum(means)/10 < 5:
                                #ax["squig"].text(idx+0.5, ypos-0.05, "**", fontsize="medium", color="black", fontweight="bold", ha="center")
                                ax["squig"].add_patch(mpl.patches.Rectangle((idx, ypos), 1, 1, color="darkgrey", alpha=0.5, zorder=1, linewidth=0, edgecolor=None))
                            elif 5 <= sum(means)/10:
                                #ax["squig"].text(idx+0.5, ypos-0.05, "***", fontsize="medium", color="black", fontweight="bold", ha="center")
                                ax["squig"].add_patch(mpl.patches.Rectangle((idx, ypos), 1, 1, color="darkgrey", alpha=0.75, zorder=1, linewidth=0, edgecolor=None))
                        else:
                            if 1 <= sum(means)/10 < 3:
                                ax["squig"].text(xtick_locs[-1], ypos-0.05, "*", fontsize="medium", color="black", fontweight="bold", ha="center")
                            elif 3 <= sum(means)/10 < 5:
                                ax["squig"].text(xtick_locs[-1], ypos-0.05, "**", fontsize="medium", color="black", fontweight="bold", ha="center")
                            elif 5 <= sum(means)/10:
                                ax["squig"].text(xtick_locs[-1], ypos-0.05, "***", fontsize="medium", color="black", fontweight="bold", ha="center")

                    # ax["squig"].text(idx+1.1, ypos+dc, nsegs[ypos], fontsize="medium")
                # ax["squig"].text(idx+1.1, len(samples)+dc, "nsegs", fontsize="medium")
                # ax["squig"].colorbar(ucmap, label="gains")
                # ax["squig"].colorbar(lcmap, label="loss")
                ax["squig"].set_ylim(0, len(samples))
                if squiggle_scale == True:
                    ax["squig"].set_xlim(0, len(chromosomes))
                    ax["squig"].set_xticks(ticks=[i+0.5 for i in range(len(chromosomes))], labels=chromosomes)
                else:
                    ax["squig"].set_xlim(0, 1)
                    ax["squig"].set_xticks(ticks=xtick_locs, labels=chromosomes)
                    ax["squig"].set_xticklabels(chromosomes, rotation = 45, ha="center")
                
                if scaleB == False:
                    ax["squig"].set_yticks(ticks=[i-hc for i in range(1, len(samples)+1)], labels=ytick_labels)
                else:
                    ax["squig"].set_yticks([])
                    scale_l = [to_sigfig(tumour_tl[s], ctx.obj['sigfig']) for s in subset_samples] 
                    scale_x = [0]*len(scale_l)
                    scale_y = [i-hc for i in range(1, len(subset_samples)+1)]
                    scale_heat = cm.jet_r
                    scale_norm = mpl.colors.Normalize(vmin=min(scale_l), vmax=max(scale_l))
                    scale_map = cm.ScalarMappable(norm=scale_norm, cmap=scale_heat)
                    scale_col = [scale_map.to_rgba(i) for i in scale_l]

                    ax["scale"].scatter(scale_x, scale_y, c = scale_col, edgecolor="none", marker="s")
                    ax["scale"].spines[["top", "bottom", "left", "right"]].set_visible(False)
                    ax["scale"].set_ylabel("Telomer length")
                    ax["scale"].set_xticks([])
                    ax["scale"].set_yticks([])
                    ax["scale"].text(0, 0, f"{round(min(scale_l), 2)} - ", ha="right", va="center", zorder=100)
                    ax["scale"].text(0, max(scale_y), f"{round(max(scale_l), 2)} - ", ha="right", va="center")
                    for i in range(int(min(scale_l))+1, int(max(scale_l))):
                        ax["scale"].text(0, np.argmin(np.abs(np.array(scale_l)-i))+1-hc, f"{i} - ", ha="right", va="center")

                ax["squig"].set_xlabel("chromosome")
                ax["squig"].set_ylabel("Relative Copy Number")
                ax["squig"].tick_params(top=False, bottom=False, left=True, right=False,
                                labelleft=True, labelbottom=True)
                for spine in ["top", "bottom", "left", "right"]:
                    ax["squig"].spines[spine].set_visible(False)
                if squiggle_score_type == "lowpass":
                    scores = lowpass_score(df[df["sampleID"].isin(sampels)])
                print("scores", len(scores))
                if squiggle_perif == True:
                    ax = add_periferal_stats(ax, scores, samples, tumour_tl, ctx.obj["cat"], ctx.obj["thresh"], ss=1, boty=0, offy=-hc, bottom=ctx.obj["ssets"], ignore_score=scoreB, ignore_box=boxB) 
                plt.grid(False)
                for i in ax:
                    ax[i].grid(False)
                fig.tight_layout()
                plt.subplots_adjust(wspace=0.1, hspace=0.1)
                plt.savefig(os.path.join("output", ctx.obj["tag"], f"squiggle_{z}.pdf"))
                plt.close()
                # progress.update(task_squiggle, advance=1)
                # plt.show()

            if low_pass_plot == True:
                print("Plotting ridge low pass...")
                T = 5.0         # Sample Period
                fs = 30.0       # sample rate, Hz
                cutoff = 0.5    # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2Hz
                nyq = 0.5 * fs  # Nyquist Frequency
                order = 2       # sin wave can be approx represented as quadratic
                n = int(T * fs) # total number of samples

                def butter_lowpass_filter(data, cutoff, fs, order):
                    normal_cutoff = cutoff / nyq
                    # Get the filter coefficients
                    b, a = butter(order, normal_cutoff, btype='low', analog=False)
                    y = filtfilt(b, a, data, axis=0)
                    return y

                print(tumour_tl)
                # samples = list(df.columns)[2:]
                # ytick_labels = [f"{s[:5]} - {'%.3f' % tumour_tl[s]}" for s in samples]

                ss = 4
                fig = plt.figure(figsize=(8.27,11.69), constrained_layout=True)
                #ax = fig.subplot_mosaic([["ridge", "score", "rp"]], gridspec_kw={"width_ratios": [20, 2, 1]})
                if low_pass_perif == False:
                    ax = fig.subplot_mosaic([["ridge"]])
                elif low_pass_perif == True and ctx.obj["ssets"] == 1:
                    ax = fig.subplot_mosaic([["ridge", "ridge", "score", "rp"], ["chrsplit", "gensplit", "dummy", "dummy"]], gridspec_kw={"width_ratios": [10, 10, 2, 2], "height_ratios": [20, 1]})
                else:
                    ax = fig.subplot_mosaic([["ridge", "ridge", "score", "rp"]], gridspec_kw={"width_ratios": [10, 10, 2, 2]}) 

                # samples = list(df.columns)[2:]
                scores = {}
                for ypos, sample in enumerate(samples):
                    scores[sample] = {}
                    data = df[["chromosome", sample]]
                    data = data.rename({"chromosome": "chrom"}, axis=1)
                    chrom_lines = list(range(1,23))
                    #data = np.array(data[sample].tolist())
                    l = len(data)
                    #y = butter_lowpass_filter(data, cutoff, fs, order)
                    chroms = [str(i) for i in range(1, 23)] + ["X"]
                    for idx, chrom in enumerate(chroms):
                        filt_data = np.array(data[data["chrom"] == chrom][sample].tolist())
                        ## low pass
                        y = butter_lowpass_filter(filt_data, cutoff, fs, order)
                        p, properties = find_peaks(y)
                        t, tproperties = find_peaks(-y)

                        c = []
                        d = []
                        if len(p) > len(t):
                            c = [i for x in zip(p, t) for i in x] + [p[-1]]
                        else:
                            c = [i for x in zip(t, p) for i in x] + [t[-1]]
                      
                        c = list(sorted(c))
                        for didx, i in enumerate(c[:-1]):
                            d.append(abs(y[c[didx]] - y[c[didx+1]]))


                        l = len(y)
                        # ax.plot(np.linspace(idx,idx+1,l), y, color="none")
                        ax["ridge"].fill_between(np.linspace(idx,idx+1,l), [y+(ypos*ss) for y in y], ypos*ss, where=(y>0), color="red", alpha=0.5, zorder=ypos)
                        ax["ridge"].fill_between(np.linspace(idx,idx+1,l), [y+(ypos*ss) for y in y], ypos*ss, where=(y<0), color="blue", alpha=0.5, zorder=ypos)
                        scores[sample][chrom] = int(round(sum([i for i in d if i > 0.5]), 0))
                        if int(round(sum([i for i in d if i > 0.5]), 0)) > 3:
                            ax["ridge"].text(idx+0.5, (ypos*ss)-(ss/2), f"{int(round(sum([i for i in d if i > 0.5]), 0))}", ha="center")
                ax["ridge"].vlines(chrom_lines, -2, (ypos*ss)+5, color="black", linestyle="--", linewidth=0.5)
                ax["ridge"].set_xlim(0, 23)
                ax["ridge"].set_ylim(-2, (ypos*ss)+ss)
                ax["ridge"].set_xticks([i+0.5 for i in range(23)], chroms)
                ax["ridge"].set_yticks([i*ss for i in range(len(samples))], ytick_labels)#[s[:5] for s in samples])
                ax["ridge"].set_xlabel("chromosome")
                ax["ridge"].set_ylabel("relative copy number")
                ax["ridge"].set_title("Low Pass filtered coverage")
                ax["ridge"].spines[["right", "top", "left", "bottom"]].set_visible(False)
                ax["ridge"].tick_params(length=0)
                plt.tight_layout()
                # plt.show()
                plt.subplots_adjust(wspace=0.1, hspace=0.1)
                if low_pass_perif:
                    ax = add_periferal_stats(ax, scores, samples, tumour_tl, ctx.obj["cat"], ctx.obj["thresh"], ss, bottom=ctx.obj["ssets"])
                # plt.show()
                plt.grid(False)
                for i in ax:
                    ax[i].grid(False)
                plt.savefig(os.path.join("output", ctx.obj["tag"], f"low_pass_{z}.pdf"))
                plt.close()

        if ctx.obj["ssets"] > 1 or score_box_plot == True:
            # regather all samples
            if ctx.obj["cat"] == False:
                shorts = [k for k, v in tumour_tl.items() if v <= ctx.obj["thresh"]]
                longs = [k for k, v in tumour_tl.items() if v > ctx.obj["thresh"]]
            else:
                shorts = [k for k, v in tumour_tl.items() if v == True]
                longs = [k for k, v in tumour_tl.items() if v == False]

            # sort by length
            samples = list(set(segs["sampleID"]))
            samples = list({s: l for s, l in sorted(tumour_tl.items(), key=lambda item: item[1]) if s in all_samples}.keys())

            shorts = [i for i in shorts if i in samples]
            longs = [i for i in longs if i in samples]

            cn_df = pd.read_csv(os.path.join("cnpinter", ctx.obj["tag"], "copyNumber_io", "copynumbers.all.csv"), low_memory=False)
            samples = [i for i in samples if i in set(cn_df[2:].columns)]
            split = len([i for i in samples if tumour_tl[i] <= ctx.obj["thresh"]])

            low_scores = lowpass_score(cn_df)

            pcf_scores = {}
            for ypos, s in enumerate(samples):
                nsegs.append(len(tot_segs[(tot_segs["sampleID"] == s) & (abs(tot_segs["mean"]) > greyout)]))
                pcf_scores[s] = {}
                for idx, c in enumerate(chromosomes):
                    means = []
                    #sampleID	chrom	arm	start.pos	end.pos	n.probes	mean
                    pcf_segs = get_start_end_mean(tot_segs[tot_segs["sampleID"] == s], c)
                    for i in pcf_segs:
                        means.append(abs(i[0]))
                    pcf_scores[s][c] = sum(means)

            pcf_scores = pd.DataFrame.from_dict(pcf_scores)
            low_scores = pd.DataFrame.from_dict(low_scores)
            print("pcf scores")
            print(pcf_scores)
            print("low scores")
            print(low_scores)

            pcf_total_scores = []
            low_total_scores = []
            for ypos, sample in enumerate(samples):
                pcf_total_scores.append(np.sum(pcf_scores[sample]))
                low_total_scores.append(np.sum(low_scores[sample]))

            pcf_chrsplit = pd.DataFrame(columns=["scores", "group"])
            pcf_chrsplit = pd.concat([pcf_chrsplit, pd.DataFrame.from_dict({"scores": np.array(list(itertools.chain(*pcf_scores[shorts].values.tolist()))), "group": "short"})])
            pcf_chrsplit = pd.concat([pcf_chrsplit, pd.DataFrame.from_dict({"scores": np.array(list(itertools.chain(*pcf_scores[longs].values.tolist()))), "group": "long"})])

            pcf_gensplit = pd.DataFrame(columns=["scores", "group"])
            pcf_gensplit = pd.concat([pcf_gensplit, pd.DataFrame.from_dict({"scores": pcf_total_scores[:split], "group": "short"})])
            pcf_gensplit = pd.concat([pcf_gensplit, pd.DataFrame.from_dict({"scores": pcf_total_scores[split:], "group": "long"})])

            low_chrsplit = pd.DataFrame(columns=["scores", "group"])
            low_chrsplit = pd.concat([low_chrsplit, pd.DataFrame.from_dict({"scores": np.array(list(itertools.chain(*low_scores[shorts].values.tolist()))), "group": "short"})])
            low_chrsplit = pd.concat([low_chrsplit, pd.DataFrame.from_dict({"scores": np.array(list(itertools.chain(*low_scores[longs].values.tolist()))), "group": "long"})])

            low_gensplit = pd.DataFrame(columns=["scores", "group"])
            low_gensplit = pd.concat([low_gensplit, pd.DataFrame.from_dict({"scores": low_total_scores[:split], "group": "short"})])
            low_gensplit = pd.concat([low_gensplit, pd.DataFrame.from_dict({"scores": low_total_scores[split:], "group": "long"})])

            pcf_chrval = stats.mannwhitneyu(pcf_chrsplit[pcf_chrsplit['group'] == 'short']['scores'].tolist(), pcf_chrsplit[pcf_chrsplit['group'] == 'long']['scores'].tolist(), alternative='greater')
            pcf_genval = stats.mannwhitneyu(pcf_gensplit[pcf_gensplit['group'] == 'short']['scores'].tolist(), pcf_gensplit[pcf_gensplit['group'] == 'long']['scores'].tolist(), alternative='greater')
            low_chrval = stats.mannwhitneyu(low_chrsplit[low_chrsplit['group'] == 'short']['scores'].tolist(), low_chrsplit[low_chrsplit['group'] == 'long']['scores'].tolist(), alternative='greater')
            low_genval = stats.mannwhitneyu(low_gensplit[low_gensplit['group'] == 'short']['scores'].tolist(), low_gensplit[low_gensplit['group'] == 'long']['scores'].tolist(), alternative='greater')

            # pcf_chrsplit["scores"].loc[pcf_chrsplit["scores"] == 0] = 0.01
            # pcf_gensplit["scores"].loc[pcf_gensplit["scores"] == 0] = 0.01
            # low_chrsplit["scores"].loc[low_chrsplit["scores"] == 0] = 0.01
            # low_gensplit["scores"].loc[low_gensplit["scores"] == 0] = 0.01

            # fig = plt.figure(figsize=(10,4), constrained_layout=True)
            # ax = fig.subplot_mosaic([["pcs", "pts"], ["lcs", "lts"]], gridspec_kw={"width_ratios": [10, 10], "height_ratios": [10, 10]})
            # fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(10, 4), constrained_layout=True)
            fig = plt.figure(figsize=(10,4), constrained_layout=True)
            fig.suptitle("Complexity scores")

            subfigs = fig.subfigures(nrows=2, ncols=1)
            for subfig, row_title in zip(subfigs, ["PCF", "Low Pass"]):
                subfig.suptitle(row_title)
 
            
            short_patch = mpl.patches.Patch(color="C0", label="Short")
            long_patch = mpl.patches.Patch(color="C1", label="Long")

            ax = subfigs[0].subplots(nrows=1, ncols=2)

            ax[0].set_title("Chromosome")
            ax[1].set_title("Genome")

            sns.boxplot(ax=ax[0], data=pcf_chrsplit, x="scores", y="group", hue="group",
                    meanprops={'marker' : 'D', 'markeredgecolor' : 'black', 'markersize' : 5, 'label':'mean'},
                    #medianprops={'visible': False}, whiskerprops={'visible': False},
                    showmeans=True, showfliers=False, showbox=True, showcaps=False, boxprops=dict(alpha=.5), zorder=10)
            sns.stripplot(ax=ax[0], data=pcf_chrsplit, x="scores", y="group", hue="group", zorder=0, s=2, jitter=0.35)

            sns.boxplot(ax=ax[1], data=pcf_gensplit, x="scores", y="group", hue="group",
                    meanprops={'marker' : 'D', 'markeredgecolor' : 'black', 'markersize' : 5, 'label':'mean'},
                    #medianprops={'visible': False}, whiskerprops={'visible': False},
                    showmeans=True, showfliers=False, showbox=True, showcaps=False, boxprops=dict(alpha=.5), zorder=10)
            sns.stripplot(ax=ax[1], data=pcf_gensplit, x="scores", y="group", hue="group", zorder=0, s=2, jitter=0.35)
            ax[0].set_xscale("log", base=10)
            ax[1].set_xscale("log", base=10)
            ax[0].set_xlim(0, max(pcf_chrsplit["scores"])+25)
            ax[1].set_xlim(0, max(pcf_gensplit["scores"])+25)
            ax[0].set(ylabel=None, xlabel=None, yticks=[])
            ax[1].set(ylabel=None, xlabel=None, yticks=[])

            # ax["pcs"].legend(labels=["short", "long"], labelcolor=["C0", "C1"], loc=2, frameon=False, title=f"chromosome\n$p =$ {conv_mathtext(chrval.pvalue)}", facecolor="none", fontsize=5)
            # ax["pts"].legend(labels=["short", "long"], labelcolor=["C0", "C1"], loc=2, frameon=False, title=f"genome\n$p =$ {conv_mathtext(genval.pvalue)}", facecolor="none", fontsize=5)
            
            ax[0].spines[["right", "top"]].set_visible(False)
            ax[0].spines[["left", "bottom"]].set_edgecolor("black")
            ax[1].spines[["right", "top"]].set_visible(False)
            ax[1].spines[["left", "bottom"]].set_edgecolor("black")

            ax[0].legend(labels=[], frameon=False, title=f"$p =$ {conv_mathtext(pcf_chrval.pvalue)}", facecolor="none", fontsize=5)
            ax[1].legend(labels=[], frameon=False, title=f"$p =$ {conv_mathtext(pcf_genval.pvalue)}", facecolor="none", fontsize=5)

            ax[0].grid(False)
            ax[0].tick_params(axis="x", which="major", bottom=True, length=3)
            ax[0].tick_params(axis="x", which="minor", bottom=True, length=1.5)
            ax[1].grid(False)
            ax[1].tick_params(axis="x", which="major", bottom=True, length=3)
            ax[1].tick_params(axis="x", which="minor", bottom=True, length=1.5)

            ax = subfigs[1].subplots(nrows=1, ncols=2)

            sns.boxplot(ax=ax[0], data=low_chrsplit, x="scores", y="group", hue="group",
                    meanprops={'marker' : 'D', 'markeredgecolor' : 'black', 'markersize' : 5, 'label':'mean'},
                    #medianprops={'visible': False}, whiskerprops={'visible': False},
                    showmeans=True, showfliers=False, showbox=True, showcaps=False, boxprops=dict(alpha=.5), zorder=10)
            sns.stripplot(ax=ax[0], data=low_chrsplit, x="scores", y="group", hue="group", zorder=0, s=2, jitter=0.35)

            sns.boxplot(ax=ax[1], data=low_gensplit, x="scores", y="group", hue="group",
                    meanprops={'marker' : 'D', 'markeredgecolor' : 'black', 'markersize' : 5, 'label':'mean'},
                    #medianprops={'visible': False}, whiskerprops={'visible': False},
                    showmeans=True, showfliers=False, showbox=True, showcaps=False, boxprops=dict(alpha=.5), zorder=10)
            sns.stripplot(ax=ax[1], data=low_gensplit, x="scores", y="group", hue="group", zorder=0, s=2, jitter=0.35)
            ax[0].set_xscale("log", base=10)
            ax[1].set_xscale("log", base=10)
            ax[0].set_xlim(0, max(low_chrsplit["scores"])+25)
            ax[1].set_xlim(0, max(low_gensplit["scores"])+25)
            ax[0].set(ylabel=None, xlabel=None, yticks=[])
            ax[1].set(ylabel=None, xlabel=None, yticks=[])
            
            ax[0].spines[["right", "top"]].set_visible(False)
            ax[0].spines[["left", "bottom"]].set_edgecolor("black")
            ax[1].spines[["right", "top"]].set_visible(False)
            ax[1].spines[["left", "bottom"]].set_edgecolor("black")

            ax[0].legend(labels=[], frameon=False, title=f"$p =$ {conv_mathtext(low_chrval.pvalue)}", facecolor="none", fontsize=5)
            ax[1].legend(labels=[], frameon=False, title=f"$p =$ {conv_mathtext(low_genval.pvalue)}", facecolor="none", fontsize=5)

            ax[0].grid(False)
            ax[0].tick_params(axis="x", which="major", bottom=True, length=3)
            ax[0].tick_params(axis="x", which="minor", bottom=True, length=1.5)
            ax[1].grid(False)
            ax[1].tick_params(axis="x", which="major", bottom=True, length=3)
            ax[1].tick_params(axis="x", which="minor", bottom=True, length=1.5)

            # ax[0,0].set_title("pcf chromosome")
            # ax[0,1].set_title("pcf genome")
            # ax[1,0].set_title("low pass chromosome")
            # ax[1,1].set_title("low pass genome")

            fig.legend(handles=[short_patch, long_patch], frameon=False, facecolor="none")
            plt.grid(False)
            plt.savefig(os.path.join("output", ctx.obj["tag"], f"score_box.pdf"))
            plt.close()






@cli_tools.command(name="gainloss")
@click.pass_context
def gainloss(ctx):
    """
    Plot gains vs losses (4th)

    Some of the cells escaping crisis share copy number variants and are derived form a common subclone in the
    original parental line.
    Therefore need a method to remove the preexisting CN variants in the original subclone.
    The CN variants have already been normalized to the parental line, so now need to try and find out what the
    subclone background is.
    Try and group the samples with similar break points.
    Find a new backgound using wahtever means possible.
    Subtract this background.
    """
    sns.set_style("white")
    cn = pd.read_csv(os.path.join("cnpinter", ctx.obj["tag"], "copyNumber_io", "segmented.copynumber.csv"), index_col=0)
    cn["sampleID"] = cn["sampleID"].astype(str)
    #cn["sampleID"] = ["DB" + str(i) for i in cn["sampleID"]]  # Add back in the DB information
    analysed_samples = os.listdir(os.path.join("cnpinter", ctx.obj["tag"], "normalized_gc_map_wav"))
    analysed_samples = [str(i.split(".")[0].replace("_cov", "")) for i in analysed_samples]
    cn_samp = {samp: v for samp, v in cn.groupby("sampleID")}

    # Background is tumor: normal
    backgrounds, tumour_tl = get_backgrounds(ctx.obj["bgfile"], ctx.obj["bgtc"], ctx.obj["bgts"], ctx.obj["bgnc"])
    backgrounds = {k: backgrounds[k] for k in backgrounds if k in analysed_samples}
    tumour_tl = {k: tumour_tl[k] for k in tumour_tl if k in analysed_samples}
    # print(backgrounds)
    # print(tumour_tl)

    samples_all, tel_lengths = zip(*sorted(tumour_tl.items(), key=lambda x: x[1]))
    print(tumour_tl)

    all_bg = backgrounds

    d_all = pd.concat(cn_samp.values())

    # print(all_bg)
    # print(d_all)

    d_all["Background"] = [all_bg[i] for i in d_all["sampleID"] if i in all_bg]

    data = d_all.copy()

    # quit()
    # all_bg = dict(zip(samples_all, backgrounds))
    # backgrounds = defaultdict(list)
    # for samp, backg in all_bg.items():
    #     if "(p)" in backg:
    #         backgrounds["Parental"].append(samp)
    #         continue
    #     backgrounds[backg].append(samp)

    # Use these samples to get heterozygous SNPs
    #control_samples = set(backgrounds["Parental"])
    #control_samples.remove("DB108")  # FOR NOW ONLY

    #parent_names = {"WT": "DB67", "LIG4": "DB68", "LIG3": "DB69", "LIG3-4": "DB70", "LIG3comp": "DB71", "PARP": "DB72",
    #                "p53-LIG3": "DB73", "ATRX": "DB74"}
    #parent_names_r = {v: k for k, v in parent_names.items()}

    # Get the q and p arm lengths from the hg19 gaps databse
    ref = ctx.obj["ref"]
    if ref.lower() in {"hg38", "grch38"}:
        df_gaps = pd.read_csv("./hg38.gaps.txt", sep="\t", header=0, index_col=0)
    elif ref.lower() in {"hg19", "grch37"}:
        df_gaps = pd.read_csv("./hg19.gaps.txt", sep="\t", header=0, index_col=0)

    df_gaps = df_gaps[df_gaps["type"] == "centromere"]
    df_gaps["chrom"] = [i.strip("chr") for i in df_gaps["chrom"]]

    p_sizes = dict(zip(df_gaps["chrom"], df_gaps["chromStart"]))
    q_sizes = dict(zip(df_gaps["chrom"], df_gaps["chromEnd"]))

    # Use these sizes to calculate the size of the q arm
    normal_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12',
                   'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chrX', 'chrY']
    if ref.lower() in {"hg38", "grch38"}:
        chrom_lengths_dict={'chr1': 248956422, 'chr2': 242193529, 'chr3': 198295559, 'chr4': 190214555, 'chr5': 181538259,
                        'chr6': 170805979, 'chr7': 159345973, 'chr8': 145138636, 'chr9': 138394717, 'chr10': 133797422, 
                        'chr11': 135086622, 'chr12': 133275309, 'chr13': 114364328, 'chr14': 107043718, 'chr15': 101991189, 
                        'chr16': 90338345, 'chr17': 83257441, 'chr18': 80373285, 'chr19': 58617616, 'chr20': 64444167, 
                        'chr21': 46709983, 'chr22': 50818468, 'chrX': 156040895, 'chrY': 57227415} 
    elif ref.lower() in {"hg19", "grch37"}:
        chrom_lengths_dict={'chr1': 249250621, 'chr2': 243199373, 'chr3': 198022430, 'chr4': 191154276, 'chr5': 180915260, 
                        'chr6': 171115067, 'chr7': 159138663, 'chr8': 146364022, 'chr9': 141213431, 'chr10': 135534747, 
                        'chr11': 135006516, 'chr12': 133851895, 'chr13': 115169878, 'chr14': 107349540, 'chr15': 102531392, 
                        'chr16': 90354753, 'chr17': 81195210, 'chr18': 78077248, 'chr19': 59128983, 'chr20': 63025520, 
                        'chr21': 48129895, 'chr22': 51304566, 'chrX': 155270560, 'chrY': 59373566}
    chrom_lengths = [chrom_lengths_dict[i] for i in normal_list]


    chrom_sizes = dict(zip(normal_list, chrom_lengths))

    for k, l in chrom_sizes.items():
        if k in q_sizes:
            q_sizes[k] = l - q_sizes[k]

    genome_size_no_X = sum(chrom_lengths[:-1])


    pd.options.mode.chained_assignment = None  # Ignore annoying pandas warning


    def plot_g_v_l(uni, tag, gl):
        uni = uni[uni["sampleID"].isin(analysed_samples)]
        sl = dict(zip(samples_all, tel_lengths))
        recs = []
        l = []
        for idx, r in uni.groupby("sampleID"):
            sample_tel_length = "{} - {}".format(r.sampleID.iloc[0], format(round(sl[r.sampleID.iloc[0]], 2), ".2f"))
            l.append(float(format(round(sl[r.sampleID.iloc[0]], 2), ".2f")))
            recs.append({"Sample": sample_tel_length,#r.sampleID.iloc[0],
                        "CNV": r.Gain.sum(), "Type": "Gain"})
            recs.append({"Sample": sample_tel_length, #r.sampleID.iloc[0],
                        "CNV": r.Loss.sum(), "Type": "Loss"})

        uni = pd.DataFrame.from_records(recs)
        print(uni)

        # uni_all = uni.copy()

        mpl.rc_context({"lines.linewidth":1})
        sns.set_palette("RdBu", 2, 1)


        fig = plt.figure(figsize=(12,3.2), constrained_layout=True)
        # set axes depending on ignore flags
        if ctx.obj["yaxis_scale"] == False:
            ax = fig.subplot_mosaic([["gl"]])
        else:
            ax = fig.subplot_mosaic([["gl"], ["scale"]], height_ratios=[20, 1])
            ax["scale"].grid(False)

        ax["gl"].grid(False)
        
        #plt.subplots_adjust(left=0.05, bottom=0.4, right=0.9)
        names = ["{} - {}".format(i, format(round(j, 2), ".2f")) for i, j in zip(samples_all, tel_lengths)]
        n_short = len([i for i in tel_lengths if i <= ctx.obj["thresh"]])

        sns.barplot(x="Sample", y="CNV", hue="Type", data=uni, ax=ax["gl"],
                    order=names, # + ["All"],
                    linewidth=0, edgecolor="black", palette={"Gain": "r", "Loss": "b"}, alpha=0.6)
        ax["gl"].set_ylabel("CNV count")
        ax["gl"].set_xlabel("")
        if ctx.obj["yaxis_scale"] == False:
            ax["gl"].set_xticks(ax["gl"].get_xticks(), labels=[i.get_text().replace(ctx.obj["prefix"], "") for i in ax["gl"].get_xticklabels()], rotation=90)
        else:
            ax["gl"].set_xticks([])
        ax["gl"].axvline(x=n_short-0.5, linestyle="--", linewidth=0.2, color="black")

        #plt.yticks(np.arange(0, 21, 2))
        ax["gl"].legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.)
        # plt.legend(loc=2)
        #print([i for i in ax.spines.items()])
        [i[1].set_linewidth(1) for i in ax["gl"].spines.items()]

        if ctx.obj["yaxis_scale"] == True:
            x = [i for i in range(int(len(uni)/2))]
            y = [0]*int(len(uni)/2)
            l = sorted(l)
            heat = cm.jet_r
            norm = mpl.colors.Normalize(vmin=min(l), vmax=max(l))
            m = cm.ScalarMappable(norm=norm, cmap=heat)
            col = [m.to_rgba(i) for i in l]

            ax["scale"].scatter(x, y, c = col, edgecolor="none", marker="s")
            ax["scale"].spines[["top", "bottom", "left", "right"]].set_visible(False)
            ax["scale"].set_xticks([])
            ax["scale"].set_yticks([])
            ax["scale"].text(0, -0.1, f"|\n{round(min(l), 2)}", ha="right", va="center", zorder=100)
            ax["scale"].text(max(x), -0.1, f"|\n{round(max(l), 2)}", ha="right", va="center")
            for i in range(int(min(l))+1, int(max(l))):
                ax["scale"].text(np.argmin(np.abs(np.array(l)-i)), -0.1, f"|\n{i}", ha="right", va="center")

        # Print some stats
        for k, d in uni.groupby("Sample"):
            gains = [count for t, count in zip(d.Type, d.CNV) if t == "Gain"]
            loss = [count for t, count in zip(d.Type, d.CNV) if t == "Loss"]

            print(gains, loss, k, stats.mannwhitneyu(gains, loss))

        spearman_rank = {}
        recursive_partition = {"gain": [], "loss": [], "total": []}
        line_color = {"gain": "red", "loss": "blue", "total": "green"}
        for t in ["gain", "loss", "total"]:
            spearman_rank[t] = stats.spearmanr(gl["tl"], gl[t])
            for i in range(1, len(gl)):
                recursive_partition[t].append(stats.mannwhitneyu(gl[t][:i], gl[t][i:], alternative="greater", axis=0).pvalue)
            ax["gl"].axvline(x=np.argmin(recursive_partition[t])+0.5, color=line_color[t], linestyle="--")

        # statsstr = fr"Spearmans Rank:" + "\n"
        # for t in ["gain", "loss", "total"]:
        #     statsstr += fr"{t}: {round(spearman_rank[t].statistic, 2)} {round(spearman_rank[t].pvalue, 3)}" + "\n"
        # statsstr += fr"Mann Whitney U:" + "\n"
        # text = ax["gl"].text(len(gl)+0.1, max(gl["gain"])/2, statsstr)
        text = ax["gl"].text(len(gl)-0.3, max(gl["gain"])/1.5, "Spearman Rank:\n")
        for t in ["gain", "loss", "total"]:
            text = ax["gl"].annotate(f"{t}: {round(spearman_rank[t].statistic, 2)} {round(spearman_rank[t].pvalue, 3)}\n", xycoords=text, xy=(0,-0.5), verticalalignment="bottom", color=line_color[t])
        text = ax["gl"].annotate("Mann Whitney U:\n", xycoords=text, xy=(0,-0.5), verticalalignment="bottom")
        for t in ["gain", "loss", "total"]:
            text = ax["gl"].annotate(f"{t}: {conv_mathtext(min(recursive_partition[t]))} ({round(gl['tl'].iloc[np.argmin(recursive_partition[t])], 2)})\n", xycoords=text, xy=(0,-0.5), verticalalignment="bottom", color=line_color[t])
            # statsstr += r'\textcolor{red}{gain}: ' + f"{round(min(recursive_partition['gain']), 3)} ({round(gl['tl'].iloc[np.argmin(recursive_partition['gain'])], 3)})" + "\n"
            # statsstr += r'\textcolor{red}{loss}: ' + f"{round(min(recursive_partition['loss']), 3)} ({round(gl['tl'].iloc[np.argmin(recursive_partition['loss'])], 3)})" + "\n"
            # statsstr += r'\textcolor{red}{total}: ' + f"{round(min(recursive_partition['total']), 3)} ({round(gl['tl'].iloc[np.argmin(recursive_partition['total'])], 3)})"
            # statsstr = r'{}'.format(statsstr)



        plt.savefig(os.path.join("output", tag, f"Gains_v_loss.{tag}.pdf"))
        plt.close()

    #data = pd.read_csv("output/CNVs_all.csv") 
    ignore_thresh = 0.05
    print(data.columns)
    data["Gain"] = [1 if i > ignore_thresh else 0 for i in list(data["mean"])]
    data["Loss"] = [1 if i < -ignore_thresh else 0 for i in list(data["mean"])]

    data = data.sort_values("Background")
    #uni = data[(data.start_unique == 1) & (data.end_unique == 1)]

    #common = data[(data.start_unique == 0) & (data.end_unique == 0)]

    recs = []
    for idx, r in data.groupby("sampleID"):
        recs.append({"sample": r.sampleID.iloc[0],
                     "gain": r.Gain.sum(), "loss": r.Loss.sum()})
    gl = pd.DataFrame.from_records(recs)
    gl["total"] = gl["gain"]+gl["loss"]
    ttl = pd.DataFrame.from_records(list(tumour_tl.items()), columns=["sample", "tl"])
    gl = gl.merge(ttl, on="sample", how="left")
    gl = gl.sort_values("tl")
    print("Spearman rank:")
    res = stats.spearmanr(gl["tl"], gl["gain"])
    print("Gains: ", res)
    res = stats.spearmanr(gl["tl"], gl["loss"])
    print("Losses: ", res)
    res = stats.spearmanr(gl["tl"], gl["total"])
    print("Total: ", res)
    gl["group"] = np.where(gl["tl"] <= ctx.obj["thresh"], True, False)
    print("gl df:\n", gl)

    plot_g_v_l(data, ctx.obj["tag"], gl)
    
    def quantify(df, tumour_tl, gl):
        mapd = pd.Series(gl.group.values, index=gl["sample"]).to_dict()
        print(mapd)
        df["group"] = df["sampleID"].map(mapd)
        fig, ax = plt.subplots(1)
        ys, xs, _ = ax.hist(df[df["group"] == True]["mean"], bins=np.arange(min(df["mean"]), max(df["mean"]) + 0.01, 0.01), color="red", alpha=0.5, label=f"below (n={len(gl[gl['group'] == True])})")
        yl, xl, _ = ax.hist(df[df["group"] == False]["mean"], bins=np.arange(min(df["mean"]), max(df["mean"]) + 0.01, 0.01), color="blue", alpha=0.5, label=f"above (n={len(gl[gl['group'] == False])})")

        ax.set_xlim(-2, 2)
        ax.vlines(0.05, 0, max([ys.max(), yl.max()]), color="black")
        ax.vlines(-0.05, 0, max([ys.max(), yl.max()]), color="black")
        # plt.show()
        ax.grid(False)
        ax.legend(title=f"Threshold {ctx.obj['thresh']}:")
        plt.savefig(os.path.join("output", ctx.obj["tag"], f"Gains_v_loss_hist.{ctx.obj['tag']}.pdf"))
        plt.close()

        for i in range(1, len(gl)):
            print(i)
            res = stats.mannwhitneyu(gl[:i]["gain"], gl[i:]["gain"])
            print("Gains: ", res)
            res = stats.mannwhitneyu(gl[:i]["loss"], gl[i:]["loss"])
            print("Loss: ", res)
            res = stats.mannwhitneyu(gl[:i]["total"], gl[i:]["total"])
            print("Total: ", res)

    quantify(data, tumour_tl, gl)
    


@cli_tools.command(name="mergeout")
@click.option("-i", "--indir", default="output")
@click.option("-s", "--skip", is_flag=True, default=False)
@click.pass_context
def mergeout(ctx, indir, skip):
    if skip == False:
        if indir == "output":
            indir = os.path.join("output", ctx.obj["tag"])
        samples = os.listdir(os.path.join(indir, "graphs"))
        samples = [i for i in samples if ".pdf" not in i]
       
        if ctx.obj["plot_inter"] == True:
            total_files = len(samples)*3
            for s in samples:
                total_files += len(os.listdir(os.path.join(indir, "graphs", s)))
                total_files += len(os.listdir(os.path.join(indir, "norm_graphs_no_bg_subtract", s)))
            print(f"Merging {total_files} results pdfs...")

            with Progress() as progress:
                task_merge = progress.add_task("Merging: ", total=total_files)
                #  map_gc_median_diff
                merger = PdfMerger()
                for f in list(sorted(os.listdir(os.path.join(indir, "map_gc_median_diff")))):
                    merger.append(os.path.join(indir, "map_gc_median_diff", f))
                    progress.update(task_merge, advance=1) 
                merger.write(os.path.join(indir, "map_gc_median_diff", f"map_gc_median_diffs.pdf"))
                merger.close

                # merge sample chromosomes and genomes
                for s in samples:
                    merger = PdfMerger()
                    for f in list(sorted(os.listdir(os.path.join(indir, "graphs", s)))):
                        merger.append(os.path.join(indir, "graphs", s, f))
                        progress.update(task_merge, advance=1) 
                    merger.write(os.path.join(indir, "graphs", s, f"{s}_merge.pdf"))
                    merger.close
                
                # merge all chromosomes and genomes
                merger = PdfMerger()
                for s in samples:
                    merger.append(os.path.join(indir, "graphs", s, f"{s}_merge.pdf"))
                    progress.update(task_merge, advance=1) 
                merger.write(os.path.join(indir, "all_merged.pdf"))
                merger.close()


                for s in samples:
                    merger = PdfMerger()
                    for f in list(sorted(os.listdir(os.path.join(indir, "norm_graphs_no_bg_subtract", s)))):
                        merger.append(os.path.join(indir, "norm_graphs_no_bg_subtract", s, f))
                        progress.update(task_merge, advance=1) 
                    merger.write(os.path.join(indir, "norm_graphs_no_bg_subtract", s, f"{s}_no_bg_subtract_merge.pdf"))
                    merger.close

                
                merger = PdfMerger()
                for s in samples:
                    merger.append(os.path.join(indir, "norm_graphs_no_bg_subtract", s, f"{s}_no_bg_subtract_merge.pdf"))
                    progress.update(task_merge, advance=1) 
                merger.write(os.path.join(indir, "all_no_bg_subtract_merged.pdf"))
                merger.close()

        else:
            total_files = len(samples)
            for s in samples:
                total_files += len(os.listdir(os.path.join(indir, "graphs", s)))

            total_files += len(list(glob.glob(f"{indir}/*.pdf")))
            total_files += 1 # all_merged.pdf as not yet made

            with Progress() as progress:
                task_merge = progress.add_task("Merging: ", total=total_files)

                files_existing = 0
                for s in samples:
                    files_existing += len(os.listdir(os.path.join(indir, "graphs", s)))

                if files_existing > 0:
                    # merge sample chromosomes and genomes
                    for s in samples:
                        merger = PdfMerger() 
                        for f in list(sorted(os.listdir(os.path.join(indir, "graphs", s)))):
                            merger.append(os.path.join(indir, "graphs", s, f))
                            progress.update(task_merge, advance=1) 
                        merger.write(os.path.join(indir, "graphs", s, f"{s}_merge.pdf"))
                        merger.close
                    
                    # merge all chromosomes and genomes
                    merger = PdfMerger()
                    for s in samples:
                        merger.append(os.path.join(indir, "graphs", s, f"{s}_merge.pdf"))
                        progress.update(task_merge, advance=1) 
                    merger.write(os.path.join(indir, "all_merged.pdf"))
                    merger.close()

                # merge chromosomes, genomes, and other graphs
                merger = PdfMerger()
                for g in glob.glob(f"{indir}/*.pdf"):
                    merger.append(g)
                    progress.update(task_merge, advance=1) 
                merger.write(os.path.join(indir, "all_cnp_out.pdf"))
                merger.close()
    else:
        if indir == "output":
            indir = os.path.join("output", ctx.obj["tag"])
        total_files = len(list(glob.glob(f"{indir}/*.pdf")))
        with Progress() as progress:
            # merge chromosomes, genomes, and other graphs
            merger = PdfMerger()
            for g in glob.glob(f"{indir}/*.pdf"):
                merger.append(g)
                progress.update(task_merge, advance=1) 
            merger.write(os.path.join(indir, "all_cnp_out.pdf"))
            merger.close()





if __name__ == "__main__":
    cli_tools()
