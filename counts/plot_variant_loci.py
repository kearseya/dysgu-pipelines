#!/usr/bin/python3
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
from sklearn.cluster import KMeans
import click
import collections
from collections import defaultdict
import operator
import re
import os
import pysam

class DictParamType(click.ParamType):
    name = "dictionary"
    def convert(self, cli_value, param, ctx):
        if isinstance(cli_value, dict):
            return cli_value
        try:
            keyvalue_pairs = cli_value.rstrip(",").split(",")
            result_dict = {}
            for pair in keyvalue_pairs:           
                key, value = [item.strip() for item in pair.split(":")]
                if "." in str(value):
                    result_dict[str(key).upper()] = float(value)
                else:
                    result_dict[str(key).upper()] = int(value)
            return result_dict
        except ValueError:
            self.fail(
                "All key-value pairs must be in python dictionary format key1: value1, key2: value2 etc. "
                f"Key-value: {pair}.",
                param,
                ctx
            )

@click.group(name="main", invoke_without_command=True)
@click.argument("indir", type=click.Path(exists=True))
@click.option("-o", "--outdir", default="count_output")
@click.option("-c", "--cyto", default="hg38_cytoBand.txt")
@click.option("-t", "--tl", default="sample_pairs.csv")
@click.option("-s", "--sep-thresh", default=3.81)
@click.option("-e", "--sample-col", default="tumor_db")
@click.option("-a", "--length-col", default="tumor_stela")
@click.option("--size-thresh", default=10000, help="SV size filter")
@click.option("--var-prob", default="INS: 0.45, DEL: 0.45, INV: 0.45, DUP: 0.45, TRA: 0.45", type=DictParamType(), help="SV type prob thresholds")
@click.option("--var-su", default="INS: 9, DEL: 9, INV: 6, DUP: 6, TRA: 8", type=DictParamType(), help="SV type supporting reads thresholds")
@click.option("--ignore-types", default=[], multiple=True, help="var types to ignore in length plot")
@click.pass_context
def main(ctx, indir, outdir, cyto, tl, sep_thresh, sample_col, length_col, size_thresh, var_prob, var_su, ignore_types):
    ctx.ensure_object(dict)
    ctx.obj["indir"] = indir
    ctx.obj["outdir"] = outdir
    ctx.obj["cyto"] = cyto
    ctx.obj["tlfn"] = tl
    ctx.obj["sep_thresh"] = sep_thresh
    ctx.obj["sample_col"] = sample_col
    ctx.obj["length_col"] = length_col
    ctx.obj["size_thresh"] = size_thresh
    ctx.obj["ignore_types"] = ignore_types

    ## allow filtering of vcf without creating new vcf file
    def_prob_thresholds = {"INS": 0.2, "DEL": 0.3, "INV": 0.15, "DUP": 0.15, "TRA": 0.4}
    def_su_thresholds = {"INS": 9, "DEL": 9, "INV": 6, "DUP": 6, "TRA": 8}

    if var_prob == None: # use default
        prob_thresholds = def_prob_thresholds
    else:
        prob_thresholds = var_prob
        for i in ["INS", "DEL", "INV", "DUP", "TRA"]: # if SV type not defined, use default
            if i not in prob_thresholds.keys():
                prob_thresholds[i] = def_prob_thresholds[i]
    if var_su == None:
        su_threhsolds = def_su_thresholds
    else:
        su_thresholds = var_su
        for i in ["INS", "DEL", "INV", "DUP", "TRA"]:
            if i not in su_thresholds.keys():
                su_thresholds[i] = def_su_thresholds[i]

    ctx.obj["var_prob"] = prob_thresholds
    ctx.obj["var_su"] = su_thresholds
    
    if ctx.invoked_subcommand is None:
        if indir == None:
            click.echo(ctx.get_help())
            ctx.exit()
            pass
        else:
            click.echo("Running pipeline")
            ctx.invoke(plot_loci)
            # ctx.invoke(plot_lengths)
    else:
        pass


def extract_chrom(s):
    n = s.replace("chr", "")
    if n == "X":
        n = 23
    elif n == "Y":
        n = 24
    return int(n)


def make_ideogram(inf, order):
    # Mostly stolen from: https://www.biostars.org/p/147364/#147637
    # Path to hg19 cytoband table from ucsc genome browser
    ideo = pd.read_csv(inf, sep="\t")
    ideo.columns = ["chrom", "start", "end", "name", "gieStain"]
    ideo = ideo[~ideo.chrom.str.contains("_")]
    ideo = ideo[~ideo.chrom.str.contains("M")]
    ideo["chrom"] = ideo["chrom"].str.replace("chr", "")
    # ideo[ideo["chrom"] == "X"]["chrom"] = 23
    # ideo[ideo["chrom"] == "Y"]["chrom"] = 24
    ideo["chrom"] = ideo["chrom"].str.replace("X", "23")
    ideo["chrom"] = ideo["chrom"].str.replace("Y", "24")

    ideo["chrom"] = ideo["chrom"].astype(int)

    # Colors for different chromosome stains
    color_lookup = {'gneg': (1., 1., 1.),
                    'gpos25': (.6, .6, .6),
                    'gpos50': (.4, .4, .4),
                    'gpos75': (.2, .2, .2),
                    'gpos100': (0., 0., 0.),
                    'acen': ("red"),
                    'gvar': (.8, .8, .8),
                    'stalk': (.9, .9, .9)}

    # Add a new column for colors
    ideo['colors'] = ideo['gieStain'].apply(lambda x: color_lookup[x])

    # Add a width column
    ideo['width'] = ideo['end'] - ideo['start']
    ideo.to_csv("ideogram.csv")

    chrom = list(ideo["chrom"])
    start = list(ideo["start"])
    width = list(ideo["width"])
    clr = list(ideo["colors"])

    collections = defaultdict(list)
    colors = defaultdict(list)

    for i, chr in enumerate(chrom):
        if chr not in order:
            continue
        collections[chr].append([start[i], width[i]])
        colors[chr].append(clr[i])

    return collections, colors


@main.command("loci")
@click.option("-s", "--save", default=False, is_flag=True, help="save results to file")
@click.pass_context
def loci(ctx, save):
    indir = ctx.obj["indir"]
    prob_thresholds = ctx.obj["var_prob"]
    su_thresholds = ctx.obj["var_su"]
    if os.path.isfile(indir):
        if indir.lower().endswith((".vcf", ".csv")):
            files = [indir]
        else:
            print("Cannot process file format")
            exit
    elif os.path.isdir(indir):
        files = os.listdir(indir)
        files = [i for i in files if i.lower().endswith((".vcf", ".csv"))]
        files = [os.path.join(indir, i) for i in files]
    failed = []
    total = 0
    loci_arr = []
    types = ["DEL", "INS", "INV", "DUP", "TRA"]
    chroms = [i for i in range(1, 24)]
    for v in files:
        total += 1
        fmt = v[-3:]
                
        if fmt == "vcf":
            f = pysam.VariantFile(v)
            s = os.path.basename(v).replace(".vcf", "")

            for line in f.fetch():
                #if list(line.filter.keys())[0] != "lowProb" or float(line.samples[s]["PROB"]) > 0.2 or float(line.samples[s]["SU"]) > 4: 
                if "_" in line.chrom:
                    continue
                sv_type = line.info["SVTYPE"]
                # chroms.add(extract_chrom(line.chrom))
                if float(line.samples[s]["PROB"]) >= prob_thresholds[sv_type]:
                    if (sv_type != "TRA" and line.stop-line.pos >= ctx.obj["size_thresh"]) or sv_type in {"TRA", "INS"}:
                        loci_arr.append((line.pos, line.stop, sv_type, extract_chrom(line.chrom)))
        
        # elif fmt == "csv":
        #     counts[s] = {}
        #     df = pd.read_csv(i)
        #     # counts[s]["lowProb"] = dict(df[df["filter"] == "lowProb"]["svtype"].value_counts())
        #     # counts[s]["PASS"] = dict(df[df["filter"] == "PASS"]["svtype"].value_counts())
        #     counts[s] = dict(df["svtype"].value_counts())
        #     chroms = [i for i in range(1, 25)]

        else:
            print(f"cannot process {i}, invalid format")

    df = pd.DataFrame.from_records(loci_arr)
    df.columns = ["start", "end", "type", "chrom"]
    print(df)

    return df, chroms

    # df = pd.DataFrame.from_dict(counts).transpose()
    # df["total"] = df.sum(axis=1, numeric_only=True)
    # df["sample"] = df.index
    # df = df.reset_index(drop=True)
    # df = df[["sample"] + types + ["total"]]
    # if save:
    #     df.to_csv(out, index=False)
    # return df, list(chroms)


def process_overlap(df, chroms):
    bars = {}
    points = {}
    bart = ["DUP", "INV", "DEL"]
    df["size"] = df["end"]-df["start"]
    ymaxs = {}
    for c in chroms:
        t = df[df["chrom"] == c]
        t = t.sort_values("start")
        t = t.reset_index(drop=True)
        b = t[t["type"].isin(bart)]
        p = t[~t["type"].isin(bart)]
        # print(b)
        #a = list(b.itertuples(index=False, name=None))
        a = list(b.to_numpy())
        a = [list(i) for i in a]
        if len(a) == 0:
            print("NONE ON CHR: ", c)
            bars[c] = []
            ymaxs[c] = 1
        else:
            a[0].append(1)
            y = 1
            ymax = 1
            for i, x in enumerate(a[:-1]):
                if a[i+1][0] < x[1]:
                    y += 1
                    if y > ymax:
                        ymax = y
                    if c == 17:
                        print(y, x, a[i+1])
                else:
                    y = 1
                a[i+1].append(y)

        ups = list(p.to_numpy())
        ups = [list(i) for i in ups]

        for i, _ in enumerate(ups):
            ups[i].append(ymax+1)

        # print(c, len(ups))
        
        ps = {"INS": [i for i in ups if i[2] == "INS"], "TRA": [i for i in ups if i[2] == "TRA"]}

        # print(c, len(ps))

        bars[c] = list(a)
        ymaxs[c] = ymax
        points[c] = ps

    print(len(points))

    return bars, points, ymaxs



@main.command("plot_loci")
@click.option("-s", "--save", default=False, is_flag=True, help="plot show or plot save")
@click.option("--scale", default="linear",
    type=click.Choice(["linear", "log", "symlog", "logit"]))
@click.option("-k", "--kgroups", default=5)
@click.option("-r", "--real-groups", default=None)
@click.option("--ksort", default=False, is_flag=True)
@click.option("--lsort", default=True, is_flag=True, help="sort by length col")
@click.option("--bar", default=True, is_flag=True)
@click.option("--bar-lines", default=[None], multiple=True)#["length", "kmeans"], multiple=True)
@click.option("--bar-numbers", default=[None], multiple=True)
@click.option("--bar-scale", default=True, is_flag=True)
@click.option("--box", default=True, is_flag=True)
@click.option("--hist", default=False, is_flag=True)
@click.pass_context
def plot_loci(ctx, save, scale, kgroups, real_groups, ksort, lsort, bar, bar_lines, bar_numbers, bar_scale, box, hist):
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
    df, chroms = ctx.invoke(loci)
    var = ["DEL", "INS", "INV", "DUP", "TRA"] # df.columns[2:].tolist()
    ubars, points, ymaxs = process_overlap(df, chroms)
    
    type_cols = {"INS": "green", "DEL": "red", "INV": "blue", "DUP": "magenta", "TRA": "black"}
    point_shape = {"INS": ".", "TRA": "|"}
    
    ideograms, ideo_colors = make_ideogram(ctx.obj["cyto"], chroms)
    # create subplots
    fig, ax = plt.subplots(ncols=1, nrows=len(chroms))
    
    for c in chroms:
        i = c-1
        idg = ideograms[c]
        idc = ideo_colors[c]
        
        ymax = ymaxs[c]
        ubar = ubars[c]
        point = points[c]

        ax[i].broken_barh(idg, (-ymax, ymax), facecolors=idc)
        
        for y in range(1, ymax+1):
            yarr = [i for i in ubar if i[-1] == y]
            bars = [(i[0], i[4]) for i in yarr]
            colors = [type_cols[i[2]] for i in yarr]
            if len(bars)> 0:
                ax[i].broken_barh(bars, (y, y+1), facecolors=colors, alpha=0.25)
        for t in point:
            if len(point[t]) > 0:
                ax[i].scatter([i[0] for i in point[t]], [i[-1] for i in point[t]], color=type_cols[t], marker=point_shape[t])

        ax[i].set_ylabel(c)
        ax[i].set_xlim(0, idg[-1][0]+idg[-1][1])
    plt.show()


if __name__ == "__main__":
    main()
