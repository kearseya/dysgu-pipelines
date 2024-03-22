#!/usr/bin/python3
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import scipy
from sklearn.cluster import KMeans
import click
import collections
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
@click.option("-t", "--tl", default="sample_pairs.csv")
@click.option("-s", "--sep-thresh", default=3.81)
@click.option("-e", "--sample-col", default="tumor_db")
@click.option("-a", "--length-col", default="tumor_stela")
@click.option("--size-thresh", default=10000, help="SV size filter")
@click.option("--var-prob", default="INS: 0.2, DEL: 0.3, INV: 0.15, DUP: 0.15, TRA: 0.4", type=DictParamType(), help="SV type prob thresholds")
@click.option("--var-su", default="INS: 9, DEL: 9, INV: 6, DUP: 6, TRA: 8", type=DictParamType(), help="SV type supporting reads thresholds")
@click.option("--ignore-types", default=[], multiple=True, help="var types to ignore in length plot")
@click.pass_context
def main(ctx, indir, outdir, tl, sep_thresh, sample_col, length_col, size_thresh, var_prob, var_su, ignore_types):
    ctx.ensure_object(dict)
    ctx.obj["indir"] = indir
    ctx.obj["outdir"] = outdir
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
            ctx.invoke(plot_counts)
            # ctx.invoke(plot_lengths)
    else:
        pass


@main.command("count")
@click.option("-s", "--save", default=False, is_flag=True, help="save results to file")
@click.pass_context
def count(ctx, save):
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
    counts = {}
    types = ["DEL", "INS", "INV", "DUP", "TRA"]
    for v in files:
        total += 1
        fmt = v[-3:]
                
        if fmt == "vcf":
            f = pysam.VariantFile(v)
            s = os.path.basename(v).replace(".vcf", "")
            counts[s] = {}
            for t in types:
                counts[s][t] = 0

            for line in f.fetch():
                #if list(line.filter.keys())[0] != "lowProb" or float(line.samples[s]["PROB"]) > 0.2 or float(line.samples[s]["SU"]) > 4: 
                sv_type = line.info["SVTYPE"]
                if float(line.samples[s]["PROB"]) >= prob_thresholds[sv_type]:
                    if (sv_type != "TRA" and line.stop-line.pos >= ctx.obj["size_thresh"]) or sv_type in {"TRA", "INS"}:
                        counts[s][sv_type] += 1
        
        elif fmt == "csv":
            counts[s] = {}
            df = pd.read_csv(i)
            # counts[s]["lowProb"] = dict(df[df["filter"] == "lowProb"]["svtype"].value_counts())
            # counts[s]["PASS"] = dict(df[df["filter"] == "PASS"]["svtype"].value_counts())
            counts[s] = dict(df["svtype"].value_counts())

        else:
            print(f"cannot process {i}, invalid format")

    df = pd.DataFrame.from_dict(counts).transpose()
    df["total"] = df.sum(axis=1, numeric_only=True)
    df["sample"] = df.index
    df = df.reset_index(drop=True)
    df = df[["sample"] + types + ["total"]]
    if save:
        df.to_csv(out, index=False)
    return df


@main.command("length")
@click.option("-s", "--save", default=False, is_flag=True, help="save results to file")
@click.pass_context
def length(ctx, save):
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
    lengths = {}
    types = ["DEL", "INS", "INV", "DUP"]
    types = [i for i in types if i not in ctx.obj["ignore_types"]]
    for v in files:
        total += 1
        fmt = v[-3:]
                
        if fmt == "vcf":
            f = pysam.VariantFile(v)
            s = os.path.basename(v).replace(".vcf", "")
            lengths[s] = {}
            for t in types:
                lengths[s][t] = []

            for line in f.fetch():
                #if list(line.filter.keys())[0] != "lowProb" or float(line.samples[s]["PROB"]) > 0.2 or float(line.samples[s]["SU"]) > 4: 
                sv_type = line.info["SVTYPE"]
                if float(line.samples[s]["PROB"]) >= prob_thresholds[sv_type]:
                    if (sv_type != "TRA" and max([line.stop-line.pos, line.info["SVLEN"]]) >= ctx.obj["size_thresh"]):
                    #if sv_type not in {"TRA"}:
                        lengths[s][sv_type].append(line.info["SVLEN"])
        
        elif fmt == "csv":
            lengths[s] = {}
            df = pd.read_csv(i)
            # lengths[s]["lowProb"] = dict(df[df["filter"] == "lowProb"]["svtype"].value_lengths())
            # lengths[s]["PASS"] = dict(df[df["filter"] == "PASS"]["svtype"].value_lengths())
            lengths[s] = dict(df["svtype"].value_lengths())

        else:
            print(f"cannot process {i}, invalid format")

    print(lengths)

    return lengths

@main.command("plot_counts")
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
def plot_counts(ctx, save, scale, kgroups, real_groups, ksort, lsort, bar, bar_lines, bar_numbers, bar_scale, box, hist):
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
    df = ctx.invoke(count)
    var = ["DEL", "INS", "INV", "DUP", "TRA"] # df.columns[2:].tolist()
    ## reorder table
    if ctx.obj["tlfn"] == None:
        df = df.sort_values(by="total")
    else:
        ldf = pd.read_csv(ctx.obj["tlfn"])
        if ctx.obj["sample_col"] != "sample":
            ldf = ldf.rename({ctx.obj["sample_col"]: "sample"}, axis=1)
        if ctx.obj["length_col"] != "length":
            ldf = ldf.rename({ctx.obj["length_col"]: "length"}, axis=1)
        ldf = ldf[["sample", "length"]]
        df = pd.merge(df, ldf, on="sample")
        print(df)
        print(df["length"].dtype)
        
        if df["length"].dtype != bool:
            df["length"] = df["length"].astype(float)
            df["short"] = np.where(df["length"] <= float(ctx.obj["sep_thresh"]), True, False)
        #print(df)
        if df["length"].dtype == bool:
            df["short"] = df["length"]
            df = df.sort_values(by=["length", "total"])
        else:
            if lsort == False:
                df = df.sort_values(by=["short", "total"])
            else:
                df = df.sort_values(by="length", ascending=False)

        if "length" in bar_lines:
            if df["length"].dtype != bool:
                for j, z in enumerate(df["length"]):
                    if z <= ctx.obj["sep_thresh"]:
                        lsplit = j
                        print(lsplit)
                        break
            else:
                for j, z in enumerate(df["length"]):
                    if z == True:
                        lsplit = j
                        print("bool", lsplit)
                        break


    df = df.reset_index()

    ## add kmeans labels ##
    kmeans = KMeans(n_clusters=kgroups, random_state=0).fit(df[var])
    df["kmeans"] = kmeans.labels_
    #print(kmeans.labels_)
    #print(df["kmeans"].tolist())
    #print(df["total"].tolist())
    #print(df["kmeans"])
    #with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
    #   print(df)
    if ksort == True:
        df = df.sort_values(by=["total", "kmeans"])

    if "kmeans" in  bar_lines:
        for i, x in enumerate(df["kmeans"].tolist()):
            if x > df["kmeans"].iloc[0]:
                ksplit = i
                break
    #if "length" in bar_lines:
    #   if length != None:
    #       for j, z in enumerate(df["length"]):
    #           if z <= threshold:
    #               lsplit = j
    #               break
    #   else:
    #       lsplit = -10


    ## stacked bar plot ##
    if bar == True:
        xticks = np.arange(0, len(df))
        labels = df["sample"].tolist()#[i.replace("_filtered.vcf", "") for i in df["sample"].tolist()]
        fig, ax = plt.subplots()
        for i, v in enumerate(var):
            if i == 0:
                ax.bar(xticks, df[v].tolist(), 0.5, label = v, log=True)
            else:
                ax.bar(xticks, df[v].tolist(), 0.5, label = v, bottom=df[var[:i]].sum(axis=1), log=True)
        #df.plot(kind="bar", stacked=True)
        plt.ylim([0, max(df["total"])])
        plt.tight_layout()
        plt.xlim([xticks[0]-1, xticks[-1]+1])
        plt.yscale(scale)
        plt.xticks(ticks = xticks, labels = labels, rotation = 90)
        plt.legend(var, loc = "upper left")
        if "kmeans" in  bar_lines:
            plt.axvline(x = ksplit-0.5, linestyle = "--", color="black")
        if "length" in bar_lines:
            plt.axvline(x = lsplit-0.5, linestyle = "--", color="red") 
        #with pd.option_context('display.max_rows', None, 'display.max_columns', None):  # more options can be specified also
        #   print(df)
        
        if "kmeans" in  bar_numbers:
            for la, xp, yp in zip(df["kmeans"], xticks, df["total"]):
                plt.annotate(la, (xp-0.125, yp+1))
        
        if ctx.obj["tlfn"] != None:
            if bar_scale == True:
                print(max(df["total"]))
                sbh = max(df["total"])/20
                scale_bar = ax.barh(-5/2, max(xticks), sbh)
                gradientbars(scale_bar, df["length"])
                norm=mpl.colors.Normalize(vmin=min(df["length"]),vmax=max(df["length"]))
                plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap="RdYlGn"), ax=ax)
                plt.ylim(-5/2, max(df["total"])+3)
        if save == False:
            plt.show()
        else:
            plt.savefig(f"{ctx.obj['outdir']}/variant_type_stacked_bat_plot.png")
            plt.clf()
    
    ## box plot ##
    if box == True:
        var_types = {"long": [], "short": []}
        mann_whitney = []
        is_bool = False
        if df["length"].dtype == bool:
            is_bool = True
        var.append("total")
        for v in var:
            if is_bool == True:
                var_types["long"].append(df[df["length"] == False][v].tolist())
                var_types["short"].append(df[df["length"] == True][v].tolist())
                #var_types[v] = [df[df["length"] == False][v].tolist(), df[df["length"] == True][v].tolist()]
            else:
                var_types["long"].append(df[df["short"] == False][v].tolist())
                var_types["short"].append(df[df["short"] == True][v].tolist())

        for i in range(len(var_types["long"])):
            mann_whitney.append(scipy.stats.mannwhitneyu(var_types["long"][i], var_types["short"][i]))
        #print(mann_whitney)
                #var_types[v] = [df[df["length"] < threshold][v].tolist(), df[df["length"] >= threshold][v].tolist()]
        #print(var_types)
        lbox = plt.boxplot(var_types["long"], positions=np.array(range(len(var_types["long"])))*2.0-0.4, sym='', widths=0.6)
        plt.setp(lbox["boxes"], color = "blue")
        plt.setp(lbox["whiskers"], color = "blue")
        plt.setp(lbox["caps"], color = "blue")
        plt.setp(lbox["medians"], color = "blue")
        sbox = plt.boxplot(var_types["short"], positions=np.array(range(len(var_types["short"])))*2.0+0.4, sym='', widths=0.6)
        plt.setp(sbox["boxes"], color = "red")
        plt.setp(sbox["whiskers"], color = "red")
        plt.setp(sbox["caps"], color = "red")
        plt.setp(sbox["medians"], color = "red")
        plt.plot([], c="blue", label="long")
        plt.plot([], c="red", label="short")
        plt.legend()
        plt.xticks(ticks=range(0, len(var_types["short"])*2, 2), labels=var)
        plt.ylim([0, max(df["total"])])
        for x, i in enumerate(mann_whitney):
            pval = i[1]
            if pval > 0.05:
                continue
            elif 0.01 < pval <=0.5:
                plt.annotate("*", (x*2, max(var_types["long"][x]+var_types["short"][x])), ha='center', fontsize="large")
            elif 0.001 < pval <= 0.01:
                plt.annotate("**", (x*2, max(var_types["long"][x]+var_types["short"][x])), ha='center', fontsize="large")
            elif pval <= 0.001:
                plt.annotate("***", (x*2, max(var_types["long"][x]+var_types["short"][x])), ha='center', fontsize="large")
        if save == False:
            plt.show()
        else:
            plt.savefig(f"{ctx.obj['outdir']}/variant_type_boxplot.png")
            plt.clf()

    ## histogram
    if hist == True:
        for v in var:
            plt.hist(df[v], alpha=0.5)
        plt.legend(var)
        if save == False:
            plt.show()
        else:
            plt.savefig(f"{ctx.obj['outdir']}/variant_type_histogram.png")
            plt.clf()


def gradientbars(bars, data):
    ax = bars[0].axes
    lim = ax.get_xlim()+ax.get_ylim()
    for bar in bars:      
        bar.set_zorder(1)
        bar.set_facecolor("none")
        x,y = bar.get_xy()
        w, h = bar.get_width(), bar.get_height()
        #grad = np.atleast_2d(np.linspace(min(data),1*w/max(data),256))
        ax.imshow(np.atleast_2d(data), extent=[x,x+w,y,y+h], aspect="auto", zorder=0, norm=mpl.colors.Normalize(vmin=min(data),vmax=max(data)), cmap="RdYlGn")
    ax.axis(lim)   
    #print(df)


@main.command("plot_length")
@click.option("-s", "--save", default=False, is_flag=True, help="show or save plot")
@click.pass_context
def plot_lengths(ctx, save):
    mpl.rcParams['axes.spines.right'] = False
    mpl.rcParams['axes.spines.top'] = False
    l = ctx.invoke(length)
    var = ["DEL", "INS", "INV", "DUP"]
    var = [i for i in var if i not in ctx.obj["ignore_types"]]
    if ctx.obj["tlfn"] == None:
        print("Requires length file to split groups")
    else:
        ldf = pd.read_csv(ctx.obj["tlfn"])
        if ctx.obj["sample_col"] != "sample":
            ldf = ldf.rename({ctx.obj["sample_col"]: "sample"}, axis=1)
        if ctx.obj["length_col"] != "length":
            ldf = ldf.rename({ctx.obj["length_col"]: "length"}, axis=1)
        ldf = ldf[["sample", "length"]]
    
    if ldf["length"].dtype != bool:
        ldf["length"] = ldf["length"].astype(float)
        ldf["short"] = np.where(ldf["length"] <= float(ctx.obj["sep_thresh"]), True, False)
    #print(ldf)
    if ldf["length"].dtype == bool:
        ldf["short"] = ldf["length"]
    ldf = ldf[ldf["sample"].isin(l.keys())]

    short_samples = ldf[ldf["short"] == True]["sample"]
    long_samples = ldf[ldf["short"] == False]["sample"]

    type_lengths = {"long": {}, "short": {}}
    for v in var:
        type_lengths["long"][v] = []
        type_lengths["short"][v] = []

    for s in long_samples:
        for v in var:
            type_lengths["long"][v] += l[s][v]
    for s in short_samples:
        for v in var:
            type_lengths["short"][v] += l[s][v]


    #var_lengths = {"long": [], "short": []}
    
    mann_whitney = []
    
    # for c in ["long", "short"]:
    #     for v in var:
    #         var_lengths[c].append(type_lengths[c][v])

    for v in var:
        try:
            mann_whitney.append(scipy.stats.mannwhitneyu(var_lengths["long"][v], var_lengths["short"][v]))
        except:
            mann_whitney.append([1, 1])
    #print(mann_whitney)
            #var_lengths[v] = [df[df["length"] < threshold][v].tolist(), df[df["length"] >= threshold][v].tolist()]
    #print(var_lengths)
    fig, ax = plt.subplots(1, len(var))
    for i, v in enumerate(var):
        lbox = ax[i].boxplot(type_lengths["long"][v], positions=[0], sym='', widths=0.6, boxprops={"color": "blue"}, whiskerprops={"color": "blue"}, flierprops={"color": "blue"}, medianprops={"color": "blue"})
        # ax[i].set(lbox["boxes"], color = "blue")
        # ax[i].set(lbox["whiskers"], color = "blue")
        # ax[i].set(lbox["caps"], color = "blue")
        # ax[i].set(lbox["medians"], color = "red")
        sbox = ax[i].boxplot(type_lengths["short"][v], positions=[1], sym='', widths=0.6, boxprops={"color": "red"}, whiskerprops={"color": "red"}, flierprops={"color": "red"}, medianprops={"color": "red"})
        # ax[i].set(sbox["boxes"], color = "red")
        # ax[i].set(sbox["whiskers"], color = "red")
        # ax[i].set(sbox["caps"], color = "red")
        # ax[i].set(sbox["medians"], color = "red")
        ax[i].spines[["right", "top"]].set_visible(False)
        if i == 0:
            ax[i].plot([], c="blue", label="long")
            ax[i].plot([], c="red", label="short")
        #ax[i].legend()
        ax[i].set_xticks(ticks=[0.5], labels=[v])
        pval = mann_whitney[i][1]
        if pval > 0.05:
            continue
        elif 0.01 < pval <=0.5:
            ax[i].annotate("*", (0.5, max(type_lengths["long"][v]+type_lengths["short"][v])), ha='center', fontsize="large")
        elif 0.001 < pval <= 0.01:
            ax[i].annotate("**", (0.5, max(type_lengths["long"][v]+type_lengths["short"][v])), ha='center', fontsize="large")
        elif pval <= 0.001:
            ax[i].annotate("***", (0.5, max(type_lengths["long"][v]+type_lengths["short"][v])), ha='center', fontsize="large")
    fig.legend(loc="upper right")

    if save == False:
        plt.show()
    else:
        plt.savefig(f"{ctx.obj['outdir']}/variant_length_boxplot.png")
        plt.clf()






if __name__ == "__main__":
    main()
