import pycircos
import pysam
import matplotlib as mlp
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#from collections import defaultdict
import collections
import operator
import re
import os
import patchworklib as pw
from math import ceil
import cairosvg
import PIL
import svgutils
import svgutils.transform as sg
from svgutils.compose import *
from rich.progress import Progress
from rich import print
import shutil
import random

import click

# def color_selector(tl):
#     if float(tl) <= 3.81:
#         return "red"
#     else:
#         return "green"

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

def read_file(fp, outdir, tl, cytoband, cnp, prob_thresholds, su_thresholds, size_threshold):
    """
    Read (dysgu) vcf file and convert to pycircos Gcircle
    """
    color_dict   = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777", "gpos100":"#000000",
                    "gvar":"#FFFFFF00", "stalk":"#C01E27", "acen":"#D82322"}
    values_all   = []
    arcdata_dict = collections.defaultdict(dict)

    low_size = 800
    up_size = 850

    aof = 45
    dof = 44
    #print(fp)

    # data from INFO string
    #reend      = re.search("END=(\d+)")
    #relen      = re.search("SV_LEN=(\d+)")
    #rename2    = re.search("CHR2=(\d+)")
    #restart2   = re.search("CHR2_POS=(\d+)")
    #reend2     = re.search("CHR2_POS=(\d+)")

    ##Set chromosomes
    circle = pycircos.Gcircle(figsize=(5, 5))
    #circle = pycircos.Gcircle(fig=pw.Brick._figure, figsize=(5, 5))
    ## add tl to label in middle of circle
    #arc = pycircos.Garc(label=str(tl), label_visible=True, size=0, facecolor=color_selector(tl), edgecolor=color_selector(tl), labelcolor=color_selector(tl), labelsize=20, raxis_range=(500,1000))
    #circle.add_garc(arc)
    #circle.fillplot(garc_id=None, data=[0, low_size],raxis_range=(0, low_size), facecolor=colour_selector(tl))
    nl = {}
    with open(fp) as f:
        f.readline()
        for line in f:
            if line.startswith("##contig"):
                #line   = line.split("<")
                #print(line)
                name    = re.search("ID=([a-zA-Z0-9_]*)", line).group(1).replace("chr", "")
                length  = int(re.search("length=(\d+)", line).group(1))
                #print(name, length)
                if "_" not in name:
                    #print(name, length)
                    nl[name] = length
    #print(nl)
    #sorted_nl = dict( sorted(nl.items(), key=operator.itemgetter(1),reverse=True))
    order_of_chr = [f"{i}" for i in range (1, 23)] + ["X", "Y"]
    for i in order_of_chr:
        name    = i
        length  = nl[i]
        arc     = pycircos.Garc(arc_id=name, size=length, interspace=2, raxis_range=(low_size, up_size), 
                       labelposition=100, label_visible=True, facecolor="white", labelsize=17)
        circle.add_garc(arc)
    circle.set_garcs()

    arcdata_dict = collections.defaultdict(dict)
    ## add cytoband bars to arcs
    if cytoband != None:
        with open(cytoband) as f:
            f.readline()
            for line in f:
                line  = line.rstrip().split("\t")
                name  = line[0].replace("chr", "")
                start = int(line[1])-1
                width = int(line[2])-(int(line[1])-1)
                if name in order_of_chr:
                    if name not in arcdata_dict:
                        arcdata_dict[name]["positions"] = []
                        arcdata_dict[name]["widths"]    = []
                        arcdata_dict[name]["colors"]    = []
                    arcdata_dict[name]["positions"].append(start)
                    arcdata_dict[name]["widths"].append(width)
                    arcdata_dict[name]["colors"].append(color_dict[line[-1]])

        for key in arcdata_dict:
            circle.barplot(key, data=[1]*len(arcdata_dict[key]["positions"]), positions=arcdata_dict[key]["positions"],
                           width=arcdata_dict[key]["widths"], raxis_range=[low_size, up_size], facecolor=arcdata_dict[key]["colors"])

    ## add copynumber data
    if cnp != None:
        low_size = 750
        up_size = 800
        cndf = pd.read_csv(cnp)
        sample = os.path.basename(fp).replace(".vcf", "")
        cn = cndf[cndf["sampleID"] == sample]
        cn["len"] = cn["end.pos"] - cn["start.pos"]
        cn["clip"] = cn["mean"].clip(upper=3)
        cn["chrom"] = cn["chrom"].astype(str)
        for index, row in cn.iterrows():
            circle.lineplot(row["chrom"], positions=[row["start.pos"], row["end.pos"]], data=[row["clip"], row["clip"]], 
                    rlim=(-2, 3), raxis_range=(700, 800), linecolor="black")
        for c in set(cn["chrom"].tolist()):
            plen = cn[cn["chrom"] == c]["len"].tolist()
            totm = sum([abs(i) for i in cn[cn["chrom"] == c]["mean"].tolist()])
            # if len(plen) > 8: # highlight complex copy number chromosome segments
            if totm >= 10 or len(plen) > 8:
                plen = [i/max(arcdata_dict[c]["positions"]) for i in plen]
                starts = cn[cn["chrom"] == c]["start.pos"].tolist()
                ends = cn[cn["chrom"] == c]["end.pos"].tolist()
                frag = np.where(np.array(plen) <= 10)[0].tolist()
                nfrag = len(frag)
                fseg = [starts[min(frag)], ends[max(frag)]]
                if 10 <= totm < 30:
                    circle.barplot(c, data=[1]*2, positions=[fseg[0]], width=[fseg[1]-fseg[0]],
                            rlim=(0,1), raxis_range=(700, 800), facecolor="yellow")
                elif 30 <= totm < 50:
                    circle.barplot(c, data=[1]*2, positions=[fseg[0]], width=[fseg[1]-fseg[0]],
                            rlim=(0,1), raxis_range=(700, 800), facecolor="orange")
                elif 50 <= totm:
                    circle.barplot(c, data=[1]*2, positions=[fseg[0]], width=[fseg[1]-fseg[0]],
                            rlim=(0,1), raxis_range=(700, 800), facecolor="red")
                else:
                #elif len(plen) > 8:
                    circle.barplot(c, data=[1]*2, positions=[fseg[0]], width=[fseg[1]-fseg[0]],
                            rlim=(0,1), raxis_range=(700, 800), facecolor="#ffcf91")




    ## colour scheme
    #type_cols = {"INS": "magenta", "DEL": "red", "INV": "blue", "DUP": "green", "TRA": "black"}
    # pastel
    #type_cols = {"INS": "#FFCE54", "DEL": "#ED5564", "INV": "#4FC1E8", "DUP": "#ADD568", "TRA": "#AC92EB"}
    type_cols = {"INS": "yellow", "DEL": "red", "INV": "deepskyblue", "DUP": "lime", "TRA": "darkorchid"}
    
    f = pysam.VariantFile(fp)
    tra_list = []
    ins_dict = collections.defaultdict(dict)
    for i in arcdata_dict.keys():
        ins_dict[i]["positions"] = []
        ins_dict[i]["widths"] = []
    for l in f.fetch():
        # if l.filter.keys()[0] == "lowProb": ## swapped for individual prob thresh for each type
        #     continue
        name = l.chrom.replace("chr", "")
        sv_type = l.info["SVTYPE"]
        sample_name = list(l.samples.keys())[0]
        if l.samples[sample_name]["PROB"] < prob_thresholds[sv_type]: # or l.samples[sample_name]["SU"] < type_su_thresholds[sv_type]:
            if l.samples[sample_name]["SU"] < su_thresholds[sv_type]:
                continue
        start = l.pos
        end = l.stop #info["END"]
        if sv_type in {"DEL", "INV", "DUP"}: #, "INS"}:
            if abs(end-start) <= size_threshold: # remove small deletions and inversions (default 500bp)
                continue
        chr1 = l.contig.replace("chr", "")
        chr2 = l.info["CHR2"].replace("chr", "")
        if chr1 == chr2 and name in order_of_chr:
            source = (name, start, start, low_size-aof)
            destination = (name, end, end, low_size-aof)
            #print(source, destination)
            if abs(end-start) < 50000: # required for line to show
                destination = (name, end+50000, end+50000, low_size-aof) 
            if sv_type != "INS":
                circle.chord_plot(source, destination, facecolor=type_cols[sv_type], edgecolor=type_cols[sv_type], linewidth=1.0)
            else:
                ins_dict[str(name)]["positions"].append(start)
                #print(f"end-start: {end-start}, rlen: {l.rlen}")
                #ins_dict[str(name)]["widths"].append(end-start)
        elif chr1 != chr2 and name in order_of_chr and chr2 in order_of_chr: 
            name2 = l.info["CHR2"].replace("chr", "")
            start2 = l.info["CHR2_POS"]
            end2 = start2
            tra_list.append(((name, start, start, low_size-aof), (name2, start2, end2, low_size-aof)))
            #source = (name, start, start, low_size)
            #destination = (name2, start2, end2, low_size)
            #circle.chord_plot(source, destination, facecolor="#000000", linewidth=0.5)
        else:
            continue
    for t in tra_list:
        circle.chord_plot(t[0], t[1], facecolor=type_cols["TRA"], edgecolor=type_cols["TRA"], linewidth=1.0)
    for key in ins_dict:
        if len(ins_dict[key]["positions"]) == 0:
            continue
        circle.scatterplot(key, data=[1]*len(ins_dict[key]["positions"]), positions=ins_dict[key]["positions"], 
                raxis_range=(low_size-aof, low_size-dof), facecolor=type_cols["INS"], edgecolor="green",  markershape=".", markersize=20)
        # circle.barplot(key, positions=ins_dict[key]["positions"], width=ins_dict[key]["widths"], 
        #         data=[1]*len(arcdata_dict[key]["positions"]),
        #         raxis_range=(812, 838), facecolor="green")#, edgecolor="green")#,  markershape=".", markersize=28)

    return circle


def arrange(indir, outdir, outname, ncols, lengths, col, col_thresh, progress, subsets, subsamples, split):
    """
    Combine all circos plots together
    """
    svgs = os.listdir(indir)
    svgs = [i for i in svgs if i.endswith(".svg")]
    svgs = sorted(svgs)
    ntotfiles = len(svgs)
    
    if lengths != None:
        l = pd.read_csv(lengths)
        samples = [i.replace(".svg", "") for i in svgs]
        l = l[l["sample"].isin(samples)]
        if col_thresh != None:
            l["short"] = l.get("short", np.where(l[col] <= float(col_thresh), True, False))
        else:
            l["col"] = l.get(col, l["short"])
        l = l.sort_values(col)
        tumour_tl = l.set_index("sample")[col].to_dict()

    else:
        l = pd.DataFrame({"sample": [i.replace(".svg", "") for i in svgs], "short": True})
        tumour_tl = {}
        for i in svgs:
            tumour_tl[i.replace(".svg", "")] = 1

    if f"{col}_conf" in l.columns: # if output from ml model, use confidence of boolean prediction to order
        l = l.sort_values(f"{col}_conf", ascending=False)
    shorts = l[l["short"] == True]["sample"].tolist()
    longs = l[l["short"] == False]["sample"].tolist()
    all_shorts = shorts.copy()
    all_longs = longs.copy()

    # Create subsampled groups
    short_set = []
    long_set = []
    nsubsets = 0
    nsamples = len(svgs)
    if subsets == 1:
        short_set.append(list({s: l for s, l in sorted(tumour_tl.items(), key=lambda item: item[1]) if s in shorts}.keys()))
        long_set.append(list({s: l for s, l in sorted(tumour_tl.items(), key=lambda item: item[1]) if s in longs}.keys()))
        nsubsets = 1
    else:
        for z in range(subsets):
            if nsamples > subsamples:
                short_samples = [shorts[i] for i in sorted(random.sample(range(len(shorts)), min(len(shorts), int(subsamples*split))))]
                long_samples = [longs[i] for i in sorted(random.sample(range(len(longs)), min(len(longs), int(subsamples*(1-split)))))] 
                short_set.append(list({s: l for s, l in sorted(tumour_tl.items(), key=lambda item: item[1]) if s in short_samples}.keys()))
                long_set.append(list({s: l for s, l in sorted(tumour_tl.items(), key=lambda item: item[1]) if s in long_samples}.keys()))
                shorts = list(set(shorts).difference(short_samples))
                longs = list(set(longs).difference(long_samples))
                nsamples = nsamples - len(short_samples) - len(long_samples)
            elif 0 < nsamples <= subsamples:
                short_set.append(shorts)
                long_set.append(longs)
                nsamples = 0

            nsubsets += 1

            if nsamples == 0:
                break
    
    task_arrange = progress.add_task("[red]Arranging svgs...", total=ntotfiles)
    
    for z in range(nsubsets):
        short_files = [i+".svg" for i in short_set[z]]
        long_files = [i+".svg" for i in long_set[z]]
        svgs = short_files + long_files
        original_size = 450
        nrows = ceil(len(svgs) / ncols)
        size = 300
        sf = size/original_size
        total_width = ncols*size
        total_height = (total_width/ncols)*nrows
        
        nfiles = len(svgs)
        nshort = len(short_files)
        nlong = len(long_files)
        srows = ceil(nshort/ncols)
        lrows = ceil(nlong/ncols)

        ## short (True) plotting (all if not ordered)
        short = sg.SVGFigure(f"{ncols*size}px", f"{srows*size}px")
        shorts = []
        r, c = 0, 0
        for i in svgs[:nshort]:
            tmp = sg.fromfile(os.path.join(indir, i))#.set_size("300pt")#.getroot()
            if col_thresh != None: # add order values to plot
                tmp.append(svgutils.compose.Text("%.2f"%round(float(l[l["sample"] == i.replace(".svg", "")][col].values[0]), 2), x=275, y=30, size=25, font="sans", weight="lighter"))

            tmp = tmp.getroot()#.set_size(300)
            shorts.append(tmp)
            tmp.moveto(r, c, scale_x=sf*1.25)
            r += size
            if r % total_width == 0:
                c += size
                r = 0
            progress.update(task_arrange, advance=1)
        #print(shorts)
        short.append(shorts)
        if lengths != None:
            short.save(os.path.join(outdir, f"short_{z}.svg"))
        else:
            Figure(f"{total_width}px", f"{total_height}px",
                    Panel(short.getroot())).save(os.path.join(outdir, f"{outname}_{z}.svg"))

        ## long (False) plotting (if order supplied)
        if lengths != None:
            long = sg.SVGFigure(f"{ncols*size}px", f"{lrows*size}px")
            longs = []
            r, c, lrows = 0, 0, 1
            for i in svgs[nshort:]:
                tmp = sg.fromfile(os.path.join(indir, i))
                if col_thresh != None: # add order values to plot
                    tmp.append(svgutils.compose.Text("%.2f"%round(float(l[l["sample"] == i.replace(".svg", "")][col].values[0]), 2), x=275, y=30, size=25, font="sans", weight="lighter"))

                tmp = tmp.getroot()

                longs.append(tmp)
                tmp.moveto(r, c, scale_x=sf*1.25)
                r += size
                if r % total_width == 0:
                    c += size
                    r = 0
                    lrows += 1
                progress.update(task_arrange, advance=1)

                long.append(longs)
                long.save(os.path.join(outdir, f"long_{z}.svg"))
                total_height = (total_width/ncols)*(srows+lrows)

                short = sg.fromfile(os.path.join(outdir, f"short_{z}.svg"))
                # stick long and short together with short on top
                Figure(f"{total_width}px", f"{total_height}px",
                    Panel(short.getroot()).move(0, 0),
                    Panel(long.getroot()).move(0, srows*size)).save(os.path.join(outdir, f"{outname}_{z}.svg"))

                a = sg.fromfile(os.path.join(outdir, f"{outname}_{z}.svg"))
                a.append(svgutils.compose.Line([(0, srows*size), (total_width, srows*size)], width=3, color="red"))
                a.save(os.path.join(outdir, f"{outname}_{z}.svg"))
        
        ## convert SVG to PNG
        task_convert = progress.add_task("[cyan]Converting to png...", total=None)
        cairosvg.svg2png(url=os.path.join(outdir, f"{outname}_{z}.svg"), write_to=os.path.join(outdir, f"{outname}_{z}.png"), scale=2)
        img = PIL.Image.open(os.path.join(outdir, f"{outname}_{z}.png"))
        wsize = int(float(img.size[0])/2)
        hsize = int(float(img.size[1])/2)
        progress.update(task_convert)
        ## better line visability
        task_resize = progress.add_task("[purple]Resizing...", total=None)
        img = img.resize((wsize, hsize), 1)
        img.save(os.path.join(outdir, f"{outname}_{z}.png"))
        progress.update(task_resize)


@click.command(context_settings={'show_default': True})
@click.argument("indir")
@click.option("-od", "--outdir", default="tmp_svg", help="Output directory")
@click.option("-oname", "--outname", default="all", help="Dataset name (output file name)")
@click.option("-n", "--ncols", default=8, help="Number of columns in combined output plot")
@click.option("--skip", is_flag=True, default=False, help="Skip making svgs, arrange existing")
## Ordering variables (optional)
@click.option("-l", "--lengths", default=None, help="Order csv (optional)")
@click.option("-c", "--col", default="short", help="Column in order dataframe to sort (continous) or split (boolean) by")
@click.option("--col-thresh", default=None, help="Threshold to split groups, required if order column is not boolean")
## Subsample variables
@click.option("--subsets", default=1, help="N subsets to plot")
@click.option("--subsamples", default=100, help="N samples per subset")
@click.option("--split", default=0.5, help="% of split group to plot per subset")
## Additional plot information
@click.option("--cytoband", type=str, default="hg19", help="Cytoband file")
@click.option("--cnp", default=None, help="Optional copy number segments")
## VCF filtering variables
@click.option("--var-prob", default="INS: 0.2, DEL: 0.3, INV: 0.15, DUP: 0.15, TRA: 0.4", type=DictParamType(), help="SV type prob thresholds")
@click.option("--var-su", default="INS: 9, DEL: 9, INV: 6, DUP: 6, TRA: 8", type=DictParamType(), help="SV type supporting reads thresholds")
@click.option("--var-size", default=500, type=int, help="Minimum SV size to plot (bp)")
def plot_all(indir, outdir, outname, ncols, skip, lengths, col, col_thresh, subsets, subsamples, split, cytoband, cnp, var_prob, var_su, var_size):
    ## create shortcut for dataset
    # if indir == "hawk":
    #     indir = "../hawk/filtered_uni_true/"
    #     lengths = "../hawk/all_lengths.csv"
    #     col = "stela"
    #     col_thresh = 3.81
    #     cytoband = "hg38"
    #     cnp = "../cnp/cnpinter/hawk_small/copyNumber_io/segmented.copynumber.csv"

    if cytoband in {"hg19", "hg38"}: # alias for included cytoband files
        cytoband = f"{cytoband}_cytoBand.txt"

    if lengths == None:
        print(f"{'='*10} In/Out vars {'='*10}\nVCF directory: {indir},\nOutput directory: {outdir},\nDataset name: {outname},\nNcols: {ncols}\n{'='*10} Circos vars {'='*10}\nCytoband file: {cytoband},\nCopy Number file: {cnp}")
    else:
        print(f"{'='*10} In/Out vars {'='*10}\nVCF directory: {indir},\nOutput directory: {outdir},\nDataset name: {outname},\nNcols: {ncols}\n{'='*10} Order vars {'='*10}\nOrder file: {lengths},\nSort/split column: {col},\nSort column split threshold: {col_thresh},\n{'='*10} Circos vars {'='*10}\nCytoband file: {cytoband},\nCopy Number file: {cnp}")
    
    files = os.listdir(indir)
    files = [str(i) for i in files if i.endswith(".vcf")]
    sample_names = [i.replace(".vcf", "") for i in files]

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

    print(f"{'='*10} Filter vars {'='*10}\nProbabilty thresholds: {prob_thresholds},\nSupporting reads thresholds: {su_thresholds},\nSV size trehsold: {var_size}")

    ## split circos plots by variable into two groups (optional)
    if lengths != None:
        l = pd.read_csv(lengths)    
        l["sample"] = l["sample"].astype(str)
        l = l[l["sample"].isin(sample_names)]
        if col not in set(list(l.columns)):
            print(f"ERROR: column {col} not in order df {lengths} (columns: {list(l.columns)})")
            exit
        l = l.sort_values(by=[col])
        l = l.reset_index()
        sample_names_sorted = l["sample"].tolist()
        files = [i+".vcf" for i in sample_names_sorted]
        nfiles = len(files)
    
    else:
        nfiles = len(files)

    if not os.path.exists(outdir):
        os.makedirs(outdir)

    with Progress() as progress:
        indiv = progress.add_task("[green]Making svgs...", total=nfiles)
        for i in files:
            if skip == True and os.path.isfile(os.path.join(outdir, i.replace(".vcf", "")+".svg")):
                continue
            c = read_file(os.path.join(indir, i), outdir, lengths, cytoband, cnp, prob_thresholds, su_thresholds, var_size)
            c.figure.savefig(os.path.join(outdir, i.replace(".vcf", "")+".svg"))
            progress.update(indiv, advance=1)

        arrange(outdir, outdir, outname, ncols, lengths, col, col_thresh, progress, subsets, subsamples, split)

    # for i in os.listdir(outdir): # prevents when rerunning code joining of two sets of plots
    #     if i == "done":
    #         continue
    #     shutil.move(os.path.join(outdir, i), os.path.join(outdir, "done", i))
    

if __name__ == "__main__":
    plot_all()
