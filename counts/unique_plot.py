#!/bin/python3

from collections import defaultdict
import sys
import os
import pickle
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.path import Path
import matplotlib.patches as patches
from matplotlib.collections import BrokenBarHCollection
from matplotlib import gridspec
import numpy as np
import pandas as pd
import glob
import pysam
import click

def load_vcfs(indir, pre):
    vcfs = glob.glob(os.path.join(indir, "*.vcf"))
    print("loading vcfs...")
    
    # Sample: svtype: chromosome: [(start, end) ... ]
    d = {key: defaultdict(lambda: defaultdict(list)) for key in [i for i in os.listdir(indir) if i.endswith("vcf")]}
    for vcf in vcfs:
        samp = vcf.split("/")[-1]#.split(".")[0]
        #sv_type = vcf.split("/")[-1].split(".")[2][:3]

        f = pysam.VariantFile(vcf)
        for l in f.fetch():
            sv_type = l.info["SVTYPE"]
            start = l.pos
            #end = l.info["END"]
            end = l.stop
            if pre == True:
                chr1 = l.contig
                chr2 = l.info["CHR2"]
            else:
                chr1 = "chr"+l.contig
                chr2 = "chr"+l.info["CHR2"]


            if chr1 == chr2:
                d[samp][sv_type][chr1].append((start, end))
            else:
                d[samp][sv_type][chr1].append((start, start + 1))  # Translocation will appear as lines
                d[samp][sv_type][chr2].append((end, end + 1))
    return d

def gen_indexes(l):
    p = 0
    t = []
    for i in l:
        t.append((p, i+p))
        p += i+1
    print(t)

def set_sizes(reference):
    if reference == "hg19":
        chrom_lengths = (249250621,243199373,198022430,191154276,180915260,171115067,159138663,146364022,
                         141213431,135534747,135006516,133851895,115169878,107349540,102531392,90354753,
                         81195210,78077248,59128983,63025520,48129895,51304566,155270560,59373566)
        chrom_indexes = ((0, 249250621), (249250622, 492449995), (492449996, 690472426), (690472427, 881626703), 
                        (881626704, 1062541964), (1062541965, 1233657032), (1233657033, 1392795696), (1392795697, 1539159719), 
                        (1539159720, 1680373151), (1680373152, 1815907899), (1815907900, 1950914416), (1950914417, 2084766312), 
                        (2084766313, 2199936191), (2199936192, 2307285732), (2307285733, 2409817125), (2409817126, 2500171879), 
                        (2500171880, 2581367090), (2581367091, 2659444339), (2659444340, 2718573323), (2718573324, 2781598844), 
                        (2781598845, 2829728740), (2829728741, 2881033307), (2881033308, 3036303868), (3036303869, 3095677435))
        # chrom_indexes = ((0, 249250621), (249250622, 492449994), (492449995, 690472424),
        #                 (690472425, 881626700), (881626701, 1062541960), (1062541961, 1233657027),
        #                 (1233657028, 1392795690), (1392795691, 1539159712),
        #                 (1539159713, 1680373143), (1680373144, 1815907890),
        #                 (1815907891, 1950914406), (1950914407, 2084766301),
        #                 (2084766302, 2199936179L), (2199936180L, 2307285719L),
        #                 (2307285720L, 2409817111L), (2409817112L, 2500171864L),
        #                 (2500171865L, 2581367074L), (2581367075L, 2659444322L), (2659444323L, 2718573305L),
        #                 (2718573306L, 2781598825L), (2781598826L, 2829728720L), (2829728721L, 2881033286L),
        #                 (2881033287L, 3036303846L), (3036303847L, 3095677412L))
        order = ('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14',
                'chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
    elif reference == "hg38":
        chrom_lengths = (248956422,242193529,198295559,190214555,181538259,170805979,159345973,145138636,138394717,
                         133979422,135086622,133275309,114364328,107043718,101991189,90338345,83257441,80373285,
                         58617616,64444167,46709983,50818468,156040895,57227415)
        chrom_indexes = ((0, 248956422), (248956423, 491149952), (491149953, 689445512), (689445513, 879660068), 
                         (879660069, 1061198328), (1061198329, 1232004308), (1232004309, 1391350282), 
                         (1391350283, 1536488919), (1536488920, 1674883637), (1674883638, 1808863060), 
                         (1808863061, 1943949683), (1943949684, 2077224993), (2077224994, 2191589322), 
                         (2191589323, 2298633041), (2298633042, 2400624231), (2400624232, 2490962577), 
                         (2490962578, 2574220019), (2574220020, 2654593305), (2654593306, 2713210922), 
                         (2713210923, 2777655090), (2777655091, 2824365074), (2824365075, 2875183543), 
                         (2875183544, 3031224439), (3031224440, 3088451855))
        order = ('chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14',
                'chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY')
    else:
        print("add lengths for this reference")

    return chrom_lengths, chrom_indexes, order

def chromosome_collections(df, y_positions=10, height=9):
    """
    Yields BrokenBarHCollection of features that can be added to an Axes
    object.
    Mostly stolen from: https://www.biostars.org/p/147364/#147637
    df : Must at least have columns ['chrom', 'start', 'end', 'color']. If no
        column 'width', it will be calculated from start/end.
    y_positions : Keys are chromosomes, values are y-value at which to anchor the
                    BrokenBarHCollection
    height : Height of each BrokenBarHCollection
    """
    for chrom, group in df.groupby('chrom'):
        xranges = group[['start', 'width']].values
        colors = group['colors']
        print(xranges)
        sys.exit()
        yield xranges, colors

def make_ideogram(inf, order):
    # Mostly stolen from: https://www.biostars.org/p/147364/#147637
    # Path to hg19 cytoband table from ucsc genome browser
    ideo = pd.read_csv(inf, sep="\t")
    ideo.columns = ["chrom", "start", "end", "name", "gieStain"]

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



def plot_SVs(bed_dict, name, order, chrom_lengths, chrom_indexes, output_dir, cytoband, show_plot):
    print("Plotting: ", name)

    # Two dicts
    ideograms, ideo_colors = make_ideogram(cytoband, order)

    cc = 0
    for ch in order:
        fig = plt.figure(figsize=(12,2))
        gs = gridspec.GridSpec(2, 1, height_ratios=[20, 1])
        axarr = [plt.subplot(gs[0]), plt.subplot(gs[1])]
        fig.subplots_adjust(hspace=0.1cyto)
        #plt.tight_layout()

        axarr[0].axes.set_xticks([])
        axarr[0].axes.set_yticks([])
        axarr[0].axis('off')
        axarr[0].set_title("{} {}".format(name.replace(".vcf", ""), ch), weight="bold")
        axarr[0].set_xlim(0, chrom_lengths[cc])
        axarr[0].set_ylim(1, 2.1)

        # Plot each SV type
        #color=iter(cm.rainbow(np.linspace(0,1,4)))
        color=iter(["#ff0040", "#00ff00", "#0080ff", "#c800fa"])


        vv = [1, 1, 1, 1]  # This is the y value to draw SVtype at
        #SVtypes = ["DEL", "INV", "DUP:TANDEM"]  # Insertions too small for pindel
        SVtypes = ["DEL", "INV", "DUP", "TRA"]

        leg_vals = []

        for idx, SV in enumerate(SVtypes):
            start_end = bed_dict[SV][ch]  # [(start, end), ..]
            if len(start_end) == 0:
                next(color)
                continue

            # First sort start_end values, in case start > end value
            srt = [sorted(ss) for ss in start_end]

            curve_height = 0.5
            clr=next(color)
            # Make a quadratic curve between start and end, midpoint = end - start / 2 + start
            for item in srt:

                vertical = vv[idx]
                start = item[0]
                if start < 0:
                    continue    # FIND THIS BUG
                end = item[1]

                x = np.linspace(-1, 1, 100)
                y = (-x*x + 2) * curve_height
                if SV=="INV" or SV=="TRA":
                    # Curve goes upside down
                    y = (x*x) * curve_height

                factor = (end - start)/2
                x = (x + 1) * factor
                x += start
                y += vertical
                # For legend plotting, change to shorter version
                if SV=="DUP:TANDEM":
                    SV="DUP"

                axarr[0].plot(x, y, color=clr, lw=1, label=SV if SV not in leg_vals else "")

                if SV not in leg_vals:
                    leg_vals.append(SV)

                # Add dotted horizontal line
                x = [0, chrom_lengths[cc]]
                y = [1.5, 1.5]
                axarr[0].plot(x, y, '--', dashes=(2,3), color="grey")

        axarr[0].legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize=10, frameon=False)

        # Now add ideogram below main figure:
        idg = ideograms[ch]
        idc = ideo_colors[ch]

        axarr[1].broken_barh(idg, (0, 1), facecolors=idc)

        axarr[1].set_xlabel("Genomic position (bp)", weight="bold" )
        axarr[1].set_xlim(0, chrom_lengths[cc])
        axarr[1].set_ylim(0, 1)
        axarr[1].axes.set_yticks([])

        out_folder = os.path.join(output_dir, os.path.basename(name).replace(".vcf", ""))

        if not os.path.exists(out_folder):
            os.makedirs(out_folder)

        if show_plot == False:
            plt.savefig("{}_{}.pdf".format(os.path.join(out_folder, name.replace(".vcf", "")), ch), bbox_inches='tight', transparent=True)
        
            cc += 1
            plt.close()
        else:
            plt.show()
            cc += 1


@click.command()
@click.argument("indir")
@click.option("-r", "--reference", type=click.Choice(["hg19", "hg38", "t2t"]), default="hg19", 
        help="reference genome type")
@click.option("-c", "--cytoband", default="/home/alex/Desktop/uni/PhD/variants/processing_vcfs/hg19_cytoBand.txt",
        help="path to cytoband file")
@click.option("-o", "--output", default="unique_circos",
        help="output directory for plots")
@click.option("-s", "--show-plot", is_flag=True, default=False,
        help="show or save")
@click.option("-p", "--pre", default=True, is_flag=True,
        help="prefix remove chr from contig names")
def plot(indir, reference, cytoband, output, show_plot, pre):
    d = load_vcfs(indir, pre)
    chrom_lengths, chrom_indexes, order = set_sizes(reference)
    for sample in d:
        plot_SVs(d[sample], sample, order, chrom_lengths, chrom_indexes, output, cytoband, show_plot) 


if __name__ == "__main__":
    plot()
