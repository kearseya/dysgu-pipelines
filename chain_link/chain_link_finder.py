"""
In the style of the original chain-finder algoritm (see chromoplexy paper https://www.sciencedirect.com/science/article/pii/S0092867413003437), group breakpoints into sets of chains.
A break point belongs to a chain if the distance to the next nearest breakpoint is below a binomial threshold.
"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import collections
from collections import defaultdict, Counter
from scipy import stats
from sklearn.neighbors import KDTree
import numpy as np
import glob
import pysam
import networkx as nx
from matplotlib.pyplot import cm
import matplotlib as mpl
import os
import re
import shutil
from subprocess import call
import math
import copy
from kneed import KneeLocator

import random
import pprint
from matplotlib import ticker
from networkx.utils import dict_to_numpy_array
from itertools import product
import seaborn as sns
from matplotlib.colors import LogNorm, SymLogNorm, ListedColormap

import patchworklib as pw
from math import ceil
import cairosvg
import PIL
import svgutils
import svgutils.transform as sg
from svgutils.compose import *

import pycircos
import click

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

@click.group(name="pipe", invoke_without_command=True)
@click.argument("indir")
@click.option("-o", "--outdir", default="chain_out")
@click.option("-p", "--prob", default=0.05)
@click.option("-b", "--breaks", default=20)
@click.option("-l", "--lower-p", default=0.001)
@click.option("-u", "--upper-p", default=0.9)
@click.option("-n", "--np-tests", default=20)
@click.option("-t", "--tl", default="sample_pairs.csv")
@click.option("-s", "--sep-thresh", default=3.81)
@click.option("-e", "--sample-col", default="tumor_db")
@click.option("-a", "--stela-col", default="tumor_stela")
@click.option("-c", "--cytofile", default="hg38_cytoBand.txt", type=str)
@click.option("--chrom-len", default=None)
@click.option("--stack-unit", default="telomere length (kb)", hidden=True)
@click.option("--size-thresh", default=50000, help="SV size filter")
@click.option("--var-prob", default="INS: 0.2, DEL: 0.3, INV: 0.15, DUP: 0.15, TRA: 0.4", type=DictParamType(), help="SV type prob thresholds")
# @click.option("--var-su", default="INS: 9, DEL: 9, INV: 6, DUP: 6, TRA: 8", type=DictParamType(), help="SV type supporting reads thresholds")
@click.pass_context
def pipe(ctx, indir, outdir, prob, breaks, lower_p, upper_p, np_tests, tl, sep_thresh, sample_col, stela_col, cytofile, chrom_len, stack_unit, size_thresh, var_prob):
    ctx.ensure_object(dict)
    ctx.obj["indir"] = indir
    ctx.obj["outdir"] = outdir
    ctx.obj["pval"] = prob
    ctx.obj["bt"] = breaks
    ctx.obj["lp"] = lower_p
    ctx.obj["up"] = upper_p
    ctx.obj["np"] = np_tests
    ctx.obj["tlfn"] = tl
    ctx.obj["sep_thresh"] = sep_thresh
    ctx.obj["sample_col"] = sample_col
    ctx.obj["stela_col"] = stela_col
    ctx.obj["cyto"] = cytofile
    ctx.obj["chrom_len"] = chrom_len
    ctx.obj["stack_unit"] = stack_unit
    ctx.obj["size_thresh"] = size_thresh
    ctx.obj["var_prob"] = var_prob
    if ctx.invoked_subcommand is None:
        if indir == None:
            click.echo(ctx.get_help())
            ctx.exit()
            pass
        else:
            click.echo("Running pipeline")
            ctx.invoke(run)
            ctx.invoke(plot_results)
    else:
        pass

@pipe.command(name="run")
@click.option("-p", "--prob", default=0.04)
@click.option("-b", "--breaks", default=8)
@click.option("-l", "--chrom-len", default=None)
@click.option("--plot", is_flag=True, default=True)
@click.pass_context
def run(ctx, prob, breaks, chrom_len, plot):
    PROB = prob
    BREAKS_THRESH = breaks
    if ctx.obj["pval"] != 0.05:
        PROB = ctx.obj["pval"]
    if ctx.obj["bt"] != 8:
        BREAKS_THRESH = ctx.obj["bt"]
    indir = ctx.obj["indir"]
    OUTDIR = ctx.obj["outdir"]
    PLOT_CIRCOS = plot
    if ctx.obj["chrom_len"] != None:
        chrom_len = ctx.obj["chrom_len"]

    ## allow filtering of vcf without creating new vcf file
    def_prob_thresholds = {"INS": 0.2, "DEL": 0.3, "INV": 0.15, "DUP": 0.15, "TRA": 0.4}
    #def_su_thresholds = {"INS": 9, "DEL": 9, "INV": 6, "DUP": 6, "TRA": 8}

    if ctx.obj["var_prob"] == None: # use default
        prob_thresholds = def_prob_thresholds
    else:
        prob_thresholds = ctx.obj["var_prob"]
        for i in ["INS", "DEL", "INV", "DUP", "TRA"]: # if SV type not defined, use default
            if i not in prob_thresholds.keys():
                prob_thresholds[i] = def_prob_thresholds[i]
    # if ctx.obj["var_su"] == None:
    #     su_threhsolds = def_su_thresholds
    # else:
    #     su_thresholds = ctx.obj["var_su"]
    #     for i in ["INS", "DEL", "INV", "DUP", "TRA"]:
    #         if i not in su_thresholds.keys():
    #             su_thresholds[i] = def_su_thresholds[i]
    
    if 0 < PROB < 1:
        print("*"*20, "\nProbability: ", PROB, " Break threshold:" , BREAKS_THRESH)
    else:
        print("*"*20, "\nLow prob: ", ctx.obj["lp"], ", High prob: ", ctx.obj["up"], ", N probs: ", ctx.obj["np"], ", Break threshold:" , BREAKS_THRESH)

    OUTDIR = os.path.join(OUTDIR, "found_chains_p{}_t{}".format(PROB, BREAKS_THRESH))

    # Save the randomized output in a different directory
    # OUTDIR = "found_false_positive_chains"
    # vcfs = glob.glob("randomized_breakpoints/*.vcf")

    PLOT_CIRCOS = True

    np.random.seed(12345678)

    if chrom_len == None:
        vcfs = list(glob.glob(f"{indir}/*.vcf"))
        chrom_len = vcfs[0]

    if chrom_len[-4:] == ".bam":
        bamf = pysam.AlignmentFile(chrom_len, "rb")
        head = bamf.header["SQ"]
        cl = {}  # Get the lengths of individual chromosomes from the bam header file
        for item in head:
            cl[item["SN"]] = float(item["LN"])
    elif chrom_len[-4:] == ".csv":
        cldf = pd.read_csv(chrom_len)
        cldf.index = cldf["chromosome"]
        cl = cldf["length"].to_dict()
    elif chrom_len[-4:] == ".tsv":
        cldf = pd.read_csv(chrom_len, sep="\t")
        cldf.index = cldf["chromosome"]
        cl = cldf["length"].to_dict()
    elif chrom_len[-4:] == ".vcf":
        tmp = pysam.VariantFile(chrom_len)
        head = tmp.header
        cl = {}
        for c in head.contigs:
            for line in str(head).split("\n"):
                if line.startswith(f"##contig=<ID={c},"):
                    cl[c] = int(re.search(r"length=(\d*)", line).group(1))

    # Use for adding color to chomosome nodes
    color_nodes = {chrom: i for chrom, i in zip(cl.keys(), list(cm.rainbow(np.linspace(0, 1, len(cl)))))}

    # Use exact breakpoints
    # df = pd.read_csv("../Microhomology_mapper/microhomology_insertions.all_from_vcf.csv")

    # Group by sample, get a list of all breaks into tuple (chr1, break1, chr2, break2)
    # Make a sorted list of all the breaks on each chromosome for each sample
    all_breaks = defaultdict(list)  # sample: list of breaks
    breaks_list = defaultdict(lambda: defaultdict(list))  # Sample: Chromosome: list of breakpoints

    if os.path.isdir(indir):
        count = 0
        vcfs = glob.glob(f"{indir}/*.vcf")
        if len(vcfs) > 0 :
            print("vcfs: ", vcfs)
            for v in vcfs:
                #print(v)
                f = pysam.VariantFile(v)
                s = os.path.basename(v).replace(".vcf", "")
                for line in f.fetch():
                    #if list(line.filter.keys())[0] != "lowProb" or float(line.samples[s]["PROB"]) > 0.2 or float(line.samples[s]["SU"]) > 4: 
                    if float(line.samples[s]["PROB"]) >= prob_thresholds[line.info["SVTYPE"]]:
                        if (line.info["SVTYPE"] != "TRA" and line.stop-line.pos >= ctx.obj["size_thresh"]) or line.info["SVTYPE"] == "TRA":
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

        beds = glob.glob(f"{indir}/*.bed")
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

        df = pd.DataFrame(columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2"])
        for s in all_breaks.keys():
            #print(pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "chrom2", "start2"]))
            df = pd.concat([df, pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2"])])
        print(df)

    elif os.path.isfile(indir):
        count = 0
        if indir[-3:] == "bed":
            df = pd.read_csv(indir, sep="\t", index_col=None) 
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

        elif indir[-3:] == "vcf":
            f = pysam.VariantFile(indir)
            count = 0
            s = os.path.basename(v).replace(".vcf", "")
            for line in f.fetch():
                #if list(line.filter.keys())[0] != "lowProb" or float(line.samples[s]["PROB"]) > 0.2 or float(line.samples[s]["SU"]) > 4: 
                if float(line.samples[s]["PROB"]) >= prob_thresholds[line.info["SVTYPE"]]:
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

    else:
        assert os.path.isfile(indir), "input not valid file or directory"

    all_breaks_copy = copy.deepcopy(all_breaks)
    print("Total SVs", count)


    # Calculate frequency of breaksites per chromosome accross all samples, needed for p-value calculation
    ps = {str(k): 0.0 for k in cl.keys()}  # Sample: avg_breaks per bp across all samples
    for key, value in breaks_list.items():
        break_n_per_chrom = defaultdict(int)

        for chrom, b in value.items():  # For each break in the sample
            break_n_per_chrom[str(chrom)] += len(b)
        avg = {str(k): v / cl[str(k)] for k, v in break_n_per_chrom.items()}
        for k, v in break_n_per_chrom.items():
            ps[k] += v


    ps = {k: (v/len(breaks_list) / cl[k]) for k, v in ps.items()}

    # Make a nearest neighbour tree for all chromosomes of samples
    nn = defaultdict(lambda: defaultdict(list))  # Sample: chromosome: neighbour tree
    for samp, value in breaks_list.items():
        for chrom, bs in value.items():
            X = np.array(bs)[:, np.newaxis]
            tree = KDTree(X)
            nn[samp][chrom] = tree  # Use dist, ind = tree.query([[1.5]], k=3)


    def get_pval(fails, prob_success):
        """Finds probability that 1 or more breaksites are observed within this given distance.
           The distance in bp is the number of 'fails'. i.e. the number of bp without a break
            Success is 0, the probability that no breaks are oberved."""
        pval = 1 - stats.binom.pmf(0, fails, prob_success)

        # Comparing the chromoplexy paper method, they use essentially a two sided binomial pmf:
        # print(pval, (1 - math.pow(1 - 2*prob_success, fails)))
        return pval


    """ Loop over the list of breakpoints for each sample; For each breakpoint look at up- and downstream sites,
    if the p-value is less than threshold, make an edge on the graph """


    def find_connected(tree, data, chrom, pos, p_val):

        dist, ind = tree.query([[pos]], k=len(data))  # k=len(data) to enumerate all neighbours
        pvals = [get_pval(i, ps[chrom]) for i in dist[0]]
        dist, ind = dist.flatten(), ind.flatten()

        c = set([])
        for i, p in enumerate(pvals[1:]):  # Skip first entry which is the current SV
            if p < p_val:
                # Add 1 to i becuase first value is skipped
                c.add(("{}:{}".format(chrom, int(pos)), "{}:{}".format(chrom, int(data[ind[i + 1]]))))

        return c


    def neigh(s, samp, nn, prob):
        
        tree1 = nn[samp][s[0]]  # Tree for chromosome 1
        tree2 = nn[samp][s[3]]  # Tree for chromosome 2, (intra or inter)

        data1 = np.array(tree1.data).flatten()
        data2 = np.array(tree2.data).flatten()

        # e.g. s = ('chr2', 33092006, 'chr2', 141363341)
        # Do for first break
        c1 = find_connected(tree1, data1, s[0], s[1], prob)
        c2 = find_connected(tree2, data2, s[3], s[5], prob)

        connected = set([])
        connected |= c1  # Set union
        connected |= c2

        return connected


    def get_clusters(n, clr, count, samp):

        cluster_edges = [i for i, j in zip(n, clr) if j == "lightgrey"]  # Extract the edges which highlight clusters

        # Make a new graph to isolate the clusters
        G = nx.Graph()
        G.add_edges_from(cluster_edges)

        with open("{}/{}/{}_{}_clusters.csv".format(OUTDIR, samp, samp, count), "w") as out:
            out.write("Chrom\tBreak_number\tSize(bp)\tMedian_spacing(bp)\tBreak_pos\n")
            for k in [G.subgraph(c).copy() for c in nx.connected_components(G)]: #nx.connected_component_subgraphs(G):  # This gives each cluster in the chain

                # Sort the nodes by position, then get the start and end of each cluster
                cluster = sorted([(i.split(":")[0], int(i.split(":")[1])) for i in k.nodes()], key=lambda x: x[1])

                c_size = cluster[-1][1] - cluster[0][1]
                success = len(cluster)  # k
                faliure = c_size - success  # r
                prob_success = ps[cluster[0][0]]  # Get the probability associated with that chromosome (p)

                neg_b = stats.nbinom.pmf(k=faliure, n=success, p=prob_success)

                clust_size = cluster[-1][1] - cluster[0][1]
                clust_len = len(cluster)
                spacing = np.array([cluster[i+1][1] - cluster[i][1] for i in range(len(cluster)-1)])
                median_spacing = np.median(spacing)
                pos = ",".join([str(i[1]) for i in cluster])

                chrom = cluster[0][0]

                out.write("{c}\t{n}\t{s}\t{m}\t{p}\t{nb}\n".format(c=chrom, n=clust_len, s=clust_size, m=median_spacing,
                                                                   p=pos, nb=neg_b))


    def arrange_circos(indir, outdir, outname, lengths, real, prob_v, break_t):
        print("arranging circos plots...")
        if outdir == None:
            outdir = indir
        outname = os.path.join(outdir, f"{outname}_p{prob_v}_t{break_t}")
        svgs = os.listdir(indir)
        svgs = [i for i in svgs if i.endswith(".svg")]
        svgs = sorted(svgs)
        original_size = 900
        #nshort = 16
        ncols = 8
        #scols = nshort/ncols
        nrows = ceil(len(svgs) / ncols)
        size = 300
        sf = size/original_size
        total_width = ncols*size
        total_height = (total_width/ncols)*nrows


        l = pd.read_csv(lengths)
        l = l.rename({ctx.obj["sample_col"]: "sample", ctx.obj["stela_col"]: "stela"}, axis=1)

        samples = [i.replace(".svg", "") for i in svgs]
        l = l[l["sample"].isin(samples)]
        if real == True:
            l["short"] = l.get("short", np.where(l["stela"] <= ctx.obj["sep_thresh"], True, False))
        else:
            l["stela"] = l.get("stela", l["short"])
        l = l.sort_values("stela")#, ascending=False)

        if "short_conf" in l.columns:
            l = l.sort_values("short_conf")
        short_samples = l[l["short"] == True]["sample"].tolist()
        scols = ceil(len(short_samples)/ncols)
        long_samples = l[l["short"] == False]["sample"].tolist()
        short_files = [i+".svg" for i in short_samples]
        long_files = [i+".svg" for i in long_samples]
        svgs = short_files + long_files
        nfiles = len(svgs)
        nshort = len(short_files)

        # short
        short = sg.SVGFigure(f"{ncols*size}px", f"{scols*size}px")
        shorts = []
        r, c = 0, 0
        for i in svgs[:nshort]:
            tmp = sg.fromfile(os.path.join(indir, i))#.set_size("300pt")#.getroot()
            if real == True:
                tmp.append(svgutils.compose.Text("%.2f"%round(float(l[l["sample"] == i.replace(".svg", "")]["stela"].values[0]), 2), x=275, y=30, size=25, font="sans", weight="lighter"))

            #tmp.set_size(f"{50000}px") # doesnt work
            tmp = tmp.getroot()#.set_size(300)
            shorts.append(tmp)
            tmp.moveto(r, c, scale_x=sf*1.25)
            r += size
            if r % total_width == 0:
                c += size
                r = 0
        #print(shorts)
        short.append(shorts)
        short.save(os.path.join(outdir, f"short_p{prob_v}_t{break_t}.svg"))

        # long
        long = sg.SVGFigure(f"{ncols*size}px", f"{((ncols-nshort)/ncols)*size}px")
        longs = []
        r, c = 0, 0
        for i in svgs[nshort:]:
            tmp = sg.fromfile(os.path.join(indir, i))
            if real == True:
                tmp.append(svgutils.compose.Text("%.2f"%round(float(l[l["sample"] == i.replace(".svg", "")]["stela"].values[0]), 2), x=275, y=30, size=25, font="sans", weight="lighter"))

            tmp = tmp.getroot()

            longs.append(tmp)
            tmp.moveto(r, c, scale_x=sf*1.25)
            r += size
            if r % total_width == 0:
                c += size
                r = 0

        #print(longs)
        long.append(longs)
        long.save(os.path.join(outdir, f"long_p{prob_v}_t{break_t}.svg"))

        Figure(f"{total_width}px", f"{total_height}px",
            Panel(short.getroot()),
            Panel(long.getroot()).move(0, scols*size)).save(f"{outname}.svg")

        a = sg.fromfile(f"{outname}.svg")
        a.append(svgutils.compose.Line([(0, scols*size), (total_width, scols*size)], width=3, color="red"))
        a.save(f"{outname}.svg")

        cairosvg.svg2png(url=f"{outname}.svg", write_to=f"{outname}.png", scale=2)
        img = PIL.Image.open(f"{outname}.png")
        wsize = int(float(img.size[0])/2)
        hsize = int(float(img.size[1])/2)
        img = img.resize((wsize, hsize), 1)
        img.save(f"{outname}.png")

    def find_neighbours(sample, probability, plot_circos):
        # Clear old results
        if os.path.exists(OUTDIR + "/" + sample):
            shutil.rmtree(OUTDIR + "/" + sample)
        # Make a folder to save results in
        if not os.path.exists(OUTDIR + "/" + sample):
            os.makedirs(OUTDIR + "/" + sample)

        real_edges = set([])
        links = set([])

        for s in all_breaks[sample]:
            sv = {("{}:{}".format(s[0], int(s[1])), "{}:{}".format(s[3], int(s[5])))}
            real_edges |= sv
            grey_link = neigh(s, sample, nn, probability)
            links |= grey_link

        G = nx.Graph()
        G.add_edges_from(links, color="lightgrey")
        G.add_edges_from(real_edges, color="black")

        # Add some color to nodes for different chromosomes
        for node in G.nodes():
            chrom = node.split(":")[0]
            cchrom = mpl.colors.rgb2hex(color_nodes[chrom])
            G.nodes[node]["cchrom"] = cchrom
            G.nodes[node]["label"] = chrom[3:]

        # Save out a graph in GraphML format
        nx.write_graphml(G, "{}/{}/{}_graph.graphml".format(OUTDIR, sample, sample))

        sub_graph_sizes = sorted([(len(i.nodes()), i) for i in [G.subgraph(c).copy() for c in nx.connected_components(G)]], #list(nx.connected_component_subgraphs(G))],
                                 key=lambda x: x[0], reverse=True)

        #
        nodes_in_chains = set([])  # Keep a list of all the breaks which end up in chains
        svs_in_chains = dict([])  # set([])  # List of the black edges
        count = 0

        color_iter = iter(["red", "blue", "green", "purple", "yellow", "orange", "darkred", "lightgreen", "black", "lightblue",
                           "darkblue", "wheat", "peru", "teal", "cyan", "yellowgreen", "greenyellow", "olive", "gold",
                           "rosybrown", "lime", "navy", "pink", "hotpink", "plum", "royalblue", "azure", "khaki",
                           "gold", "coral"]*100)
        for size, g in sub_graph_sizes:
            if size > BREAKS_THRESH:

                clr = next(color_iter)  # The color can be used as the chain-id

                nodes_in_chains |= set(g.nodes())
                svs_in_chains.update({k: clr for k, v in nx.get_edge_attributes(g, 'color').items() if v == "black"})
                
                if plot_circos == True:
                    # Make a figure of the chain structure
                    fig = plt.figure()
                    pos = nx.spring_layout(g, k=0.1, iterations=30)  # Increase spacign of nodes, so they dont overlap
                    plt.title(sample + ", " + str(size), weight="bold")

                    edges, colors = zip(*nx.get_edge_attributes(g, 'color').items())

                    get_clusters(edges, colors, count, sample)  # Save a copy of all the clusters within the chain
                    node, labels = zip(*nx.get_node_attributes(g, 'label').items())
                    labs = {k: v for k, v in zip(node, labels)}
                    nodes_colors = {i: g.degree(i) for i in g.nodes()}

                    nx.draw_networkx_edges(g, pos, edgelist=edges, edge_color=colors, width=1, node_size=50, alpha=0.8)  # font_size=12
                    nx.draw_networkx_nodes(g, pos, nodelist=nodes_colors.keys(), node_color=list(nodes_colors.values()), cmap='rainbow', node_size=200, alpha=0.9, linewidths=1)  # font_size=12,
                    nx.draw_networkx_labels(g, pos, labs,  font_size=12, alpha=0.4)  # label_pos=0.3,

                    plt.axis('off')

                    ax1 = fig.add_axes([0.88, 0.22, 0.03, 0.6])
                    cmap = mpl.cm.rainbow
                    norm = mpl.colors.Normalize(vmin=min(nodes_colors.values()), vmax=max(nodes_colors.values()))
                    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm, orientation='vertical', ticks=range(0, max(set(nodes_colors.values())) + 1, 2))
                    cb1.set_label('Breakpoint degree')
                    plt.savefig("{}/{}/{}_{}_{}.pdf".format(OUTDIR, sample, sample, count, size))
                    plt.close()

                count += 1  # Chain found, keep count

                # Get the number of nodes in chains or not in chains
                edges, colors = zip(*nx.get_edge_attributes(g, 'color').items())

                chain = set([])
                chain_rev = set([])
                for n in sorted(zip(edges, colors), key=lambda x: x[0][0]):
                    if n[1] == 'black':  # True edge not a 'closeness edge'

                        n1, n2 = n[0][0].split(":"), n[0][1].split(":")

                        bn = (n1[0], int(n1[1]), n2[0], int(n2[1]))
                        bn_rev = (n2[0], int(n2[1]), n1[0], int(n1[1]))

                        chain.add(bn)
                        chain_rev.add(bn_rev)  # The edge may be the other way around

                #print(chain)
                # Save a copy of the chain in circos friendly format
                with open("{}/{}/{}_{}_{}_breakpoints.csv".format(OUTDIR, sample, sample, count, size), "w") as out:
                    for c in sorted(chain, key=lambda x: (x[0], x[1])):
                        out.write("{chrom} {start1} {end1} {chrom2} {start2} {end2} z=0,color={color}\n".format(chrom=c[0],
                                                                                                            start1=c[1],
                                                                                                            start2=c[1],
                                                                                                            chrom2=c[2],
                                                                                                            end1=c[3],
                                                                                                            end2=c[3],
                                                                                                            color=clr))

        nodes_in_chains = [i.split(":") for i in nodes_in_chains]
        nodes_in_chains = set([(i, int(j)) for i, j in nodes_in_chains])  # Convert into usable format

        not_in_chains = []
        result = []
        #print(svs_in_chains)
        for sv in all_breaks[sample]:
            #print(sv)
            d1 = ("{}:{}".format(sv[0], int(sv[1])), "{}:{}".format(sv[3], int(sv[5])))
            d2 = ("{}:{}".format(sv[3], int(sv[5])), "{}:{}".format(sv[0], int(sv[1])))

            if d1 not in svs_in_chains and d2 not in svs_in_chains:
                #print("non-chained")
                not_in_chains.append(sv)
                result.append((False, ""))
                sv.append(False)
                sv.append("")
            else:
                if d1 in svs_in_chains:
                    #print("d1")
                    result.append((True, svs_in_chains[d1]))
                    sv.append(True)
                    sv.append(svs_in_chains[d1])
                else:
                    #print("d2")
                    result.append((True, svs_in_chains[d2]))
                    sv.append(True)
                    sv.append(svs_in_chains[d2])
        # Sanity check
        # print(len(all_breaks[sample]), len(not_in_chains), len(svs_in_chains))
        # assert len(all_breaks[sample]) - len(not_in_chains) == len(svs_in_chains)
        print("Breaks in chains", len(nodes_in_chains))

        with open("{}/{}/{}_non-chained.csv".format(OUTDIR, sample, sample), "w") as out:
            for c in sorted(not_in_chains, key=lambda x: (x[0], x[1])):
                out.write("{chrom} {start1} {end1} {chrom2} {start2} {end2} z=0,color=grey\n".format(chrom=c[0],
                                                                                                 start1=c[1],
                                                                                                 start2=c[4],
                                                                                                 chrom2=c[3],
                                                                                                 end1=c[2],
                                                                                                 end2=c[5]))
        if plot_circos:
            if os.path.exists(f"{OUTDIR}/circos") == False:
                os.makedirs(f"{OUTDIR}/circos")
            color_dict   = {"gneg":"#FFFFFF00", "gpos25":"#EEEEEE", "gpos50":"#BBBBBB", "gpos75":"#777777", "gpos100":"#000000",
                            "gvar":"#FFFFFF00", "stalk":"#C01E27", "acen":"#D82322"}
            values_all   = []
            arcdata_dict = collections.defaultdict(dict)

            var_prob_threshold = 0.4
            low_size = 800
            up_size = 850

            aof = 45
            dof = 44
            circle = pycircos.Gcircle(figsize=(10, 10))
            nl = {}
            for i in cl:
                nl[i.replace("chr", "")] = cl[i]
           
            ## create chromosomes
            order_of_chr = [f"{i}" for i in range (1, 23)]
            order_of_chr.append("X")
            order_of_chr.append("Y")
            for i in order_of_chr:
                name    = i
                length  = nl[i]
                arc     = pycircos.Garc(arc_id=name, size=length, interspace=2, raxis_range=(low_size, up_size), 
                               labelposition=100, label_visible=True, facecolor="white", labelsize=17)
                circle.add_garc(arc)
            circle.set_garcs()

            ## plot cytoband
            arcdata_dict = collections.defaultdict(dict)
            with open(ctx.obj["cyto"]) as f:
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

            ## plot not in links - first so they get covered
            tra_list = []
            color = "grey"
            for var in not_in_chains:
                #print(var)
                chr1, start1, end1, chr2, start2, end2, chained, color = var
                color = "grey"
                chained = False
                name = chr1.replace("chr", "")
                name2 = chr2.replace("chr", "")
                if chr1 == chr2 and name in order_of_chr:
                    source = (name, int(start1), int(start1), low_size-aof)
                    destination = (name, int(end1), int(end1), low_size-aof)
                    #print(source, destination)
                    circle.chord_plot(source, destination, facecolor=color, edgecolor=color, linewidth=0.1)#facecolor=type_cols[sv_type], edgecolor=type_cols[sv_type], linewidth=0.4)
                elif chr1 != chr2 and name in order_of_chr and name2 in order_of_chr: 
                    #start2 = start #l.info["CHR2_POS"]
                    #end2 = start2
                    tra_list.append(((name, start1, end1, low_size-aof), (name2, end2, end2, low_size-aof)))
                else:
                    continue
            for t in tra_list:
                circle.chord_plot(t[0], t[1], facecolor="grey", linewidth=0.1)

            ## plot links
            tra_list = []
            breakpoint_files = glob.glob(f"{OUTDIR}/{sample}/*_breakpoints.csv")
            for bfp in breakpoint_files:
                with open(bfp) as f:
                    for line in f:
                        #print(line)
                        chr1, start1, end1, chr2, start2, end2, color = line.split()
                        name = chr1.replace("chr", "")
                        name2 = chr2.replace("chr", "")
                        color = color.split(",")[1].split("=")[1]
                        if chr1 == chr2 and name in order_of_chr:
                            source = (name, int(start1), int(start2), low_size-aof)
                            destination = (name, int(end1), int(end2), low_size-aof)
                            #print(source, destination)
                            circle.chord_plot(source, destination, facecolor=color, edgecolor=color, linewidth=0.1)#facecolor=type_cols[sv_type], edgecolor=type_cols[sv_type], linewidth=0.4)
                        elif chr1 != chr2 and name in order_of_chr and name2 in order_of_chr: 
                            #end2 = start2
                            tra_list.append(((name, int(start1), int(start2), low_size-aof), (name2, int(end1), int(end2), low_size-aof), color))
                        else:
                            continue
                    for t in tra_list:
                        circle.chord_plot(t[0], t[1], facecolor=t[2], edgecolor=t[2], linewidth=0.1)
            
            circle.figure.savefig(f"{OUTDIR}/circos/{sample}.svg")
       
        #print("FIND NEIGHBOURS FUNCTION")
        #print(all_breaks[sample])
        return len(nodes_in_chains) / 2, count, len(all_breaks[sample]), result


    def process_data(probability):

        res = []
        results = []
        for sample, dx in df.groupby("Sample"):
            print("Finding for", sample)
            nic, nc, nb, chained_or_not = find_neighbours(sample, probability, plot)
            res.append((sample, nic, nc, nb))

            assert len(chained_or_not) == len(dx)
            results += chained_or_not
    
        return res, results


    if 0 < PROB < 1:
        # Process once with fixed pval
        r, results = process_data(PROB)

        df = pd.DataFrame(columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "chained", "color"])
        for s in all_breaks.keys():
            df = pd.concat([df, pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "chained", "color"])])
        print(df)
        
        arrange_circos(os.path.join(OUTDIR, "circos"), os.path.join(OUTDIR, "circos"), "all", ctx.obj["tlfn"], True, PROB, BREAKS_THRESH)
        #print( len(results), len(df))
        df["pval"] = PROB
        #df["chained"] = [i[0] for i in results]
        #df["color"] = [i[1] for i in results]
        df.to_csv(os.path.join(OUTDIR, "all_svs.unique.chains.csv"), index=None)
        plot_chained_vs_non_chained(OUTDIR, ctx.obj["tlfn"], ctx.obj["sample_col"], ctx.obj["stela_col"], ctx.obj["sep_thresh"], ctx.obj["stack_unit"])

    else:
        # Process with a range of p-values and collect the results
        data = []
        knee_x_ticks = []
        chained_p_vals = []
        if os.path.exists(os.path.join(ctx.obj["outdir"], f"t{BREAKS_THRESH}_p{ctx.obj['lp']}-{ctx.obj['up']}")) == False:
            os.makedirs(os.path.join(ctx.obj["outdir"], f"t{BREAKS_THRESH}_p{ctx.obj['lp']}-{ctx.obj['up']}"))
        
        for prob_v in np.linspace(ctx.obj["lp"], ctx.obj["up"], ctx.obj["np"]):
            all_breaks = copy.deepcopy(all_breaks_copy)

            OUTDIR = os.path.join(ctx.obj["outdir"], f"t{BREAKS_THRESH}_p{ctx.obj['lp']}-{ctx.obj['up']}", f"found_chains_p{prob_v}_t{BREAKS_THRESH}")
            tmp_df = df.copy()
            print(prob_v, "-"*80)
            res, rr = process_data(prob_v)
            if plot == True:
                arrange_circos(os.path.join(OUTDIR, "circos"), os.path.join(OUTDIR, "circos"), "all", ctx.obj["tlfn"], True, prob_v, BREAKS_THRESH)

            df = pd.DataFrame(columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "chained", "color"])
            for s in all_breaks.keys():
                df = pd.concat([df, pd.DataFrame([[s]+list(x) for x in all_breaks[s]], columns=["Sample", "chrom1", "start1", "end1", "chrom2", "start2", "end2", "chained", "color"])])
            df["pval"] = prob_v
            df.to_csv(os.path.join(OUTDIR, "all_svs.unique.chains.csv"), index=None)

            #for sample, nic, nc, nb in res:
            #    data.append({"sample": sample, "pval": prob_v, "Chained_breaks": nic, "Total_breaks": nb, "Chains": nc})
            #tmp_df["pval"] = prob_v
            #tmp_df["chained"] = [i[0] for i in rr]
            #tmp_df["color"] = [i[1] for i in rr]
            #tmp_df.to_csv(os.path.join(OUTDIR, "all_svs.unique.chains.csv"), index=None)
            chained_stat = plot_chained_vs_non_chained(OUTDIR, ctx.obj["tlfn"], ctx.obj["sample_col"], ctx.obj["stela_col"], ctx.obj["sep_thresh"], ctx.obj["stack_unit"])
            chained_p_vals.append(chained_stat.pvalue)
            knee_x_ticks.append(prob_v)

        if ctx.obj["np"] > 2:
            plot_elbow(knee_x_ticks, chained_p_vals, os.path.join(OUTDIR, ".."))
        pd.DataFrame.from_records(data).to_csv(os.path.join(OUTDIR, "..", "variable_pval.csv"), index=None)



### For plotting ###
def blockmodel(G, partitions):
    """Stolen from NetworkX 1.9. Returns a reduced graph constructed using the generalized block modeling
    technique. Faster than the quotient graph method, but not by much.

    The blockmodel technique collapses nodes into blocks based on a
    given partitioning of the node set.  Each partition of nodes
    (block) is represented as a single node in the reduced graph.

    Edges between nodes in the block graph are added according to the
    edges in the original graph.  If the parameter multigraph is False
    (the default) a single edge is added with a weight equal to the
    sum of the edge weights between nodes in the original graph
    The default is a weight of 1 if weights are not specified.  If the
    parameter multigraph is True then multiple edges are added each
    with the edge data from the original graph.

    Parameters
    ----------
    G : graph
        A networkx Graph or DiGraph
    partitions : list of lists, or list of sets
        The partition of the nodes.  Must be non-overlapping.
    multigraph : bool, optional
        If True return a MultiGraph with the edge data of the original
        graph applied to each corresponding edge in the new graph.
        If False return a Graph with the sum of the edge weights, or a
        count of the edges if the original graph is unweighted.

    Returns
    -------
    blockmodel : a Networkx graph object

    Examples
    --------
    >>> G=nx.path_graph(6)
    >>> partition=[[0,1],[2,3],[4,5]]
    >>> M=nx.blockmodel(G,partition)

    References
    ----------
    .. [1] Patrick Doreian, Vladimir Batagelj, and Anuska Ferligoj
    	"Generalized Blockmodeling",Cambridge University Press, 2004.
    """
    # Create sets of node partitions, frozenset makes them hashable
    part = [frozenset(i) for i in partitions]

    # Check for overlapping node partitions
    u = set()
    for p1, p2 in zip(part[:-1], part[1:]):
        u.update(p1)
        if len(u.intersection(p2)) > 0:
            raise nx.NetworkXException("Overlapping node partitions.")

    # Initialize blockmodel graph
    M = nx.Graph()

    # Add nodes and properties to blockmodel
    # The blockmodel nodes are node-induced subgraphs of G
    # Label them with integers starting at 0
    for i, p in zip(range(len(part)), part):
        M.add_node(p)
        # The node-induced subgraph is stored as the node 'graph' attribute
        SG = G.subgraph(p)
        M.nodes[p]['graph'] = SG

    # Create mapping between original node labels and new blockmodel node labels
    block_mapping = {}
    for n in M:
        nodes_in_block = M.nodes[n]['graph'].nodes()
        block_mapping.update(dict.fromkeys(nodes_in_block, n))

    # Add edges to block graph
    for u, v, d in G.edges(data=True):

        bmu = block_mapping[u]
        bmv = block_mapping[v]
        if bmu == bmv:  # no self loops
            continue

        # For graphs and digraphs add single weighted edge
        weight = d.get('weight', 1.0)  # default to 1 if no weight specified
        if M.has_edge(bmu, bmv):
            M[bmu][bmv]['weight'] += weight
        else:
            M.add_edge(bmu, bmv, weight=weight)
    return M


def analyse_clutser_distribution(g):

    np.random.seed(1000)
    random.seed(101)
    samp = g.split("/")[1]
    G = nx.read_graphml(g)

    # Get biggest subgraph, the chain
    g = sorted(nx.connected_component_subgraphs(G), key=lambda x: len(x.nodes()), reverse=True)[0]

    ori = g.copy()
    clrs = nx.get_edge_attributes(ori, "color")  # Used to find black or grey edges
    be = [k for k, v in clrs.items() if v == "black"]
    # Now remove black edges to get the clusters
    e = nx.get_edge_attributes(g, "color")
    g.remove_edges_from([i for i, clr in e.items() if clr == "black"])

    # Now find the ranges for the blockmodel partition
    partition = {index: sorted([i for i in i.nodes()])
                 for index, i in
                 enumerate(sorted([g.subgraph(c) for c in nx.connected_components(g)], key=lambda x: len(x.nodes()), reverse=True))}

    # print(be)
    # For each partition get the number of intra cluster SVs
    intra_links = {k: [] for k in partition.keys()}
    for p, s in partition.items():
        for sv in s:
            edge = [i for i in be if sv in i]
            if edge:
                if edge[0][0] in s and edge[0][1] in s:
                    intra_links[p].append(edge[0])
    # print(intra_links)
    counts = {k: len(v) for k, v in intra_links.items()}
    counts_p_sizes = {k: (len(v), len(partition[k])) for k, v in intra_links.items()}
    # print(counts_p_sizes)

    c = Counter(counts.values())
    c = {k: (v / (float(sum(c.values())) / len(partition[k]))) for k, v in c.items()}  # The CDF
    # print(c)
    plt.figure()
    # plt.ylim(0, 1), plt.xlim(-0.05, 15)
    plt.scatter(c.keys(), c.values())

    plt.figure()
    y, x = zip(*counts_p_sizes.values())

    plt.scatter(x, y)
    # plt.show()

    partpos = {k: [] for k in partition.keys()}
    dall = []
    for k, v in partition.items():

        d = []
        for i in range(len(v)-1):
            d.append(v[i+1] - v[i])
        # print(d)
        # print(v - v.min())
        dall += d

    plt.hist(sorted(dall), bins=100)
    #plt.yscale("log")
    # plt.show()

    # Now plot the distribution
    return


def generate_block_edge_matrix(partition, ori, M):
    """Need to find the adjacency matrix for the block model, also the number of internal edges for each
    partition needs to be added to this."""

    adj = nx.to_numpy_array(M)

    # Go through nodes and check if any link to within this partition (intra) from original graph
    clrs = nx.get_edge_attributes(ori, "color")
    black_edges = [k for k, v in clrs.items() if v == "black"]
    internal_edges = []
    gb = nx.Graph()
    gb.add_edges_from(black_edges)
    for index, p in partition.items():  # Sorted by number of nodes within partition
        internal_edges.append(len([i for i in gb.edges(p) if (i[0] in p and i[1] in p)]))

    for index, v in enumerate(internal_edges):
        adj[index, index] = v

    return adj.astype(float)


def calc_modularity(G, communities, weight='weight'):
    """
    https://networkx.github.io/documentation/stable/_modules/networkx/algorithms/community/quality.html
    Returns the modularity of the given partition of the graph.

    Modularity is defined in [1]_ as

    .. math::

        Q = \frac{1}{2m} \sum_{ij} \left( A_{ij} - \frac{k_ik_j}{2m}\right)
            \delta(c_i,c_j)

    where $m$ is the number of edges, $A$ is the adjacency matrix of
    `G`, $k_i$ is the degree of $i$ and $\delta(c_i, c_j)$
    is 1 if $i$ and $j$ are in the same community and 0 otherwise.

    Parameters
    ----------
    G : NetworkX Graph

    communities : list
        List of sets of nodes of `G` representing a partition of the
        nodes.

    Returns
    -------
    Q : float
        The modularity of the paritition.

    Raises
    ------
    NotAPartition
        If `communities` is not a partition of the nodes of `G`.

    Examples
    --------
    >>> G = nx.barbell_graph(3, 0)
    >>> nx.algorithms.community.modularity(G, [{0, 1, 2}, {3, 4, 5}])
    0.35714285714285704

    References
    ----------
    .. [1] M. E. J. Newman *Networks: An Introduction*, page 224.
       Oxford University Press, 2011.

    """

    multigraph = G.is_multigraph()
    directed = G.is_directed()
    m = G.size(weight=weight)
    if directed:
        out_degree = dict(G.out_degree(weight=weight))
        in_degree = dict(G.in_degree(weight=weight))
        norm = 1 / m
    else:
        out_degree = dict(G.degree(weight=weight))
        in_degree = out_degree
        norm = 1 / (2 * m)

    def val(u, v):
        try:
            if multigraph:
                w = sum(d.get(weight, 1) for k, d in G[u][v].items())
            else:
                w = G[u][v].get(weight, 1)
        except KeyError:
            w = 0
        # Double count self-loops if the graph is undirected.
        if u == v and not directed:
            w *= 2
        return w - in_degree[u] * out_degree[v] * norm

    Q = sum(val(u, v) for c in communities for u, v in product(c, repeat=2))
    return Q * norm


def modularity(g, partition, samp, random):
    Q = calc_modularity(g, partition.values())
    return {"samp": samp, "Q": Q, "random": random, 'nodes': len(g.nodes())}


def calc_prob_mat(mat, ens):
    prob = np.zeros_like(mat)
    for i in range(len(prob)):
        for j in range(len(prob)):
            cnt = Counter(ens[i, j, :])
            x, y = zip(*[(k, ky / float(sum(cnt.values()))) for k, ky in cnt.items()])
            probs = dict(zip(x, y))
            if int(mat[i, j]) in probs:
                prob[i, j] = probs[int(mat[i, j])]  # This is the pmf, the probability of observing exactly x links
            else:
                prob[i, j] = 1.0 / len(ens[0, 0, :])  # The number of simulations
    return prob

def community_link_prob(mat, partition, ens, samp, y_labels, outdir):

    prob = calc_prob_mat(mat, ens)

    from matplotlib.colors import LogNorm
    fig = plt.figure()
    ax = plt.subplot()
    plt.subplots_adjust(left=0.2, right=0.8)
    plt.rcParams.update({'font.size': 16})
    # plt.title("Block incidence permutation-test, {}".format(samp))
    mask = np.zeros_like(prob)
    mask[np.triu_indices_from(mask)] = True
    mask = np.roll(mask, 1, axis=1)
    mask[:, 0] = 0
    for a in ['top', 'bottom', 'left', 'right']:
        ax.spines[a].set_visible(False)

    im = plt.imshow(np.ma.array(prob, mask=mask), interpolation="none", cmap="GnBu_r", norm=LogNorm(vmin=prob.min(), vmax=prob.max()))
    plt.xticks(np.arange(0, len(mat), 1.0))
    plt.xticks(np.arange(0, len(mat)), [len(partition[i]) for i in range(len(partition.keys()))], rotation=90)
    plt.yticks(np.arange(0, len(mat)), y_labels)

    cb = plt.colorbar(fraction=0.026, pad=0.04, label="Probability")

    plt.clim(1e-3, 0.05)
    # Only plot for clusters with > 2 breaks
    lim = len([i for i in range(len(partition.keys())) if len(partition[i]) > 2]) - 1

    plt.xlim(-0.5, lim + 0.5), plt.ylim(lim + 0.5, -0.5)
    plt.savefig(outdir + "/{}/{}_block_prob_matrix_pval.pdf".format(samp.split("_")[0], samp))


def calc_heterogeneity_index(y_labels, mat, ens):

    def calc_het(y_labels, mat, ens):
        prob = calc_prob_mat(mat, ens)
        val = [0]
        norm_y = np.array(y_labels).astype(float) / sum(y_labels)
        for i in range(len(mat)):
            for j in range(len(mat)):
                pval = prob[i, j]
                nodes = norm_y[i]  # min((norm_y[i], norm_y[j]))  # The weight of nodes between these clusters
                v = nodes / pval
                #if v > 1:
                val.append(-2 * (np.log(pval))) # - (np.log(nodes))))
                if i == j:
                    break
        return sum(val)

    datav = calc_het(y_labels, mat, ens)
    print("Data h", datav)
    # Calculate a distribution from random graphs
    randv = [calc_het(y_labels, ens[:, :, i], ens) for i in range(100)] #len(ens[0, 0, :]))]
    print("Rand h", np.array(randv).mean())

    # Number of random values greater than the observed data
    s = len([i for i in randv if i > datav])
    prob = 1.
    if s > 0:
        prob = s / float(len(randv))  # The fraction of simulations greater than the data == probability
    else:
        prob = 1. / len(randv)
    return prob


def random_black_edges_graph(ori_edges, ori_nodes, grey_edges):
    # Sometimes there are additional edges, due to limited SV resolution, some SV calls appear to share the same
    # breakpoint - this can happen when short fragments are involved and so does not violate the infinite sites
    # model.
    ori_nodes = list(ori_nodes)
    G2 = nx.Graph()
    G2.add_nodes_from(ori_nodes)
    nodes = list(G2.nodes())
    random.shuffle(nodes)
    if len(nodes) % 2 > 0:
        nodes = nodes[:-1]  # Must be even
    nit = iter(nodes)
    # No self edges, only one edge per pair
    random_edges = [(i, next(nit)) for i in nit]
    G2.add_edges_from(random_edges, color="black")
    num_to_add = len(ori_edges) - len(G2.edges())
    # Pick random nodes to add edges between
    while num_to_add > 0:
        u, v = random.choice(ori_nodes), random.choice(ori_nodes)
        if u != v and not G2.has_edge(u, v):
            G2.add_edge(u, v, color="black")
            num_to_add -= 1
    return G2  # Return anyway if a connected graph cant be made, should happen very rarely


def plot_block_matrix(ori, block_matrix_black, partition, samp, samp_idx, outdir):
    # Generate a list of random block model arrays to get a p-value for each partition in terms of the number of
    # nodes involved in SVs
    ensemble_blocks = []
    ori_nodes = ori.nodes()
    ori_edges = [(u, v) for u, v, d in ori.edges(data=True) if d["color"] == "black"]
    grey_edges = [(u, v, d) for u, v, d in ori.edges(data=True) if d["color"] == "lightgrey"]
    obs = 500  # 1000
    for i in range(obs):
        G2 = random_black_edges_graph(ori_edges, ori_nodes, grey_edges)
        M2 = blockmodel(G2, [partition[v] for v in sorted(partition.keys())])
        bm_black = generate_block_edge_matrix(partition, G2, M2)
        ensemble_blocks.append(bm_black)
    ens = np.dstack(ensemble_blocks)

    # Make some row labels for the graphs
    row_labels = []
    for k in range(len(partition)):
        row = sorted([i for i in partition[k]], key=lambda x: (x.split(":")[0], int(x.split(":")[1])))[0].split(":")
        row[1] = str(int(round(float(row[1]) / 1e6)))
        # row = '${}_{{{}}}$'.format(row[0][3:].replace("000", "").replace("n_gl", ""), row[1])
        row = '{}:{}'.format(row[0][3:].replace("000", "").replace("n_gl", ""), row[1])
        row_labels.append(row)

    y_labels = [len(partition[i]) for i in range(len(partition.keys()))]

    community_link_prob(block_matrix_black, partition, ens, samp, y_labels, outdir)
    community_link_prob(ens[:, :, 1], partition, ens, samp + "_rand", y_labels, outdir)
    print("-"*20)

    h_value = calc_heterogeneity_index(y_labels, block_matrix_black, ens)
    print(h_value)

    fig = plt.figure()
    ax = plt.subplot()
    plt.subplots_adjust(left=0.2, right=0.8)
    # plt.title("Block incidence for SVs, {}".format(samp))
    plt.rcParams.update({'font.size': 16})

    mask = np.zeros_like(ensemble_blocks[0])
    mask[np.triu_indices_from(mask)] = True
    mask = np.roll(mask, 1, axis=1)
    mask[:, 0] = 0
    cmap = "magma_r"
    im = plt.imshow(np.ma.array(block_matrix_black, mask=mask), interpolation="none", cmap=cmap)
    plt.xticks(np.arange(0, len(block_matrix_black), 1.0))
    plt.xticks(np.arange(0, len(block_matrix_black)), [len(partition[i]) for i in range(len(partition.keys()))], rotation=90)
    plt.yticks(np.arange(0, len(block_matrix_black)), y_labels)
    cb = plt.colorbar(im, fraction=0.026, pad=0.04)
    for a in ['top', 'bottom', 'left', 'right']:
        ax.spines[a].set_visible(False)
    cb.set_label('Breakpoint adjacency')
    tick_locator = ticker.MaxNLocator(nbins=np.max(block_matrix_black) if block_matrix_black.max() < 10 else 10, integer=True, min_n_ticks=2)
    cb.locator = tick_locator
    cb.update_ticks()

    lim = len([i for i in range(len(partition.keys())) if len(partition[i]) > 2]) - 1
    plt.xlim(-0.5, lim + 0.5), plt.ylim(lim + 0.5, -0.5)

    plt.savefig(outdir + "/{}/{}_{}_block_prob_matrix_black.pdf".format(samp, samp, samp_idx))
    #plt.show()
    plt.close()


    fig = plt.figure()
    ax = plt.subplot()
    plt.subplots_adjust(left=0.2, right=0.8)
    # plt.title("Block matrix random, {}".format(samp))
    plt.rcParams.update({'font.size': 16})
    im = plt.imshow(np.ma.array(ensemble_blocks[0], mask=mask), interpolation="none", cmap=cmap)
    for a in ['top', 'bottom', 'left', 'right']:
        ax.spines[a].set_visible(False)
    plt.xticks(np.arange(0, len(block_matrix_black), 1.0))
    plt.xticks(np.arange(0, len(block_matrix_black)), [len(partition[i]) for i in range(len(partition.keys()))], rotation=90)
    plt.yticks(np.arange(0, len(block_matrix_black)), y_labels)
    cb = plt.colorbar(im, fraction=0.026, pad=0.04)

    cb.set_label('Breakpoint adjacency')
    tick_locator = ticker.MaxNLocator(nbins=np.max(block_matrix_black) if block_matrix_black.max() < 10 else 10, integer=True, min_n_ticks=2)
    plt.clim(0, np.max(block_matrix_black))
    cb.locator = tick_locator
    cb.update_ticks()

    lim = len([i for i in range(len(partition.keys())) if len(partition[i]) > 2]) - 1
    plt.xlim(-0.5, lim + 0.5), plt.ylim(lim + 0.5, -0.5)

    plt.savefig(outdir + "/{}/{}_{}_block_prob_matrix_mean_ensemble.pdf".format(samp, samp, samp_idx))
    plt.close()

    return G2  # Return a random graph


def plot_block_model_nodes(M, M_rand, random_graph, partition, ori, samp_idx, samp, outdir):
    # Generate a block model plot ----->
    clrs = nx.get_edge_attributes(ori, "color")  # Used to find black or grey edges

    # For each partition, find how many true SVs link to within the cluster. The name is the "chromosome: start in Mb"
    node_names_intra_cluster_links = {}  # Partiton index: [name, intra edges, partition]
    node_names_intra_cluster_links_random = {}  # Use G2 for random, should be in scope from above
    node_labels = {}  # Index: name
    random_clrs = nx.get_edge_attributes(random_graph, "color")

    for index, p in partition.items():  # Sorted by number of nodes within partition
        p = list(p)
        name = '${}_{{{}}}$'.format(p[0].split(":")[0][3:].replace("000", "").replace("n_gl", ""),
                                    int(round(float(p[0].split(":")[1])) / 1e6))
        node_labels[index] = name

        # Go through nodes and check if any link to within this partition (intra) from original/random graphs
        # Only use black edges, these are the TRUE SVs
        connecting_edges = len([i for i in ori.edges(p) if i[0] in p and i[1] in p and
                                ((i in clrs and clrs[i] == "black") or (i[::-1] in clrs and clrs[i[
                                                                                                 ::-1]] == "black"))])  # The edge may be the wrong way round so need to check both 'directions'
        node_names_intra_cluster_links[index] = [name, connecting_edges]

        # Do for random graph
        connecting_edges_rand = len([i for i in random_graph.edges(p) if i[0] in p and i[1] in p and (
        (i in random_clrs and random_clrs[i] == "black") or (
        i[::-1] in random_clrs and random_clrs[i[::-1]] == "black"))])
        node_names_intra_cluster_links_random[index] = [name, connecting_edges_rand]

    M.colors, M_rand.colors = {}, {}  # The degree of each node will become the color
    M.edge_labels, M_rand.edge_labels = {}, {}  # The weight of each edge
    M.intra = {k: 100+v[1]*100 for k, v in node_names_intra_cluster_links.items()}  # Controls size of node
    M.colors = {i: M.degree(i) for i in M.nodes()}

    M_rand.intra = {k: 100+v[1]*100 for k, v in node_names_intra_cluster_links_random.items()}
    M_rand.colors = {i: M_rand.degree(i) for i in M_rand.nodes()}

    # Find out maximum and minimum color values
    highest = max([max(set(M.colors.values())), max(set(M_rand.colors.values()))])
    lowest = min([min(set(M.colors.values())), min(set(M_rand.colors.values()))])
    print("Lowest, Highest", lowest, highest)
    # Convert to RGBarray
    cmap = mpl.cm.GnBu
    norm = mpl.colors.Normalize(vmin=min(M.colors.values()), vmax=max(M.colors.values())) #clip=True)
    # Plot for actual data
    fig = plt.figure(figsize=(4,3))
    plt.subplots_adjust(right=0.7)
    plt.title("{}".format(samp))
    params = {'mathtext.default': 'regular'}
    plt.rcParams.update(params)
    pos = nx.spring_layout(M, iterations=20)
    weights = [M[u][v]['weight']*1.5 for u, v in M.edges()]

    nx.draw(M, pos=pos, node_size=list(M.intra.values()), node_color=list(M.colors.values()), cmap=cmap, #norm=norm,
            width=weights, linewidths=1, vmin=lowest, vmax=highest)

    ax = plt.gca()  # to get the current axis
    ax.collections[0].set_edgecolor("black")

    nx.draw_networkx_edges(M, pos=pos, edge_color="grey", width=weights)

    pos_higher = {}
    y_off = 0.0  # offset on the y axis
    x_off = 0.02
    for k, v in pos.items():
        pos_higher[k] = (v[0]+x_off, v[1]+y_off)

    #nx.draw_networkx_labels(M, pos=pos_higher, labels=node_labels, font_color="black", font_size=12)

    ax1 = fig.add_axes([0.78, 0.22, 0.03, 0.6])

    norm2 = mpl.colors.Normalize(vmin=lowest, vmax=highest)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm2, orientation='vertical', ticks=range(0, highest + 1, 2))
    cb1.set_label('Cluster degree')

    plt.savefig(outdir + "/{}/{}_{}_block_model_largest_subgraph.pdf".format(samp, samp, samp_idx), transparent=True)


    # Plot for random graph --->
    norm = mpl.colors.Normalize(vmin=min(M_rand.colors.values()), vmax=max(M_rand.colors.values()))

    fig = plt.figure(figsize=(4, 3))
    plt.subplots_adjust(right=0.7)

    plt.title("Random SV edges {}".format(samp))
    params = {'mathtext.default': 'regular'}
    plt.rcParams.update(params)
    pos = nx.spring_layout(M_rand, iterations=20)
    weights = [M_rand[u][v]['weight']*1.5 for u, v in M_rand.edges()]

    nx.draw(M_rand, pos=pos, node_size=list(M_rand.intra.values()), node_color=list(M_rand.colors.values()), cmap=cmap,
            #norm=norm, 
            width=weights, linewidths=1)
    ax = plt.gca()  # to get the current axis
    ax.collections[0].set_edgecolor("black")

    nx.draw_networkx_edges(M_rand, pos=pos, edge_color="grey", width=weights)

    pos_higher = {}
    y_off = 0.0  # offset on the y axis
    x_off = 0.02
    for k, v in pos.items():
        pos_higher[k] = (v[0]+x_off, v[1]+y_off)

    # nx.draw_networkx_labels(M_rand, pos=pos_higher, labels=node_labels, font_color="black", font_size=12)
    ax1 = fig.add_axes([0.78, 0.22, 0.03, 0.6])

    norm2 = mpl.colors.Normalize(vmin=lowest, vmax=highest)
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap=cmap, norm=norm2, orientation='vertical', ticks=range(0, highest + 1, 2))

    cb1.set_label('Cluster degree')

    plt.savefig(outdir + "/{}/{}_{}_block_model_random_sv_edges.pdf".format(samp, samp, samp_idx), transparent=True)
    plt.close()


def modularity_density(partitions, black_edges):
    partitions = map(set, partitions.values())
    Gq = nx.Graph()
    Gq.add_edges_from(black_edges)
    Qds = 0
    m = float(len(Gq.edges()))
    for idx, part in enumerate(partitions):
        nc = float(len(part))
        # Find out how many self-community links and non-self links
        mc = 0.
        ec = 0.
        for i in part:
            # Get neighbours
            if not Gq.has_node(i):
                nc -= 1
                continue
            neigh = Gq.neighbors(i)
            for j in neigh:
                if j in part:
                    mc += 1
                else:
                    ec += 1
        if nc == 1:
            pc = 0.
        elif nc == 0:  # No nodes left in partition
            continue
        else:
            pc = 2.*mc / (nc * (nc - 1))

        first_sum = ((mc / m) * pc) - ( ((2 * mc + ec) / (2 * m)) * pc)**2

        # Get other partitions
        other_com = [p for ii, p in enumerate(partitions) if ii != idx]

        inner_sum = 0.
        for other in other_com:
            nc2 = len(other)
            # Get number of links between current and other
            mc_c2 = 0.
            for node in other:
                if not Gq.has_node(node):
                    nc2 -= 1
                    continue
                mc_c2 += len([1 for neigh in Gq.neighbors(node) if neigh in part])

            if nc2 == 0 or mc_c2 == 0:
                continue

            inner_sum += mc_c2**2 / (2. * m * nc * nc2)

        Qds += (first_sum - inner_sum)
    return Qds


def process_chain(g, samp, samp_idx, outdir):
    ori = g.copy()
    g = nx.Graph(g)
    # Now remove black edges to get the clusters
    e = nx.get_edge_attributes(g, "color")
    black_edges = [i for i, clr in e.items() if clr == "black"]
    grey_edges = [i for i, clr in e.items() if clr == "lightgrey"]
    g.remove_edges_from(black_edges)
    # Now find the ranges for the blockmodel partition
    partition = {index: i.nodes() for index, i in enumerate(sorted([g.subgraph(c) for c in nx.connected_components(g)], key=lambda x:
                 len(x.nodes()), reverse=True))}
    print("partition: ", partition)
    # Cant anaylse modularity with only 1 partition
    if len(partition) < 3:
        return []
    print(samp)
    #
    M = blockmodel(ori, [partition[v] for v in sorted(partition.keys())])
    # Save number of nodes in each partition
    # with open(in_folder + "/{}/{}_partition.txt".format(samp, samp), "w") as pout:
    #    pout.write(",".join([str(len(i)) for i in sorted(partition.values(), key=lambda x: len(x), reverse=True)]))
    print("Number of breakpoints, partitions", len(ori.nodes()), len(partition))

    block_matrix_black = generate_block_edge_matrix(partition, ori, M)
    random_graph = plot_block_matrix(ori, block_matrix_black, partition, samp, samp_idx, outdir)

    black_edges_random = [(u, v) for u, v, d in random_graph.edges(data=True) if d["color"] == "black"]

    Qds = modularity_density(partition, black_edges)
    Qds_rand = modularity_density(partition, black_edges_random)

    M_rand = blockmodel(random_graph, [partition[v] for v in sorted(partition.keys())])

    plot_block_model_nodes(M, M_rand, random_graph, partition, ori, samp_idx, samp, outdir)

    partition_key = {n: key for key, val in partition.items() for n in val}
    # print("black_edges: ", black_edges)
    # print("partition_key: ", partition_key)

    nx.set_node_attributes(random_graph, partition_key, "partition")

    # Show node assoritvity
    black_only = nx.Graph()
    black_only.add_edges_from(black_edges)
    nx.set_node_attributes(black_only, partition_key, "partition")
    # print("black_only: ", black_only)

    # Keep only clusters with more than 3 nodes, otherwise not informative

    def calc_assrt(graph):
        md = nx.attribute_mixing_dict(graph, "partition")
        # print("md: ", md)
        if None in md:  # remove None (not in any partition)
            del md[None]
        #md = {k: v for k, v in md.items() if sum(v.values()) > -1}
        mm = dict_to_numpy_array(md)
        # print("mm: ", mm)
        #mm = mm / mm.sum()
        #if len(mm) < 1:  # Only 1 cluster
        #    return False
        return nx.algorithms.assortativity.correlation.attribute_ac(mm)

    assort = calc_assrt(black_only)
    assort_rand = calc_assrt(random_graph)

    records = [modularity(black_only, partition, samp, False),
               modularity(random_graph, partition, samp, True)]

    #if not assort or not assort_rand or assort is np.nan or assort_rand is np.nan:
    #    return []

    # Keep only largest connected component for diameter
    random_graph.add_edges_from(grey_edges)
    # Rerely random graph is disconnected
    lbefore = len(random_graph.nodes())
    random_graph = sorted(list([random_graph.subgraph(c) for c in nx.connected_components(random_graph)]), key=lambda x: len(x.nodes()), reverse=True)[0]
    print("Number of nodes missing from random", lbefore - len(random_graph.nodes()))

    lbefore = len(M_rand.nodes())
    M_rand = sorted(list([M_rand.subgraph(c) for c in nx.connected_components(M_rand)]), key=lambda x: len(x.nodes()), reverse=True)[0]
    print("Number of nodes missing from block Model", lbefore - len(M_rand.nodes()))

    return [{"Sample": samp,
             "edges": len(black_edges),
             "edges": len(black_edges_random),
             "Q_data": records[0]["Q"],
             "Q_rand": records[1]["Q"],
             "Assortativity_data": assort,
             "Assortativity_rand": assort_rand,
             "Qds_data": Qds,
             "Qds_rand": Qds_rand,
             "Diameter_data": nx.diameter(ori),
             "Diameter_rand": nx.diameter(random_graph),
             "samp_idx": samp_idx,
             }]


def label_propagation_communities(G):
    """
    https://networkx.github.io/documentation/stable/_modules/networkx/algorithms/community/label_propagation.html
    Generates community sets determined by label propagation

    Finds communities in `G` using a semi-synchronous label propagation
    method[1]_. This method combines the advantages of both the synchronous
    and asynchronous models. Not implemented for directed graphs.

    Parameters
    ----------
    G : graph
        An undirected NetworkX graph.

    Yields
    ------
    communities : generator
        Yields sets of the nodes in each community.

    Raises
    ------
    NetworkXNotImplemented
       If the graph is directed

    References
    ----------
    .. [1] Cordasco, G., & Gargano, L. (2010, December). Community detection
       via semi-synchronous label propagation algorithms. In Business
       Applications of Social Network Analysis (BASNA), 2010 IEEE International
       Workshop on (pp. 1-8). IEEE.
    """
    coloring = _color_network(G)
    # Create a unique label for each node in the graph
    labeling = {v: k for k, v in enumerate(G)}
    while not _labeling_complete(labeling, G):
        # Update the labels of every node with the same color.
        for color, nodes in coloring.items():
            for n in nodes:
                _update_label(n, labeling, G)

    for label in set(labeling.values()):
        yield set((x for x in labeling if labeling[x] == label))


def _color_network(G):
    """Colors the network so that neighboring nodes all have distinct colors.

       Returns a dict keyed by color to a set of nodes with that color.
    """
    coloring = dict()  # color => set(node)
    colors = nx.coloring.greedy_color(G)
    for node, color in colors.items():
        if color in coloring:
            coloring[color].add(node)
        else:
            coloring[color] = set([node])
    return coloring


def _labeling_complete(labeling, G):
    """Determines whether or not LPA is done.

       Label propagation is complete when all nodes have a label that is
       in the set of highest frequency labels amongst its neighbors.

       Nodes with no neighbors are considered complete.
    """
    return all(labeling[v] in _most_frequent_labels(v, labeling, G)
               for v in G if len(G[v]) > 0)


def _most_frequent_labels(node, labeling, G):
    """Returns a set of all labels with maximum frequency in `labeling`.

       Input `labeling` should be a dict keyed by node to labels.
    """
    if not G[node]:
        # Nodes with no neighbors are themselves a community and are labeled
        # accordingly, hence the immediate if statement.
        return {labeling[node]}

    # Compute the frequencies of all neighbours of node
    freqs = Counter(labeling[q] for q in G[node])
    max_freq = max(freqs.values())
    return {label for label, freq in freqs.items() if freq == max_freq}


def _update_label(node, labeling, G):
    """Updates the label of a node using the Prec-Max tie breaking algorithm

       The algorithm is explained in: 'Community Detection via Semi-Synchronous
       Label Propagation Algorithms' Cordasco and Gargano, 2011
    """
    high_labels = _most_frequent_labels(node, labeling, G)
    if len(high_labels) == 1:
        labeling[node] = high_labels.pop()
    elif len(high_labels) > 1:
        # Prec-Max
        if labeling[node] not in high_labels:
            labeling[node] = max(high_labels)


def detect_comm(g):
    """Unused, but might be worth investigating in future"""
    black_edgeset = set([(u, v) for u, v, d in g.edges(data=True) if d['color'] == 'black'])
    results = []
    # For each connected sub graph in the original graph
    for sub in sorted(nx.connected_component_subgraphs(g), key=lambda x: len(x.nodes()), reverse=True):
        if len(sub.nodes()) > break_n_thresh:
            # Find the minimum spanning forest between grey edges (remove unnecessary edges)
            # A forest is a union of the sub graphs, if they are unconnected
            ori = sub.copy()
            ori.remove_edges_from(black_edgeset)

            spanning_forest = nx.minimum_spanning_tree(ori)

            data = spanning_forest.copy()
            data.add_edges_from(n for n in sub if n in black_edgeset)

            partitions = list(label_propagation_communities(data))
            Q = calc_modularity(data, partitions)

            G2 = random_black_edges_graph(sub.nodes())
            rand = spanning_forest.copy()
            rand.add_edges_from(G2.edges())
            partitions = list(label_propagation_communities(rand))
            Q_rand = calc_modularity(rand, partitions)

            print("Q: ", Q, " Qrand: ",  Q_rand)
            results.append({"Q": Q, "Q_rand": Q_rand})

    return results


def proc_blockmodel(G, samp, break_n_thresh, outdir):
    """
    Go through all the chains in each sample
    :param g:
    :param samp:
    :return:
    """
    subs = [G.subgraph(c) for c in nx.connected_components(G)]
    subs = sorted(subs, key=lambda x: len(x.nodes()), reverse=True)
    chain_results = []
    for index, item in enumerate(subs):  # Process for each chain seperately
        if len(item.nodes()) < break_n_thresh:
            continue
        chain_results += process_chain(item, samp, index, outdir)

    pprint.pprint(chain_results)

    return chain_results


def nearest_neighour_graph(f):
    """Reformat the graph, so that grey edges are only allowed between nearest neighbours on the reference.
    By extention, both nodes of an SV link to any nearest neighbours at both ends; this is a result of an SV
    effectively representing 0 distance in bp between two bits of the reference. SO for two SVs: S1a S1b and S2a and S2b
    are in proximity at S1b and S2a, edges include (S1a, S2a), (S1b, S2a), (S1b, S2b), (S2b, S1a)"""


def analyse_graphs(graphs, break_n_thresh, outdir):
    t = []  #["DB58"]
    results = []
    for i in graphs:
        samp = i.split("/")[-1].split("_")[0]
        if samp in t:
            continue
        print(i)
        G = nx.read_graphml(i)
        #if len(G.nodes()) < break_n_thresh:
        #    continue
        if len([i for i in G.edges(data=True) if i[2]['color'] == 'lightgrey']) == 0:
            continue
        r = proc_blockmodel(G, samp, break_n_thresh, outdir)
        results += r
    df = pd.DataFrame.from_records(results)
    print(df)
    df.to_csv(os.path.join(outdir, "graph_attributes.csv"))


def plot_attributes(df, outdir):
    pvalues = {}
    for item in ("Assortativity", "Q", "Qds", "Diameter", ): #"Max_BC", "Max_CF", "Clustering", "Trans",
                 #"B_Diameter", "B_Max_BC", "B_Max_CF", "B_Clustering", "B_Trans"):
        fig, ax = plt.subplots(layout="constrained")

        #plt.subplots_adjust(bottom=0.2, left=0.45)
        data = df[["{}_data".format(item), "{}_rand".format(item)]]
        p = sns.boxplot(data=data, palette=["dodgerblue", "lightgrey"])
        sns.stripplot(data=data)
        ax.set_xticklabels(["Data", "Random"])  # [item, "{} rand".format(item)])
        plt.ylabel(item)
        #plt.xticks(rotation=45)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        import matplotlib.patches as mpatches

        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.1,
                         box.width, box.height * 0.9])

        # red_patch = mpatches.Patch(color='dodgerblue', label='Data')
        # grey_patch = mpatches.Patch(color='lightgrey', label='Random')
        # plt.legend(handles=[red_patch, grey_patch], loc='upper center', frameon=False, ncol=1, bbox_to_anchor=(0.5, -0.05))

        plt.tight_layout()
        plt.savefig(f"{outdir}/general_graphs/{item}_all.pdf", transparent=True)
        plt.close()
        # plt.show()

    for item in ("Assortativity", "Q", "Qds", "Diameter"): #"Max_BC", "Max_CF", "Clustering", "Trans",
                 #"B_Diameter", "B_Max_BC", "B_Max_CF", "B_Clustering", "B_Trans"):
        print(item)
        pval = stats.ttest_rel(df[item+"_data"], df[item+"_rand"])
        pvalues[item] = pval
        print("Mean {} data, median:".format(item), df[item + "_data"].mean(), np.median(df[item + "_data"]))
        print("Mean {} random:".format(item), df[item + "_rand"].mean(), np.median(df[item + "_rand"]))
        print("P-value", pval)
        print("")
    return pvalues 


### side bars ###
def plot_dfr(dfr, idf, f, thresh):
    chained_below = dfr[(dfr[" "] == "Chained") & (dfr[f"TL < {thresh} kb"])]["SVs"]
    chained_above = dfr[(dfr[" "] == "Chained") & (~dfr[f"TL < {thresh} kb"])]["SVs"]

    non_chained_below = dfr[(dfr[" "] != "Chained") & (dfr[f"TL < {thresh} kb"])]["SVs"]
    non_chained_above = dfr[(dfr[" "] != "Chained") & (~dfr[f"TL < {thresh} kb"])]["SVs"]

    fig, ax = plt.subplots(figsize=(4, 4), layout="constrained")
    #plt.subplots_adjust(left=0.25, right=0.75, top=0.75)
    ax = sns.boxplot(x=" ", hue=f"TL < {thresh} kb", y="SVs", data=dfr)
    ax.set_xticks([0,1], ["Chained", "Non-chained"])
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.75, box.height])

    try:
        a = stats.mannwhitneyu(chained_above, chained_below)
        b = stats.mannwhitneyu(non_chained_above, non_chained_below)
        print("Chained above vs chained below", a)
        print("Non chained above vs non-chained below", b)
        #plt.title("pval_chained={}, \npval_non-chained={}".format(round(a.pvalue, 5), round(b.pvalue, 5)))
        ax.text(0, max([chained_above.max(), chained_below.max()])+1, str(round(a.pvalue, 5)), ha="center", va="bottom")
        ax.text(1, max([non_chained_above.max(), non_chained_below.max()])+1, str(round(b.pvalue, 5)), ha="center", va="bottom")
        chained_stat = a
    except:
        chained_stat = 1
        pass

    # Put a legend to the right of the current axis
    ax.legend(loc='center left', title=f"TL < {thresh} kb", bbox_to_anchor=(1, 0.5), frameon=False)
    #plt.tight_layout()
    plt.savefig("{}/chained_vs_threshold_counts.{}.pdf".format(f, idf), transparent=True)
    return chained_stat


### knee/elbow plot ###
def plot_elbow(x, y, f):
    print("elbow plot: ")
    print("x = ", x)
    print("y = ", y)
    kn = KneeLocator(x, y, curve="convex", direction="decreasing")
    fig, ax = plt.subplots()
    ax.plot(x, y, zorder=10)
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax)
    ax.vlines(kn.knee, ymin, ymax, color="grey", linestyle="--", zorder=1, alpha=0.75)
    if kn.knee != None:
        ax.text(kn.knee, max(y), f"x = {round(kn.knee, 5)}\ny = {round(kn.knee_y, 5)}", va="top")
    ax.set_xlabel("Chained P-val")
    ax.set_ylabel("MannWhitneyU")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig("{}/elbow.pdf".format(f), transparent=True)
    

### stacked bar ###
def plot_chained_vs_non_chained(f, tel, sample_col, stela_col, thresh, stack_unit):
    # Collect the total number of SVs
    df = pd.read_csv(os.path.join(f, "all_svs.unique.chains.csv"), sep=",", index_col=None)
    # print(df)
    tel = pd.read_csv(tel)
    tel = tel.rename({sample_col: "Sample", stela_col: "tel_length"}, axis=1)
    # print(tel)
    df["tel_length"] = 0
    df["ID"] = 0
    
    tel.index = tel["Sample"]
    tel = tel["tel_length"].to_dict()

    df["tel_length"] = df["Sample"].apply(lambda x: tel.get(x))
    allowed = set(["chr" + str(i) for i in range(1, 23)] + ["chrX", "chrY"] + [str(i) for i in range(1, 23)] + ["X", "Y"])
    df = df[df["chrom1"].isin(allowed) & df["chrom2"].isin(allowed)]

    x = df.sort_values(["tel_length", "chained"])[["tel_length", "Sample", "chained", "ID"]]
    sample_list = x["Sample"].tolist()
    length_list = x["tel_length"].tolist()
    #x["Sample"] = [i[2:] for i in x["Sample"]]
    x = x.groupby(["tel_length", "Sample", "chained"]).agg("count")
    x = x.unstack(level=2)
    x.columns = x.columns.droplevel(level=0)
    x = x.fillna(0)

    ax = x.plot(kind="bar", stacked=True, colormap=ListedColormap(sns.color_palette("PuBu", 4)),
           figsize=(12, 6), edgecolor="grey", lw=0.25)
    #plt.subplots_adjust(bottom=0.25)
    plt.tight_layout()
    #locs, labels = plt.xticks()
    #labels = [i.get_text() for i in labels]
    #labels = ["{} - {:.2f}".format(*list(eval(i))[::-1]) for i in labels]
    #labels = ["{} - {:.2f}".format(i[1], i[0]) for i in labels]
    labels = []
    for s, l in zip(sample_list, length_list):
        labels.append(f"{s[:5]} - {'%.2f' % round(float(l), 2)}")
    labels = list(dict.fromkeys(labels))

    plt.xticks(ticks = [i for i in range(len(labels))], labels = labels)
    plt.xlabel(f"Sample - {stack_unit}")
    plt.ylabel("SV count")
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.savefig("{}/chained_vs_non-chained.pdf".format(f), transparent=True)
    #plt.show()
    plt.close()

    # Make a plot for below vs above threshold 
    df["below_thresh"] = df["tel_length"] < thresh
    #df["sample_id"] = [int(i[2:]) for i in df["Sample"]]
    #df["c"] = [1]*len(df)
    x = df.sort_values(["below_thresh"])[["below_thresh", "chained", "ID", "Sample"]]

    recs = []
    for samp, dfr in df.groupby("Sample"):
        # dfr["chained"] = [True if "True" in i else False for i in dfr["chained"]]
        recs.append({" ": "Chained", f"TL < {thresh} kb": dfr["below_thresh"].iloc[0], "SVs": dfr["chained"].sum(), "Sample": samp})
        recs.append({" ": "Non-chained", f"TL < {thresh} kb": dfr["below_thresh"].iloc[0], "SVs": (~dfr["chained"]).sum(), "Sample": samp})

    dfr = pd.DataFrame.from_records(recs)
    print(dfr)

    chained_stat = plot_dfr(dfr, "all", f, thresh)
    
    no_chained = [i for i, j, k in zip(dfr["Sample"], dfr[" "], dfr["SVs"]) if j == "Chained" and k == 0] # + ["DB163"]
    print(no_chained)
    dfr2 = dfr[~dfr["Sample"].isin(no_chained)]
    print(set(dfr2["Sample"]))

    plot_dfr(dfr2, "with_chained", f, thresh)

    # dfr3 = dfr[dfr["Sample"].isin(no_chained)]
    # print(set(dfr3["Sample"]))

    # plot_dfr(dfr3, "without_chained", f)
    return chained_stat



@pipe.command(name="plot")
@click.pass_context
def plot_results(ctx):
    break_n_thresh = ctx.obj["bt"]
    ## if plot being called from command line
    if re.match(r".*found_chains_p\d+\.\d+_t\d+", ctx.obj["indir"]):
        in_dirs = [ctx.obj["indir"]]
    elif re.match(r".*t\d+_p\d+\.\d+-\d+\.\d+", ctx.obj["indir"]):
        in_dirs = os.listdir(ctx.obj["indir"])
        in_dirs = [os.path.join(ctx.obj["indir"], i) for i in in_dirs if os.path.isdir(os.path.join(ctx.obj["indir"], i))]
    ## if plot being called in pipeline
    elif 0 < ctx.obj["pval"] < 1:
        in_dirs = [os.path.join(ctx.obj["outdir"], f"found_chains_p{ctx.obj['pval']}_t{ctx.obj['bt']}")]
    elif ctx.obj["pval"] >= 1:
        print(os.path.join(ctx.obj["outdir"], f"t{ctx.obj['bt']}_p{ctx.obj['lp']}-{ctx.obj['up']}"))
        in_dirs = os.listdir(os.path.join(ctx.obj["outdir"], f"t{ctx.obj['bt']}_p{ctx.obj['lp']}-{ctx.obj['up']}"))
        print(in_dirs)
        in_dirs = [os.path.join(ctx.obj["outdir"], f"t{ctx.obj['bt']}_p{ctx.obj['lp']}-{ctx.obj['up']}", i) for i in in_dirs if os.path.isdir(os.path.join(ctx.obj["outdir"], f"t{ctx.obj['bt']}_p{ctx.obj['lp']}-{ctx.obj['up']}", i))]

    print("in_dirs: ", in_dirs)

    for in_folder in in_dirs:
        if os.path.exists(os.path.join(in_folder, "general_graphs")) == False:
            os.makedirs(os.path.join(in_folder, "general_graphs")) 
        
        break_n_thresh = int(in_folder.split("_")[-1].replace("t", ""))

        graphs = glob.glob(in_folder + "/*/*_graph.graphml")
        print(graphs)

        np.random.seed(123)
        random.seed(123)

        analyse_graphs(graphs, break_n_thresh, in_folder)
        df = pd.read_csv(os.path.join(in_folder, "graph_attributes.csv"))
        pvalues = plot_attributes(df, in_folder)
        
        plt.figure(figsize=(6, 4))
        ax = plt.subplot(111)
        plt.subplots_adjust(bottom=0.2, left=0.35)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        bins = np.linspace(-1, 1, 20)
        pval = pvalues["Assortativity"]
        plt.title("Assortivity coeff pval={:0.3g}".format(pval.pvalue))
        # plt.hist(assortivity, color="red", label="Data", alpha=0.6, bins=bins)
        # plt.hist(assortivity_rand, color="blue", label="Random", alpha=0.6, bins=bins)

        sns.histplot(df["Assortativity_data"], color="red", bins=bins, kde=True, kde_kws={"bw_adjust": 0.15})
        sns.histplot(df["Assortativity_rand"], color="blue", bins=bins, kde=True, kde_kws={"bw_adjust": 0.15})

        plt.xlabel("Assortivity coeff")
        plt.ylabel("Number of chains")

        # Shrink current axis by 20%
        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # Put a legend to the right of the current axis
        # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(f"{in_folder}/general_graphs/assortivity_coeff.pdf")
        # plt.show()

        # df = pd.DataFrame.from_records(all_modularity)
        # df.to_csv(in_folder + "modularity.csv", index=None, sep="\t")

        plt.figure(figsize=(6, 4))
        ax = plt.subplot(111)
        plt.subplots_adjust(bottom=0.2, left=0.35)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        bins = np.linspace(-0.5, 1, 25)
        tt = stats.ttest_rel(df["Q_data"], df["Q_rand"]).pvalue

        plt.title("Modularity pval={:0.3g}".format(tt))
        print("Modularity p-value", tt)
        sns.histplot(df["Q_data"], color="purple", bins=bins, kde=True, kde_kws={"bw_adjust": 0.15})
        sns.histplot(df["Q_rand"], color="orange", bins=bins, kde=True, kde_kws={"bw_adjust": 0.15})

        plt.xlabel("Modularity")
        plt.ylabel("Normed chains")

        # Shrink current axis by 20%
        # box = ax.get_position()
        # ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

        # Put a legend to the right of the current axis
        # ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        plt.savefig(f"{in_folder}/general_graphs/modularity_coeff.pdf")


if __name__ == "__main__":
    pipe()
