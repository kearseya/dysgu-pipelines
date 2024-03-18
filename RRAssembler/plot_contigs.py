import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import numpy as np
from pyfaidx import Fasta
from dysgu_repeats import compute_rep
import os
from scipy.stats import mannwhitneyu
#from tandem_detection import tandem_finder

def set_colors():
    chromosomes = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    colors = cm.rainbow(np.linspace(0, 1, 23))
    cd = {}
    for ch, cl in zip(chromosomes, colors):
        cd[ch] = cl
    cd["other"] = "grey"
    return cd

def read_prox_bed(fp):
    prox_df = pd.read_csv(fp, sep="\t")
    mm = {}
    for c in prox_df.groupby("#chrom"):
        mm[c[0]] = []
        tup = list(c[1][["chromStart", "chromEnd"]].itertuples(index=False, name=None))
        for p, i in enumerate(tup[:-1]):
            if i[0] < tup[p+1][0]:
                mm[c[0]].append(i[1])
                mm[c[0]].append(tup[p+1][0])
            else:
                mm[c[0]].append(i[0])
                mm[c[0]].append(tup[p+1][1])
        mm[c[0]] = list(set(mm[c[0]]))
    return mm

def plot_individual(c, cd, ref, show_regions=True, show_sequence=True, show_reference=False, show_plot=False, show_prox=True):
    plt.rcParams["figure.figsize"] = (20,3)
    bs = 0
    bw = 5
    basecol = {"T": "green", "C": "red", "A": "blue", "G": "yellow", "N": "grey"}
    revscol = {"A": "green", "G": "red", "T": "blue", "C": "yellow", "N": "grey"}
    if show_prox == True:
        mm = read_prox_bed("hg38_telomeres.bed")
        pcm = cm.get_cmap("hot_r", 10000000)
    if len(c[1]) > 1: 
        # mapq = c[1]["mapq"].tolist()
        # if mapq[0] == 0 or mapq[-1] == 0:
        #     continue
        fig, ax = plt.subplots()
        ax.axis("off")
        startend = list(c[1][["qstart", "qend"]].itertuples(index=False, name=None))
        values = list(c[1][["qstart", "qlength"]].itertuples(index=False, name=None))
        colalpha = list(c[1][["chrom", "mapq"]].itertuples(index=False, name=None))
        xreg = [i[0]+(i[1]/2) for i in values]
        for v, m in zip(values, colalpha):
            #plt.broken_barh([v], (bs, bw), facecolors=[cd[i] for i in c[1]["chrom"].tolist()], ec="black", alpha=m/60)
            if m[1] < 20:
                a = 20
            else:
                a = m[1]
            try:
                ax.broken_barh([v], (bs, bw), facecolors=cd[m[0]], ec="black", alpha=a/70)
            except:
                ax.broken_barh([v], (bs, bw), facecolors="grey", ec="black", alpha=a/70) 
        if len(startend) > 1:
            for pos, i in enumerate(startend[:-1]):
                gap = startend[pos+1][0]-i[1]
                if gap > 0:
                    y = bs+(bw/2)
                    ax.plot([i[1], startend[pos+1][0]], [y, y], color="black")
                    ax.text(x=i[1]+(gap/2), y=y, s=str(gap), ha="center")
        for i, s in enumerate(list(c[1]["strand"])):
            ax.annotate(s, xy=(xreg[i], bs+bw/2), xytext=(xreg[i], bs+bw/2), ha="center", va="center")
        if show_regions == True:
            regions = [f"{c}:{s}-{e}" for c, s, e in zip(c[1]["chrom"], c[1]["rstart"], c[1]["qlength"])]
            yreg = np.empty((len(c[1]),),int)
            yreg[::2] = bs+bw
            yreg[1::2] = bs
            for i, r in enumerate(regions):
                if i % 2 == 0:
                    ax.annotate(r, xy=(xreg[i], yreg[i]), xytext=(xreg[i], yreg[i]), ha="center")
                else:
                    ax.annotate(r, xy=(xreg[i], yreg[i]), xytext=(xreg[i], yreg[i]), ha="center", va="top")
        if show_sequence == True:
            ax.axis("on")
            ## seq
            seq = list(c[1]["seq"].dropna())[0]
            contigseq = seq
            seq = [*seq]
            seqc = [basecol[i.upper()] for i in seq]
            cig = c[1]["cigar"].tolist()[0]
            # print(cig)
            ax.broken_barh([(i, 1) for i in range(len(seq))], (bs, bw/8), facecolors=seqc)
            ## ref
            if show_reference == True:
                regions = [(c, s, e, l) for c, s, e, l in zip(c[1]["chrom"], c[1]["rstart"], c[1]["rend"], c[1]["qlength"])]
                ref = Fasta(ref)
                seqr = []
                for i, r in enumerate(regions):
                    tc, s, e, l = r
                    if tc != "rep":
                        seqt = ref[tc][s-1:e].seq
                        print(compute_rep(seqt))
                        seqr.append([*seqt])
                # print(seqr)
                seqrc = []
                maind = c[1][c[1]["seq"] != np.nan]["strand"].values[0]
                # print(maind)
                for s, d in zip(seqr, list(c[1]["strand"])):
                    # seqrc.append([basecol[i.upper()] for i in s])
                    if d == "+" and maind == "+":
                        seqrc.append([basecol[i.upper()] for i in s[::-1]])
                    elif d == "-" and maind == "+":
                        seqrc.append([revscol[i.upper()] for i in s])
                    elif d == "+" and maind == "-":
                        seqrc.append([basecol[i.upper()] for i in s[::-1]])
                    else:
                        seqrc.append([revscol[i.upper()] for i in s])
                # print(startend)
                # print(seqrc)
                for r, s in zip(startend, seqrc):
                    ax.broken_barh([(i+r[0], 1) for i in range(len(s))], (bs+bw/8, bw/8), facecolors=s)

                #ax.set_xticks(ticks = [i for i in range(len(seq))], labels = [*seq])

            if show_prox == True:
                chrstart = list(c[1][["chrom", "rstart", "qstart", "qlength"]].itertuples(index=False, name=None))
                proxarr = []
                for chrom, r, s, l in chrstart:
                    mindist = 10000000
                    if chrom in mm:
                        for i in mm[chrom]:
                            if mindist > abs(r-i):
                                mindist = abs(r-i)
                                print(mindist)
                    proxarr.append(mindist)
                ax.broken_barh([(s, l) for chrom, r, s, l in chrstart], (bs+(bw/16)*15, bw/16), facecolors=pcm(proxarr)) 


        contig_name = f"{c[1]['sample'].values[0]}_{len(c[1])}_{c[1]['qlen'].values[0]}_{round(compute_rep(contigseq), 3)}"
        for pos in c[1][["chrom", "rstart", "qlength"]].itertuples(index=False, name=None):
            contig_name += f"_{pos[0]}_{pos[1]}_{pos[2]}"
        
        ax.set_title(f"{c[1]['sample'].values[0]} {c[1]['qlen'].values[0]}")
        
        #plt.legend(list(cd.keys()), list(cd.values()))
        plt.tight_layout()
        mapq = c[1]["mapq"].tolist()

        if show_plot == True:
            if mapq[0] == 0 or mapq[-1] == 0:
                print("mapq0")
            if len(c[1][~c[1]["vcf"]]) > 0:
                print("not seen")
            else:
                print("seen")
            plt.show()

        else:
            df_len = len(c[1])
            sample = c[1]["sample"].values[0]
            if mapq[0] == 0 or mapq[-1] == 0:
                if 2 <= df_len <= 3:
                    plt.savefig(f"individual_contigs/map0/{sample}/{df_len}/{contig_name}.png")
                else:
                    plt.savefig(f"individual_contigs/map0/{sample}/{contig_name}.png")
                return
            if len(c[1][~c[1]["vcf"]]) > 0:
                if 2 <= df_len <= 3:
                    plt.savefig(f"individual_contigs/vcf/{sample}/{df_len}/{contig_name}.png")
                else:
                    plt.savefig(f"individual_contigs/vcf/{sample}/{contig_name}.png")
                return
            else:
                if 2 <= df_len <= 3:
                    plt.savefig(f"individual_contigs/not_vcf/{sample}/{df_len}/{contig_name}.png")
                else:
                    #plt.savefig(f"individual_contigs/seen/{c[1]['qname'].values[0]}.png")
                    plt.savefig(f"individual_contigs/not_vcf/{sample}/{contig_name}.png")
                return



def load_data(psl_file, ref):
    df = pd.read_csv(psl_file, sep="\t")
    df["qlength"] = df["qend"] - df["qstart"]
    df[["sample", "1", "2", "3", "4", "5"]] = df['qname'].str.split('_', expand=True)
    #df = df.sort_values(["qlen"], ascending=False)
    print(df)
    cd = set_colors()

    map_path = "individual_contigs/map0"
    vcf_path = "individual_contigs/vcf"
    nvcf_path = "individual_contigs/not_vcf"
    samples = set(df["sample"].tolist())
    for path in [map_path, vcf_path, nvcf_path]:
        if os.path.exists(path) == False:
            os.makedirs(path)
        for sample in samples:
            if os.path.exists(os.path.join(path, sample)) == False:
                os.makedirs(os.path.join(path, sample))
            for size in ["2", "3"]:
                if os.path.exists(os.path.join(path, sample, size)) == False:
                    os.makedirs(os.path.join(path, sample, size))
    
    individual_plot = False
    if individual_plot == True:
        plt.rcParams["figure.figsize"] = (20,3)
        bs = 0
        bw = 5
        show_regions = True
        show_sequence = True
        basecol = {"T": "green", "C": "red", "A": "blue", "G": "yellow", "N": "grey"}
        revscol = {"A": "green", "G": "red", "T": "blue", "C": "yellow", "N": "grey"}
        for c in df.groupby("qname", sort=False):
            if len(c[1]) > 1: 
               plot_individual(c, cd, ref, show_plot=True) 
        individual_plot = False
    else:
        plt.rcParams["figure.figsize"] = (20,20)
        show_regions = False
        show_prox = True
        if show_prox == True:
            mm = read_prox_bed("hg38_telomeres.bed")
            pcm = cm.get_cmap("hot_r", 10000000)
            proxd = {}
        for tdf in df.groupby("sample"):
            bs = 0
            bw = 1
            if show_prox == True:
                proxd[tdf[0]] = []
            ttdf = tdf[1].sort_values(["qlen", "qstart"], ascending=True).reset_index()
            ttdf = ttdf.sort_values(["n_alignments", "qstart"], ascending=True).reset_index()
            print(ttdf)
            print(len(tdf[1]))
            fig, ax = plt.subplots()
            n_segs = []
            for c in ttdf.groupby("qname", sort=False):
                n_segs.append(len(c[1]))
                mapq = c[1]["mapq"].tolist()
                chromosomes_involved = set(c[1]["chrom"].tolist())
                if mapq[0] == 0 or mapq[-1] == 0 or len(chromosomes_involved) == 1:
                    continue
                if len(c[1]) > 2:
                    startend = list(c[1][["qstart", "qend"]].itertuples(index=False, name=None))
                    values = list(c[1][["qstart", "qlength"]].itertuples(index=False, name=None))
                    colalpha = list(c[1][["chrom", "mapq"]].itertuples(index=False, name=None))
                    for v, m in zip(values, colalpha):
                        #plt.broken_barh([v], (bs, bw), facecolors=[cd[i] for i in c[1]["chrom"].tolist()], ec="black", alpha=m/60)
                        if m[1] < 20:
                            a = 20
                        else:
                            a = m[1]
                        try:
                            ax.broken_barh([v], (bs, bw), facecolors=cd[m[0]], ec="black", alpha=a/70)
                        except:
                            ax.broken_barh([v], (bs, bw), facecolors="grey", ec="black", alpha=a/70)
                    if len(startend) > 1:
                        for pos, i in enumerate(startend[:-1]):
                            gap = startend[pos+1][0]-i[1]
                            if gap > 0:
                                y = bs+(bw/2)
                                ax.plot([i[1], startend[pos+1][0]], [y, y], color="black")
                                #plt.text(x=i[1]+(gap/2), y=y, s=str(gap), ha="center")
                    if show_prox == True:
                        chrstart = list(c[1][["chrom", "rstart", "qstart", "qlength"]].itertuples(index=False, name=None))
                        proxarr = []
                        for chrom, r, s, l in chrstart:
                            mindist = 10000000
                            if chrom in mm:
                                for i in mm[chrom]:
                                    if mindist > abs(r-i):
                                        mindist = abs(r-i)
                                        print(mindist)
                            proxarr.append(mindist)
                        proxd[tdf[0]].append(proxarr)
                        ax.broken_barh([(s, l) for chrom, r, s, l in chrstart], (bs+(bw/4)*3, bw/4), facecolors=pcm(proxarr)) 

                    # if show_regions == True:
                    #     regions = [f"{c}:{s}-{e}" for c, s, e in zip(c[1]["chrom"], c[1]["rstart"], c[1]["qlength"])]
                    #     print(regions)
                    #     xreg = [i[0]+(i[1]/2) for i in values]
                    #     print(list(c[1][["qstart", "qend"]].itertuples(index=False, name=None)))
                    #     print(values)
                    #     print(xreg)
                    #     yreg = np.empty((len(c[1]),),int)
                    #     yreg[::2] = bs+bw
                    #     yreg[1::2] = bs
                    #     for i, r in enumerate(regions):
                    #         if i % 2 == 0:
                    #             plt.annotate(r, xy=(xreg[i], yreg[i]), xytext=(xreg[i], yreg[i]), ha="center")
                    #         else:
                    #             plt.annotate(r, xy=(xreg[i], yreg[i]), xytext=(xreg[i], yreg[i]), ha="center", va="top")

                    #plt.title(f"{c[1]['sample'].values[0]} {c[1]['qlen'].values[0]}")
                    #plt.legend(list(cd.keys()), list(cd.values()))
                    #plt.axis("off")
                    #plt.show()
                    bs += 1
                    print(bs, end="\r")

            if show_prox == True:
                #proxd[tdf[0]] = proxarr
                sm = mpl.cm.ScalarMappable(cmap=pcm)
                sm.set_array([0, 5000000, 10000000])
                fig.colorbar(sm, ax=ax)

            #plt.axis("off")
            ax.set_title(tdf[1]["sample"].values[0])
            ax.legend(labels=list(cd.keys()), loc="center right")
            #ax = plt.gca()
            leg = ax.get_legend()
            for x, i in enumerate(list(cd.keys())):
                try:
                    leg.legendHandles[x].set_color(cd[i])
                except:
                    continue
            from mpl_toolkits.axes_grid.inset_locator import inset_axes
            histax = inset_axes(ax, width="30%", height="30%", loc="lower right")
            histax.hist(n_segs)
            plt.tight_layout() 
            #plt.show()
            plt.savefig(f"all_plots/{tdf[1]['sample'].values[0]}_all.png")
    if show_prox == True:
        print(proxd)
        for s in list(proxd.keys()):
            proxd[s] = [v for b in proxd[s] for v in b]
        order = [('DB219', 2.273732984), ('DB203', 2.407897541), ('DB229', 2.79013245), ('DB177', 2.884216374), ('DB195', 2.919859504), ('DB191', 3.176696429), ('DB197', 3.221461538), ('DB199', 3.263198795), ('DB169', 3.342906863), ('DB144', 3.356428), ('DB209', 3.411828767), ('DB223', 3.439524862), ('DB217', 3.518849315), ('DB148', 3.590629412), ('DB159', 3.599047619), ('DB185', 3.735082418), ('DB213', 3.74), ('DB150', 3.894544554), ('DB152', 4.065735294), ('DB183', 4.425712707), ('DB187', 4.464242424), ('DB156', 4.640192308), ('DB165', 4.769766169), ('DB189', 4.77673913), ('DB161', 4.896705426), ('DB154', 4.903419118), ('DB227', 4.942786047), ('DB181', 4.998), ('DB207', 5.022224719), ('DB171', 5.258006897), ('DB167', 5.426142857), ('DB225', 5.654782609), ('DB215', 5.7130625), ('DB173', 5.730572581), ('DB146', 5.817026316), ('DB205', 5.826369369), ('DB175', 5.84973913), ('DB211', 5.965103093), ('DB201', 5.98139726), ('DB179', 6.008243243), ('DB234', 6.10312), ('DB221', 6.3486), ('DB193', 6.535644444), ('DB163', 7.056571429)]
        short_samples = {'DB209', 'DB159', 'DB169', 'DB203', 'DB217', 'DB144', 'DB177', 'DB197', 'DB191', 'DB219', 'DB223', 'DB148', 'DB229', 'DB199', 'DB195'}
        long_samples = {'DB221', 'DB234', 'DB227', 'DB225', 'DB171', 'DB146', 'DB150', 'DB189', 'DB193', 'DB205', 'DB207', 'DB154', 'DB187', 'DB215', 'DB201', 'DB179', 'DB167', 'DB173', 'DB181', 'DB163', 'DB165', 'DB161', 'DB156', 'DB211', 'DB152', 'DB183', 'DB175'}
        for lim in [50000, 100000, 500000, 1000000, 5000000, 10000000]:
            print("limit: ", lim)
            short = []
            long = []
            for i in order:
                if i[0] not in list(proxd.keys()):
                        continue
                lt = 0
                for x in proxd[i[0]]:
                        if x < lim:
                                lt += 1
                if len(proxd[i[0]]) == 0:
                    print(i[0], "empty")
                    if i[0] in short_samples:
                        short.append(0)#lt/len(proxd[i[0]]))
                    if i[0] in long_samples:
                        long.append(0)#lt/len(proxd[i[0]]))

                else:
                    print(i[0], "\t", lt, "\t", len(proxd[i[0]]), "\t", lt/len(proxd[i[0]]))
                    if i[0] in short_samples:
                        short.append(lt/len(proxd[i[0]]))
                    if i[0] in long_samples:
                        long.append(lt/len(proxd[i[0]]))
            print("mannwhitneyu: ", mannwhitneyu(short, long, alternative="greater"))

    return df

#load_data("out_data/map_info.bed")
#load_data("test.out")

