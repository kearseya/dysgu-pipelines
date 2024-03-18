import pandas as pd
import numpy as np
from scipy.signal import butter, filtfilt, find_peaks
from scipy.signal import find_peaks_cwt
from scipy.fft import fft, fftfreq
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib as matplotlib
import plotly.graph_objects as go
from statsmodels.robust import mad
import pywt
import seaborn as sns
from scipy.stats import mannwhitneyu
import sounddevice as sd
import itertools

df = pd.read_csv("cnpinter/hawk_1k/copyNumber_io/winsorized.all.csv", low_memory=False)
#print(df)
#c = df[df["chrom"] == 1]
#c = df[df["chrom"] == "1"]
##print(c)
#c = c[["pos", "DB159"]]
##print(c)
#x = np.array(c["DB159"].tolist())
##print(x)
#data = x

# Filter requirements.
T = 5.0         # Sample Period
fs = 30.0       # sample rate, Hz
cutoff = 0.5      # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2Hz
nyq = 0.5 * fs  # Nyquist Frequency
order = 2       # sin wave can be approx represented as quadratic
n = int(T * fs) # total number of samples


T = 5.0         # Sample Period
fs = 30.0       # sample rate, Hz
cutoff = 0.5      # desired cutoff frequency of the filter, Hz ,      slightly higher than actual 1.2Hz
nyq = 0.5 * fs  # Nyquist Frequency
order = 2       # sin wave can be approx represented as quadratic
n = int(T * fs) # total number of samples



def butter_lowpass_filter(data, cutoff, fs, order):
    normal_cutoff = cutoff / nyq
    # Get the filter coefficients
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data, axis=0)
    return y


def color_map_color(value, cmap_name='RdBu', vmin=-2, vmax=2):
    # norm = plt.Normalize(vmin, vmax)
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    cmap = cm.get_cmap(cmap_name)  # PiYG
    rgb = cmap(norm(abs(value)))[:3]  # will return rgba, we take only first 3 so we get rgb
    color = matplotlib.colors.rgb2hex(rgb)
    return color

## test data
# t = np.arange(0, 10, 0.1)
# # sin wave
# sig = np.sin(1.2*2*np.pi*t)
# print(sig)
# # Lets add some noise
# noise = 1.5*np.cos(9*2*np.pi*t) + 0.5*np.sin(12.0*2*np.pi*t)
# print(noise)

# data = sig + noise
# print(data)

# plt.plot(np.arange(0, 10, 0.1), data)
# plt.show()
##
# plt.plot(np.arange(0, len(data), 1), data)
# plt.show()

# Filter the data, and plot both the original and filtered signals.
def plot_data(data):
    l = len(data)
    ## low_pass
    y = butter_lowpass_filter(data, cutoff, fs, order)
    #sd.play(y)
    ## wavelet
    # d = data
    # # # y = pywt.wavedec(d, 'haar', level=2, mode='per')
    # noisy_coefs = pywt.wavedec(d, 'haar', level=2, mode='per')
    # sigma = mad(noisy_coefs[-1])
    # uthresh = sigma * np.sqrt(2 * np.log(len(d)))

    # denoised = noisy_coefs[:]
    # denoised[1:] = (pywt.threshold(i, value=uthresh) for i in denoised[1:])
    # y = pywt.waverec(denoised, 'haar', mode='per')
    #y = pywt.dwt(d, 'db2', 'smooth')
    # y = y[0]
    # l = len(y)

    p, properties = find_peaks(y)
    t, tproperties = find_peaks(-y)
    #print(peaks)
    #print(properties)

    d = []
    di = []
    #pos = []
    sec = []

    # for pe, tr in zip(p, t):
    #     d.append(abs(y[pe] - y[tr]))
        #pos.append(min([pe, tr]))
    #print(d)
    c = []

    if len(p) > len(t):
        c = [i for x in zip(p, t) for i in x] + [p[-1]]
        #sec = [[i, x] for i, x in zip(p, t)]
    else:
        c = [i for x in zip(t, p) for i in x] + [t[-1]]
        #sec = [[i, x] for i, x in zip(t, p)]  
    # print(p)
    # print(t)
    # print(c)
  
    c = list(sorted(c))
    for idx, i in enumerate(c[:-1]):
        d.append(abs(y[c[idx]] - y[c[idx+1]]))
    print(d)

    ## single
    plt.scatter(np.arange(0, l, 1), data, alpha=0.6, color="grey", edgecolor="none", zorder = 5)
    plt.plot(np.arange(0, l, 1), y, zorder = 10)
    plt.plot(p, y[p], "x", zorder=15)
    plt.plot(t, y[t], "x", zorder=15)
    #plt.plot(c, y[c], ".", color="red", zorder=100)

    print(len(d), len(c))
    for idx, diff in enumerate(d):
        if diff > 0.5:
            plt.broken_barh([[c[idx], (c[idx+1]-c[idx])]], (min([y[c[idx]], y[c[idx+1]]]), diff), facecolors="red", alpha=0.3, zorder=2)
            #plt.Rectangle((c[idx], -2), int(c[idx+1]-c[idx]), 7, fc="red", alpha=0.5, zorder=2)
    # plt.vlines(x=peaks, ymin=x[peaks] - properties["prominences"],
    #            ymax = x[peaks], color = "C1")
    # plt.hlines(y=properties["width_heights"], xmin=properties["left_ips"],
    #            xmax=properties["right_ips"], color = "C1")

    plt.xlim(0, l)
    plt.ylim(-2, 5)
    plt.tight_layout()
    plt.title(f"{sample} chr{chrom} score: {round(sum([i for i in d if i > 0.25]), 2)}")
    plt.show()

    ## iter
    #ax[idx].scatter(np.arange(0, l, 1), data, alpha=0.6, color="grey", edgecolor="none", zorder = 5)
    #ax[idx].plot(np.arange(0, l, 1), y, zorder = 10)
    #ax[idx].plot(p, y[p], "x", zorder=15)
    #ax[idx].plot(t, y[t], "x", zorder=15)

    #for s, diff in zip(sec, d):
    #    ax[idx].broken_barh([s], (-2, 5), facecolors=color_map_color(diff), alpha=0.3, zorder=2)
    ## ax[idx].vlines(x=peaks, ymin=x[peaks] - properties["prominences"],
    ##            ymax = x[peaks], color = "C1")
    ## ax[idx].hlines(y=properties["width_heights"], xmin=properties["left_ips"],
    ##            xmax=properties["right_ips"], color = "C1")

    #ax[idx].set_xlim(0, l)
    #ax[idx].set_ylim(-2, 5)
    ##ax[idx].tight_layout()
    #ax[idx].set_title(f"{sample} chr{chrom} score: {round(sum([i for i in d if i > 0.25]), 2)}")
    ##ax[idx].show()


# idx = -1

# fig, ax = plt.subplots(nrows=4, ncols=3)
# ax = ax.ravel()

# for chrom in ["1", "2", "3", "4"]:
#     for sample in ["DB146", "DB159", "DB219"]:
#         idx += 1
#         c = df[df["chrom"] == chrom]
#         c = c[["pos", sample]]
#         x = np.array(c[sample].tolist())
#         data = x
#         plot_data(data)

# plt.tight_layout()
# plt.show()

#idx=0
#for chrom in [str(i) for i in range(1, 23)]:
#    for sample in list(df.columns)[3:]:
#        #fig, ax = plt.subplots(nrows=2, ncols=1)
#        #ax = ax.ravel()
#        data = np.array(df[df["chrom"] == chrom][sample].tolist())
#        plot_data(data)
        #plt.tight_layout()
        #plt.show()

def plot_raw_genome(sample):
    data = df[["chrom", sample]]
    chrom_lines = []
    for chrom in [str(i) for i in range(1, 24)] + ["x"]:
        chrom_lines.append(len(data[data["chrom"] == chrom]))
    chrom_lines = list(np.cumsum(chrom_lines))
    data = np.array(data[sample].tolist())
    sd.play(data)
    l = len(data)
    y = butter_lowpass_filter(data, cutoff, fs, order)
    plt.plot(np.arange(0,l,1), y, color="none")
    plt.fill_between(np.arange(0,l,1), y, 0, where=(y>0), color="red", alpha=0.5)
    plt.fill_between(np.arange(0,l,1), y, 0, where=(y<0), color="blue", alpha=0.5)
    plt.vlines(chrom_lines, -2, 15)
    plt.xlim(0, l)
    plt.show()
    # from scipy.io.wavfile import write
    # write("DB165.wav", 44100, y)

def plot_genome(sample):
    data = df[["chrom", sample]]
    chrom_lines = list(range(1,23))
    #data = np.array(data[sample].tolist())
    l = len(data)
    #y = butter_lowpass_filter(data, cutoff, fs, order)
    chroms = [str(i) for i in range(1, 23)] + ["X"]
    for idx, chrom in enumerate(chroms):
        #print(chrom)
        filt_data = np.array(data[data["chrom"] == chrom][sample].tolist())
        ## low pass
        #y = butter_lowpass_filter(filt_data, cutoff, fs, order)
        ## raw
        #y = filt_data
        d = filt_data
        noisy_coefs = pywt.wavedec(d, 'haar', level=2, mode='per')
        sigma = mad(noisy_coefs[-1])
        uthresh = sigma * np.sqrt(2 * np.log(len(d)))

        denoised = noisy_coefs[:]
        denoised[1:] = (pywt.threshold(i, value=uthresh) for i in denoised[1:])
        y = pywt.waverec(denoised, 'haar', mode='per')
        

        p, properties = find_peaks(y)
        t, tproperties = find_peaks(-y)
        # d = []
        # for pe, tr in zip(p, t):
        #     d.append(abs(y[pe] - y[tr]))

        c = []
        d = []
        if len(p) > len(t):
            c = [i for x in zip(p, t) for i in x] + [p[-1]]
        else:
            c = [i for x in zip(t, p) for i in x] + [t[-1]]
      
        c = list(sorted(c))
        for idx, i in enumerate(c[:-1]):
            d.append(abs(y[c[idx]] - y[c[idx+1]]))


        l = len(y)
        plt.plot(np.linspace(idx,idx+1,l), y, color="none")
        plt.fill_between(np.linspace(idx,idx+1,l), y, 0, where=(y>0), color="red", alpha=0.5)
        plt.fill_between(np.linspace(idx,idx+1,l), y, 0, where=(y<0), color="blue", alpha=0.5)
        plt.text(idx+0.5, 10, f"{int(round(sum([i for i in d if i > 0.5]), 0))}", ha="center")
    plt.vlines(chrom_lines, -2, 15, color="black", linestyle="--")
    plt.xlim(0, 23)
    plt.ylim(-2, 12)
    plt.xticks([i+0.5 for i in range(23)], chroms)
    plt.xlabel("chromosome")
    plt.ylabel("relative copy number")
    plt.title(sample)
    plt.tight_layout()
    plt.show()

# plot_raw_genome("DB219")

print(df)

def plot_stack_genome():#
    ss = 4
    fig = plt.figure(constrained_layout=True)
    ax = fig.subplot_mosaic([["ridge", "ridge", "score", "rp"], ["chrsplit", "gensplit", "dummy", "dummy"]], gridspec_kw={"width_ratios": [10, 10, 2, 1], "height_ratios": [10, 1]})
    samples = list(df.columns)[3:]
    scores = {}
    for ypos, sample in enumerate(samples):
        scores[sample] = {}
        data = df[["chrom", sample]]
        chrom_lines = list(range(1,23))
        #data = np.array(data[sample].tolist())
        l = len(data)
        #y = butter_lowpass_filter(data, cutoff, fs, order)
        chroms = [str(i) for i in range(1, 23)] + ["X"]
        for idx, chrom in enumerate(chroms):
            #print(chrom)
            filt_data = np.array(data[data["chrom"] == chrom][sample].tolist())
            ## low pass
            y = butter_lowpass_filter(filt_data, cutoff, fs, order)
            p, properties = find_peaks(y)
            t, tproperties = find_peaks(-y)
            # d = []
            # for pe, tr in zip(p, t):
            #     d.append(abs(y[pe] - y[tr]))

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
            scores[sample][chrom] = int(round(sum([i for i in d if i > 0.25]), 0))
            if int(round(sum([i for i in d if i > 0.5]), 0)) > 3:
                ax["ridge"].text(idx+0.5, (ypos*ss)-(ss/2), f"{int(round(sum([i for i in d if i > 0.5]), 0))}", ha="center")
    ax["ridge"].vlines(chrom_lines, -2, (ypos*ss)+5, color="black", linestyle="--", linewidth=0.5)
    ax["ridge"].set_xlim(0, 23)
    ax["ridge"].set_ylim(-2, (ypos*ss)+ss)
    ax["ridge"].set_xticks([i+0.5 for i in range(23)], chroms)
    ax["ridge"].set_yticks([i*ss for i in range(len(samples))], [s[:5] for s in samples])
    ax["ridge"].set_xlabel("chromosome")
    ax["ridge"].set_ylabel("relative copy number")
    ax["ridge"].set_title("Low Pass filtered coverage")
    ax["ridge"].spines[["right", "top", "left", "bottom"]].set_visible(False)
    ax["ridge"].tick_params(length=0)
    # plt.tight_layout()
    # plt.show()

    scores = pd.DataFrame.from_records(scores)
    scores = scores[samples]
    swarm_df = pd.DataFrame({"scores": [], "sample": []})
    total_scores = []
    for ypos, sample in enumerate(samples):
        swarm_df = pd.concat([swarm_df, pd.DataFrame({"scores": scores[sample].tolist(), "sample": [sample]*23})], ignore_index=True)
        total_scores.append(np.sum(scores[sample]))
    print(total_scores)
    recursive_partition = []
    for p in range(1, len(total_scores)):
        recursive_partition.append(mannwhitneyu(total_scores[:p], total_scores[p:], alternative="greater", axis=0).pvalue)
    print(recursive_partition)
    #swarm_df = pd.concat([swarm_df, pd.DataFrame({"scores": [-10], "sample": ["spacer"]})], ignore_index=True)
    print(scores)
    print(swarm_df)
    sns.boxplot(ax=ax["score"], data=swarm_df, x="scores", y="sample",
            meanprops={'marker' : 'D', 'markeredgecolor' : 'black', 'markersize' : 6, 'label':'mean'},
            #medianprops={'visible': False}, whiskerprops={'visible': False},
            showmeans=True, showfliers=False, showbox=True, showcaps=False, boxprops=dict(alpha=.5), zorder=10)
    sns.stripplot(ax=ax["score"], data=swarm_df, x="scores", y="sample", zorder=0)
    ax["score"].set_xlim(-0.5, max(swarm_df["scores"])+1)
    ax["score"].set_ylim(-2/ss, len(samples))
    #ax["score"].axis("off")
    ax["score"].tick_params(left=False, labelleft=False)
    ax["score"].set_xscale("symlog", base=2)
    #ax["score"].set_yticks([i*ss for i in range(len(samples))], samples)
    ax["score"].set(ylabel=None)
    ax["score"].grid(which="major", axis="x")
    #ax["score"].invert_yaxis()

    ax["rp"].plot(recursive_partition, [i*ss for i in range(len(samples)-1)])
    ax["rp"].set_ylim(-2, (len(samples)*ss))
    #ax["rp"].set_xlim(0, 1)
    ax["rp"].set_xticks([0.005,0.010,0.050,0.100,0.500,1],[0.005,0.010,0.050,0.100,0.500,1])
    ax["rp"].set_xscale("log")
    ax["rp"].set_xlabel("pvalue")
    
    ax["rp"].vlines(0.001, -2, (len(samples)*ss)+ss, linestyle="--", linewidth=0.5)
    ax["rp"].vlines(0.01, -2, (len(samples)*ss)+ss, linestyle="--", linewidth=0.5)
    ax["rp"].vlines(0.05, -2, (len(samples)*ss)+ss, linestyle="--", linewidth=0.5)

    #ax["rp"].hlines(np.argmin(recursive_partition)*ss, 0, 1, linestyle="--", linewidth=0.5, color="red")
    rps, _ = find_peaks(-np.array(recursive_partition))
    rps = [i for i in rps if recursive_partition[i] <= min([0.001, 0.005, 0.01, 0.05], key=lambda x:abs(x-min(recursive_partition)))]
    ax["rp"].hlines([i*ss for i in rps], 0, 1, linestyle="--", linewidth=0.5, color="red")

    ax["rp"].set(ylabel=None)
    ax["rp"].tick_params(left=False, labelleft=False)

    split = np.argmin(recursive_partition)
    #chrsplit = {"short": np.array(list(itertools.chain(*scores[scores.columns[:split]].values.tolist()))), "long": np.array(list(itertools.chain(*scores[scores.columns[split:]].values.tolist())))}
    #gensplit = {"short": total_scores[:split], "long": total_scores[split:]}
    #chrsplit = pd.DataFrame.from_dict(chrsplit, orient="index").transpose()
    #gensplit = pd.DataFrame.from_dict(gensplit, orient="index").transpose()
    chrsplit = pd.DataFrame(columns=["scores", "group"])
    chrsplit = pd.concat([chrsplit, pd.DataFrame.from_dict({"scores": np.array(list(itertools.chain(*scores[scores.columns[:split]].values.tolist()))), "group": "short"})])
    chrsplit = pd.concat([chrsplit, pd.DataFrame.from_dict({"scores": np.array(list(itertools.chain(*scores[scores.columns[split:]].values.tolist()))), "group": "long"})])

    gensplit = pd.DataFrame(columns=["scores", "group"])
    gensplit = pd.concat([gensplit, pd.DataFrame.from_dict({"scores": total_scores[:split], "group": "short"})])
    gensplit = pd.concat([gensplit, pd.DataFrame.from_dict({"scores": total_scores[split:], "group": "long"})])

    print(chrsplit)
    print(gensplit)

    sns.boxplot(ax=ax["chrsplit"], data=chrsplit, x="scores", y="group",
            meanprops={'marker' : 'D', 'markeredgecolor' : 'black', 'markersize' : 6, 'label':'mean'},
            #medianprops={'visible': False}, whiskerprops={'visible': False},
            showmeans=True, showfliers=False, showbox=True, showcaps=False, boxprops=dict(alpha=.5), zorder=10)
    sns.stripplot(ax=ax["chrsplit"], data=chrsplit, x="scores", y="group", zorder=0)
    sns.boxplot(ax=ax["gensplit"], data=gensplit, x="scores", y="group",
            meanprops={'marker' : 'D', 'markeredgecolor' : 'black', 'markersize' : 6, 'label':'mean'},
            #medianprops={'visible': False}, whiskerprops={'visible': False},
            showmeans=True, showfliers=False, showbox=True, showcaps=False, boxprops=dict(alpha=.5), zorder=10)
    sns.stripplot(ax=ax["gensplit"], data=gensplit, x="scores", y="group", zorder=0)
    ax["chrsplit"].set_xscale("symlog", base=2)
    ax["gensplit"].set_xscale("symlog", base=2)
    ax["chrsplit"].set_xlim(-0.5, max(chrsplit["scores"])+1)
    ax["gensplit"].set_xlim(min(gensplit["scores"])-0.5, max(gensplit["scores"])+1)
    ax["chrsplit"].set(ylabel=None, xlabel=None, yticks=[])
    ax["gensplit"].set(ylabel=None, xlabel=None, yticks=[])
    chrval = mannwhitneyu(chrsplit[chrsplit['group'] == 'short']['scores'].tolist(), chrsplit[chrsplit['group'] == 'long']['scores'].tolist(), alternative='greater')
    genval = mannwhitneyu(gensplit[gensplit['group'] == 'short']['scores'].tolist(), gensplit[gensplit['group'] == 'long']['scores'].tolist(), alternative='greater')

    ax["chrsplit"].legend(labels=[], labelcolor=["C0", "C1"], loc=2, title=f"chromosome\n$p = {chrval.pvalue:.3g}$")
    ax["gensplit"].legend(labels=[], labelcolor=["C0", "C1"], loc=2, title=f"genome\n$p = {genval.pvalue:.3g}$")
    ax["dummy"].axis("off")

    plt.show()


plot_stack_genome()






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

    f = Rbf(x, y, z, smooth=5)
    z = f(x, y)

    if True:
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





# fig = go.Figure()
# fig.add_trace(go.Scatter(
#             y = data,
#             line =  dict(shape =  'spline' ),
#             name = 'signal with noise'
#             ))
# fig.add_trace(go.Scatter(
#             y = y,
#             line =  dict(shape =  'spline' ),
#             name = 'filtered signal'
#             ))
# fig.show()

