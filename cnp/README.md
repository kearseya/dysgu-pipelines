# Copy Number Plotting

This pipeline takes raw coverage beds of a cohort, and by using a csv containing the tumour and normal names will normalise coverage values to the GC% and mappability of each region. Once normalised the coverage values from the normal will be subtracted from the tumour to create a relative copy number profile.

Whilst the whole process can be run via one command, there are several subcommands that split up each stage of the analysis:

## Overview

This script takes the raw coverage bed files for both the tumour and normal and normalises the values based on the GC percentage and map-ability of the region. A matrix for coverage for each combination of GC% between 30 and 60, and mappability values between 60 and 101 is calculated. The median of all these values is then used as the “genomic median” which is then subtracted from the matrix coverage values. This is then interpolated with RBFInterpolation from scipy before using these normalised values to subtract from the raw coverage values. A haar wavelet transformation is then applied to normalised coverage values. Finally the normal coverage can be subtracted from the tumour coverage to create a relative copy number profile displaying the gains and losses of the tumour sample.

The tumour values which have now been normalised by their normal counterparts are winsorized (similar to clipped) before performing piecewise constant segmentation (pcf). These segments can then be used to calculate a “complexity score” which is the sum of the absolute relative copy number of all segments. This method of scoring gives higher values to samples with more segments, whilst minimising the effect of a copy number gain of a chromosome arm for example. Another more experimental form of analysis is also performed by the script by using a low pass filter over the relative copy number values. In this case the score is calculated by gathering all peaks and troughs of the processed signal and measuring the absolute sum of all cases where the peak and trough are >0.5 apart.


## Usage

```
Usage: copy_number_pipeline.py [OPTIONS] COMMAND [ARGS]...

  Tool related commands

Options:
  -i, --indir TEXT         input dir with raw coverages
  -r, --ref TEXT
  -w, --winsize INTEGER    window size to generate files for
  -b, --bg TEXT...         pairs csv with stela [file, tumour col, tumour
                           stela, normal col]
  -s, --subsample INTEGER  n samples to plot
  -t, --threshold FLOAT    threshold to split groups
  --split FLOAT            percentage split of groups in subsample to plot
  --subsets INTEGER        n subsets of size subsample
  -c, --categorical        Use if split categorical
  -t, --tag TEXT           tag to identify dataset
  -g, --gamma INTEGER      pcf gamma (lower = more relaxed)
  -k, --kmin INTEGER       pcf kmin (smallest segment size = kmin*winsize)
  -p, --plot-inter         plot output from intermediary steps
  --yaxis-scale            replace y-tick text with color scale
  --ignore-score           remove score boxplot from perif
  --ignore-box             dont plot score box in perif
  --skip-genome            skip plotting individual genome coverages
  --help                   Show this message and exit.

Commands:
  compare     Plots normalised read counts (2nd)
  gainloss    Plot gains vs losses (4th)
  mergeout
  normalise   Generates the x.cov.bed files (1st)
  plot        Plot normalised coverage of genomes (3rd)
  preprocess

```

## Preprocessing

This stage is a walk through of how to create the prerequisit files required for running (a gc percentage bed which can be generated from the reference, and a mappability file that will have to be externally sourced). Files for a 10kb window size for hg38 and hg19 are provided, as well as a 1mb window for hg38.

Note to add more references, instances where `elif ref.lower() in {"hg38", "grch38"}:` are seen, a `chrom_lengths_dict` will need to be added to the code (eg. lines ~192, ~338, ~600).

## Normalising

Normalise to GC% and mappability.

Find the genomic median:
```
s = s[(s["GC"].between(30, 60)) & (s["map"].between(70, 101))]
median = np.median(s["coverage"])
```

Normalise to the genomic median:
```
x, y, z = [], [], []
arr = np.zeros((101, 101))
for (i, j), v in norm_array.items():
    if v < 1000:
        sub = v - median
        arr[i, j] = sub
        x.append(j), y.append(i), z.append(sub)

f = RBFInterpolator(x, y, z, smooth=5)
z = f(x, y)
...
norm_value = {(int(xx), int(yy)): zz for xx, yy, zz in zip(x, y, z)}
for g, m, v in zip(list(s["GC"]), list(s["map"]), list(s["coverage"])):
    if (m, g) in norm_value:
        normed_counts.append(v - norm_value[(m, g)])
    else:
        normed_counts.append(0)
```

Wavelet transformation (Haar):
```
d = np.array(normed_counts)

noisy_coefs = pywt.wavedec(d, 'haar', level=2, mode='per')
sigma = mad(noisy_coefs[-1])
uthresh = sigma * np.sqrt(2 * np.log(len(d)))

denoised = noisy_coefs[:]
denoised[1:] = (pywt.threshold(i, value=uthresh) for i in denoised[1:])
sig = pywt.waverec(denoised, 'haar', mode='per')

if len(sig) > len(s["GC"]):
    sig = sig[:-1]
```

## Plot

Plotting can be called after the normalising step has been run. It is recommended to use this command to generate if the pipeline has alreaady been run once and plots require editting. Whilst most plotting options can be changed at run time in the command line, booleans at lines ~1235-1246 can also be set to get the desired outcome.

Relative copy number is visualised in serveral forms, the plot class can plot values for: individual chromosomes and/or whole genome per sample, and for the cohort a heatmap, "squiggle" plot (logged segments over genome with normalisation to chromosome size), low pass filtered (squiggle and low pass w/ or w/o periferal).

"Complexity scores" are caluclated for both the pcf segmentation and low pass. For pcf, this is done by calculating the absolute sum of all segments per sample with:

```
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
```

For low pass, this is done by finding the sum of all differences between the peaks and troughs which are more than 0.5:

```
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
```

## Gains vs Losses

Here graphs relating to the gains vs losses can be generated. Gains and losses are defined as segments that are above and below a threshold (0.05).

```
ignore_thresh = 0.05
data["Gain"] = [1 if i > ignore_thresh else 0 for i in list(data["mean"])]
data["Loss"] = [1 if i < -ignore_thresh else 0 for i in list(data["mean"])]
```

## Merge out

Used for merging the generated figures into one pdf using the pypdf library (useful for when making plots inside a research environment that have an airlock type system).
