# Circos Plotting

Circos plots can be generated from just dysgu vcf files. but this script can also accept optional files such as cytoband files (hg38 and hg19 included) and the output from the copy number script to add more information.

Pycircos is the python library used to generate the circos plots within this script. The length values from the vcf header are extracted using regex to create a dictionary of the lengths which can then be used to create the “Garc” objects which are added to a “Gcircle” object. Optionally at this stage iterating a cytoband file can be used to add banding in the outer chromosome arcs, and reading of a copy number profile from the previous pipeline to add copy number data around the inner edge of each chromosome (highlighting chromosomes with high scores or number of segments). Pysam is then used to iterate over the variants within the vcf to extract the SV type, loci, length data and excluding variants that do not meet the filtering criteria (probability, length, supporting reads). Variants that are not classed as insertions are then plotted as chord plots, and insertions as scatter around the edge. The plots are then saved in SVG format, and are then be stitched together using svgutils transform commands.

## Usage

```
Options:
  -od, --outdir TEXT      Output directory  [default: tmp_svg]
  -oname, --outname TEXT  Dataset name (output file name)  [default: all]
  -n, --ncols INTEGER     Number of columns in combined output plot  [default:
                          8]
  --skip                  Skip making svgs, arrange existing
  -l, --lengths TEXT      Order csv (optional)
  -c, --col TEXT          Column in order dataframe to sort (continous) or
                          split (boolean) by  [default: short]
  --col-thresh TEXT       Threshold to split groups, required if order column
                          is not boolean
  --subsets INTEGER       N subsets to plot  [default: 1]
  --subsamples INTEGER    N samples per subset  [default: 100]
  --split FLOAT           % of split group to plot per subset  [default: 0.5]
  --cytoband TEXT         Cytoband file  [default: hg19]
  --cnp TEXT              Optional copy number segments
  --var-prob DICTIONARY   SV type prob thresholds  [default: INS: 0.2, DEL:
                          0.3, INV: 0.15, DUP: 0.15, TRA: 0.4]
  --var-su DICTIONARY     SV type supporting reads thresholds  [default: INS:
                          9, DEL: 9, INV: 6, DUP: 6, TRA: 8]
  --var-size INTEGER      Minimum SV size to plot (bp)  [default: 500]
  --help                  Show this message and exit.
```

## Quick start

```
# plot circos
python3 circos.py vcfs
# plot circos with probabilty filter and size filter
python3 circos.py --var-prob 'INS: 0.85, DEL: 0.85, INV: 0.45, DUP: 0.25, TRA: 0.35' --var-size 50000 vcfs
# plot circos with copy number ring
python3 circos.py --cnp output/segmented.copynumber.csv vcfs
# plot circos ordered/split by value (short)
python3 circos.py -l output/sample_pairs_short.csv -c short --cnp output/segmented.copynumber.csv vcfs
```

If there is a large cohort, multiple aggregated plots can also be made using the subsetting option `python3 circos.py --subsets 5 --subsamples 50 --ncols 10 vcfs` will generate 5 plots with 50 samples each (and 10 circos plots per row).

## Adapt for non-human reference

The line: `order_of_chr = [f"{i}" for i in range (1, 23)] + ["X", "Y"]` may have to be editted to include a list of the desired chromomsomes if hg38 or hg19 is not being used, note that the "chr" prefix should not be included as the string is removed later `name = line[0].replace("chr", "")`.

## How it works

Lengths are extracted from the vcf header using:
```
nl = {}
with open(fp) as f:
    f.readline()
    for line in f:
        if line.startswith("##contig"):
            name    = re.search("ID=([a-zA-Z0-9_]*)", line).group(1).replace("chr", "")
            length  = int(re.search("length=(\d+)", line).group(1))
            if "_" not in name: # can be removed if chromosomes including an underscore are wanted
                nl[name] = length
```

Vcf file parsed using pysam library, and variants filtered with input parameters:
```
f = pysam.VariantFile(fp)
tra_list = []
ins_dict = collections.defaultdict(dict)
for i in arcdata_dict.keys():
        ins_dict[i]["positions"] = []
        ins_dict[i]["widths"] = []
for l in f.fetch():
    name = l.chrom.replace("chr", "")
    sv_type = l.info["SVTYPE"]
    sample_name = list(l.samples.keys())[0]
    ## skip adding to plot if criteria not met
    if l.samples[sample_name]["PROB"] < prob_thresholds[sv_type]: # or l.samples[sample_name]["SU"] < type_su_thresholds[sv_type]:
        if l.samples[sample_name]["SU"] < su_thresholds[sv_type]:
            continue
    ## skip if variant not large enough
    start = l.pos
    end = l.stop #info["END"]
    if sv_type in {"DEL", "INV"}: #, "INS"}:
        if abs(end-start) <= size_threshold:
            continue
```

Use pycircos to generate citcos plot:
```
circle = pycircos.Gcircle(figsize=(5, 5))
...
for l in f.fetch():
    ...
    chr1 = l.contig.replace("chr", "")
    chr2 = l.info["CHR2"].replace("chr", "")
    if chr1 == chr2 and name in order_of_chr:
        source = (name, start, start, low_size-aof)
        destination = (name, end, end, low_size-aof)

        if abs(end-start) < 50000: # required for line to show
            destination = (name, end+50000, end+50000, low_size-aof) 
        if sv_type != "INS":
            circle.chord_plot(source, destination, facecolor=type_cols[sv_type], edgecolor=type_cols[sv_type], linewidth=0.4)
        else:
            ins_dict[str(name)]["positions"].append(start)
    elif chr1 != chr2 and name in order_of_chr and chr2 in order_of_chr: 
        name2 = l.info["CHR2"].replace("chr", "")
        start2 = l.info["CHR2_POS"]
        end2 = start2
        tra_list.append(((name, start, start, low_size-aof), (name2, start2, end2, low_size-aof)))
    else:
        continue
    for t in tra_list:
        circle.chord_plot(t[0], t[1], facecolor="#000000", linewidth=0.4)
    for key in ins_dict:
        if len(ins_dict[key]["positions"]) == 0:
        continue
    circle.scatterplot(key, data=[1]*len(ins_dict[key]["positions"]), positions=ins_dict[key]["positions"], 
                raxis_range=(low_size-aof, low_size-dof), facecolor="green", edgecolor="green",  markershape=".", markersize=20)
return circle
```

Stitch circos plots together using svgutils:
```
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

    short.append(shorts)
    if lengths != None:
        short.save(os.path.join(outdir, f"short_{z}.svg"))
```
