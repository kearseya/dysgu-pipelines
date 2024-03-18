# Contig assembly by recursion

## Usage

```
Usage: RRAssembler.py [OPTIONS] COMMAND [ARGS]...

Options:
  -b, --bam-glob TEXT
  -c, --chain-file TEXT
  -r, --ref TEXT
  -p, --pad INTEGER
  --help                 Show this message and exit.

Commands:
  analyse  3rd analyse contigs
  collect  5th collet read map
  fetch    1st fetch and assemble contigs
  map      2nd map reads to contigs
  plot     6th plot contigs
  remove   4th remove common sites
```

## Fetch

Fetching reads is performed with:

```
sample_bam = bam_paths[sample]
bam = pysam.AlignmentFile(sample_bam, "rb")

out_singles = pysam.AlignmentFile("./{o}/{s}/{s}_singles.bam".format(o=OUT, s=sample), "wb", template=bam)
out_pairs = pysam.AlignmentFile("./{o}/{s}/{s}_paired.bam".format(o=OUT, s=sample), "wb", template=bam)

with Progress(TextColumn("[progress.description]{task.description}"), BarColumn(), TaskProgressColumn(), TimeRemainingColumn(), TimeElapsedColumn()) as progress:
	task = progress.add_task("[red]Recursive serach...", total=len(df)*2)
	discordants = []
	for idx, r in df.iterrows():
		d1 = recursive_find.recurive_find(bam, r.chrom1, r.chr1_start, r.chr1_end, pad, max_depth=3)
		progress.update(task, advance=1)
		d2 = recursive_find.recurive_find(bam, r.chrom2, r.chr2_start, r.chr2_end, pad, max_depth=3)
		discordants += recursive_find.filter_duplicates(d1+d2)
		progress.update(task, advance=1)

unique = recursive_find.filter_duplicates(discordants)
print("Reads to assemble", len(unique))

# Sort reads by Q name and split into singles and paired reads for best assembly
singles, pairs = recursive_find.sort_disc(unique)

# Write a fastQ file
with out_singles, out_pairs:
	for read in unique:   # WRITE OUT ALL READS, NOT JUST SINGLES
		out_singles.write(read)
	for read in pairs:
		out_pairs.write(read)
```

This output is then used in assembly with spades, by calling `assemble_and_map.sh` passing the variables via the command line.

## Map

Mapping is performed in the `map_reads_to_contigs.py` script which calls `bwa mem` to perform the mapping.

## Analyse

This portition uses dodi, a library written by Kez Cleal https://github.com/kcleal/dodi. Dodi chooses an optimal set of alignments from a candidate list which can be used for improving split-read mapping, or for selectively mapping to target regions of the genome.

## Remove

Contigs may represent common SVs which are not unique to the post-crisis clones which are filtered out

## Collect

Adjacent alignments in the contig are compared for simularity. Sections of microhomology are always included in the comparison. For insertions, the insertion is only added to one of the adjacent sequences. Insertions from which fall before the first alignment or after the second alignment are ignored. The decision of which alignment the insertion is joined to is governed by the highest simularity score.

## Plot

Contigs for each sample are plotted as a stacked broken-hbar plot with colours to indicate the chromosome diferent sections of the contig are mapped to.
