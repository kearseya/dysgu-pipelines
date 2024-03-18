# Chain link

For each chromosome, the frequency of break-sites is calculated across all samples, then for each sample and each chromosome a K-nearest neighbour tree is made. Then all breakpoints are iterated over checking for links where a binomial probability mass function gives a p-value less than a threshold, if they are then an edge is made on a graph. Clusters in the graph that reach above the threshold for number of breakpoints are then processed and saved. 

Algorithm based on https://www.sciencedirect.com/science/article/pii/S0092867413003437

## Usage

```
Usage: chain_link_finder.py [OPTIONS] INDIR COMMAND [ARGS]...

Options:
  -o, --outdir TEXT
  -p, --prob FLOAT
  -b, --breaks INTEGER
  -l, --lower-p FLOAT
  -u, --upper-p FLOAT
  -n, --np-tests INTEGER
  -t, --tl TEXT
  -s, --sep-thresh FLOAT
  -e, --sample-col TEXT
  -a, --stela-col TEXT
  -c, --cytofile TEXT
  --chrom-len TEXT
  --size-thresh INTEGER   SV size filter
  --var-prob DICTIONARY   SV type prob thresholds
  --help                  Show this message and exit.

Commands:
  plot
  run
```


## Calculate frequency

```
ps = {str(k): 0.0 for k in cl.keys()}  # Sample: avg_breaks per bp across all samples
for key, value in breaks_list.items():
	break_n_per_chrom = defaultdict(int)

	for chrom, b in value.items():  # For each break in the sample
		break_n_per_chrom[str(chrom)] += len(b)
	avg = {str(k): v / cl[str(k)] for k, v in break_n_per_chrom.items()}
	for k, v in break_n_per_chrom.items():
		ps[k] += v


ps = {k: (v/len(breaks_list) / cl[k]) for k, v in ps.items()}
```

## Make nearest neighbour tree

```
nn = defaultdict(lambda: defaultdict(list))  # sample: chromosome: neighbour tree
for samp, value in breaks_list.items():
	for chrom, bs in value.items():
		X = np.array(bs)[:, np.newaxis]
		tree = KDTree(X)
		nn[samp][chrom] = tree

```

## Finding connected

The p-value is calculated with: `pval = 1 - stats.binom.pmf(0, fails, prob_success)` where hat 1 or more breaksites are observed within this given distance. The distance in bp is the number of fails. i.e. the number of bp without a break, success is 0, the probability that no breaks are oberved.


The list of breakpoints for each sample is iterated over. For each breakpoint look at up- and downstream sites, if the p-value is less than threshold, make an edge on the graph.

```
def find_connected(tree, data, chrom, pos, p_val):
	dist, ind = tree.query([[pos]], k=len(data))  # k=len(data) to enumerate all neighbours
	pvals = [get_pval(i, ps[chrom]) for i in dist[0]]
	dist, ind = dist.flatten(), ind.flatten()
	c = set([])
	for i, p in enumerate(pvals[1:]):  # Skip first entry which is the current SV
		if p < p_val:
			c.add(("{}:{}".format(chrom, int(pos)), "{}:{}".format(chrom, int(data[ind[i + 1]]))))
	return c
```

Which is used to create a set of connected SVs:

```
def neigh(s, samp, nn, prob):
	tree1 = nn[samp][s[0]]  # Tree for chromosome 1
	tree2 = nn[samp][s[3]]  # Tree for chromosome 2, (intra or inter)

	data1 = np.array(tree1.data).flatten()
	data2 = np.array(tree2.data).flatten()

	c1 = find_connected(tree1, data1, s[0], s[1], prob)
	c2 = find_connected(tree2, data2, s[3], s[5], prob)

	connected = set([])
	connected |= c1  # |= ior
	connected |= c2

	return connected
```

## Write clusters to file

```
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
```


