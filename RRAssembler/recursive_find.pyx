#!/bin/python3
#cython: language_level=3

from pysam.libcalignmentfile cimport AlignmentFile, AlignedSegment
from collections import defaultdict, deque, Counter
import intervaltree
import random
import itertools
from cpython cimport array
cimport numpy as np


def get_reads(str region, AlignmentFile b):
    disc = []
    for aln in b.fetch(region=region):
        # -2 read mapped in proper pair
        # -256 not primary alignment
        if not aln.flag & 258:
            disc.append(aln)
        # Keep if read has more than 20 base pairs soft clipped
        elif aln.cigartuples:
            if len([i for i in aln.cigartuples if i[0] == 4 and i[1] > 20]) > 0:
                disc.append(aln)

    # Down sample region if necessary
    if len(disc) > 1000:  # Really high coverage
        random.shuffle(disc)
        disc = disc[:1000]  # Random downsample
        # disc = run_uclust(disc)

    return disc


def merge_intervals(intervals):
    # https://codereview.stackexchange.com/questions/69242/merging-overlapping-intervals
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: tup[0])
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append((higher, 1))
        else:
            lower = merged[-1][0]
            lower_count = merged[-1][1]
            if higher[0] <= lower[1]:
                upper_bound = max(lower[1], higher[1])
                merged[-1] = ((lower[0], upper_bound), lower_count + 1)  # replace by merged interval
            else:
                merged.append((higher, 1))
    # merged = [i for i in merged if i[1] > 1]
    return merged


def get_connecting_regions(list disc, AlignmentFile b, int pad=350):
    # Make an interval list for reads to get other regions of interest
    chroms = defaultdict(list)
    for aln in disc:
        chroms[b.get_reference_name(aln.rnext)].append((aln.pnext - pad, aln.pnext + pad))
    # Return only the region representation of merged regions
    # Drop low coverage regions, thats is regions with only 1 read in:
    return {k: [i[0] for i in merge_intervals(v) if i[1] > 1] for k, v in chroms.items() if merge_intervals(v)}


def update_visited(visited, chrom, begin, end, regions):
    # Add a new point to the visited and remove any overlaps from the regions list
    if not visited[chrom].overlap(begin, end): # originally search, now: at(point), overlap(b, e), envelop(b, e)
        visited[chrom].addi(begin, end, None)

    regions = {k: [i for i in v if not visited[k].overlap(*i)] for k, v in regions.items()}
    return regions


def get_names(dict regions):
    # Convert the regions dict into region names for pop
    names = []
    for k, v in regions.items():
        for start, end in v:
            names.append((k, start, end))
    return names


def recurive_find(AlignmentFile bam, str chrom, int start, int end, int padding, int max_depth=3):
    all_discordants = []
    # List of intervals which have been vistited, only visit ones which do not have overlaps
    visited = defaultdict(intervaltree.IntervalTree)
    names = [(chrom, start, end)]
    cdef int depth = 0
    cdef set seen = set([])
    while depth < max_depth:  # (Originally 100 for single region only)
        try:
            chrom, start, end = names.pop()
            start = 1 if int(start) < 1 else start  # Make sure start is not negative
            r = (chrom, start, end)
        except:
            # print("No pop")
            break
        if r in seen:
            continue
        seen.add(r)

        disc = get_reads("{}:{}-{}".format(*r), bam)
        all_discordants += disc
        new_r = get_connecting_regions(disc, bam, pad=padding)  # 350
        new_r = update_visited(visited, r[0], r[1], r[2], new_r)
        new_names = get_names(new_r)
        names = new_names + names  # Keep old names at the end to do breadth first, not depth first
        depth += 1

    return all_discordants

def filter_duplicates(reads):
    cdef set seen_reads = set([])
    unique = []
    for item in reads:
        try:
            id_tup = item.qname, item.flag  # Each read should have a unique name and flag
        except:
            print("Read has improper name or flag")
            continue
        if id_tup in seen_reads:
            continue
        else:
            unique.append(item)
            seen_reads.add(id_tup)
    return unique


def sort_disc(unique):
    pairs = []
    singles = []
    s = itertools.groupby(sorted(unique, key=lambda x: x.qname), key=lambda x: x.qname)
    ls = []
    for k, v in s:
        l = list(v)
        ls.append(len(l))
        if len(l) == 2:
            if not l[1].flag & 64:  # Make sure first in pair comes first
                l = l[::-1]
            pairs += l
        else:
            singles += l

    print("Singles:", len(singles), " Pairs:", len(pairs), Counter(ls))
    return singles, pairs
