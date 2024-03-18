#!/bin/python3
#cython: language_level=3

#, boundscheck=False, wraparound=False

import array
from cpython cimport array
import numpy as np
cimport numpy as np

#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#DTYPE = np.int
#ctypedef np.int_t DTYPE_t

DTYPE = float
ctypedef float DTYPE_t


def optimal_path(np.ndarray[DTYPE_t, ndim=2] segments, float contig_length, float max_insertion=1500, float min_aln=20,
                 float max_homology=1500, float ins_cost=2.0, float hom_cost=1.5, float inter_cost=20, float intra_cost=10):
    """
    The scoring has become quite complicated.
     = Last alignment score - jump cost

    where each jump has a cost:
    jump cost = S + microhomology_length + insertion_length        # S = 10 for intra, 20 for inter

    The last score takes into account mismatched and gaps

    :param mapped:  The input dataframe with candidate alignments
    :param contig_length:   Length of the contig
    :param max_insertion: Arbitrary high number
    :param min_aln: The minimum sequence which is not part of a microhomology overlap
    :param max_homology: Arbitrary high number
    :return: A dataframe with the inferred alignments
    """

    # Start at first node then start another loop running backwards through the preceeding nodes.
    # Choose the best score out of the in edges.

    # Use a special start and end node to score the first and last alignments
    cdef array.array int_array_template = array.array('i', [])
    cdef array.array float_array_template = array.array('f', [])

    pred = array.clone(int_array_template, len(segments), zero=True)
    node_scores = array.clone(float_array_template, len(segments), zero=True)

    cdef int i, j, p, end_i, next_i
    cdef float current_chrom, current_start, current_end, next_chrom, next_start,\
        next_end, jump_cost, next_score, micro_h, ins, best_score, current_score, sc, max_s

    # Deal with first score
    node_scores[0] = segments[0, 3] - (segments[0, 1] * ins_cost)
    pred[0] = -1


    # start from segment two because the first has been scored
    for i in xrange(1, len(segments)):
        next_chrom = segments[i, 0]
        next_start = segments[i, 1]
        next_end = segments[i, 2]
        next_score = segments[i, 3]

        #next_chrom, next_start, next_end, next_score = segments[i]

        p = -1                                  # Ins_cost?
        best_score = next_score - (next_start * ins_cost)  # Implies all preceding alignments skipped!

        # Walking backwards mean the search may be terminated at some point
        for j in xrange(i-1, -1, -1):

            current_chrom = segments[j, 0]
            current_start = segments[j, 1]
            current_end = segments[j, 2]
            current_score = segments[j, 3]

            # Allow alignments with minimum sequence and max overlap
            if next_start > current_end - max_homology and next_end > current_end + min_aln and \
                                    next_start - current_start > min_aln:

                if next_start > current_end and next_start - current_end > max_insertion:
                    continue

                micro_h = current_end - next_start #seg[1]
                if micro_h < 0:
                    ins = abs(micro_h)
                    micro_h = 0
                else:
                    ins = 0

                if next_chrom == current_chrom:
                    jump_cost = intra_cost
                else:
                    jump_cost = inter_cost

                # Calculate score, last_score = node_scores[j]
                current_score = node_scores[j] - (micro_h * hom_cost) - (ins * ins_cost) - jump_cost + next_score

                if current_score > best_score:
                    best_score = current_score
                    p = j

        node_scores[i] = best_score
        pred[i] = p

    # Update the score for jumping to the end of the sequence
    max_s = -1e6
    end_i = 0
    for i in range(len(segments)):
        sc = node_scores[i] - (ins_cost * (contig_length - segments[i, 2]))
        # node_scores[i] = node_scores[i] - (ins_cost * (contig_length - segments[i][2]))
        if sc > max_s:
            max_s = sc
            end_i = i

    # Get path from max
    indexes = [end_i] #[node_scores.index(max(node_scores))]

    while True:
        # Use len(indexes) - 1 to get around wraparound constraint
        next_i = pred[indexes[len(indexes) - 1]]
        if next_i == -1:
            break
        indexes.append(next_i)

    return indexes[::-1], max_s

