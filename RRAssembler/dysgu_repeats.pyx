#!/bin/python3
#cython: language_level=3

#, boundscheck=False, wraparound=False

from libcpp.unordered_map cimport unordered_map
from math import exp
from xxhash cimport hash as xxhasher

cpdef float compute_rep(seq):
    cdef unordered_map[float, int] last_visited
    cdef float tot_amount = 0
    cdef float total_seen = 0
    cdef int k, i, diff
    cdef float decay, max_amount, amount

    cdef bytes s_bytes = bytes(seq.encode("ascii"))
    cdef const unsigned char* sub_ptr = s_bytes

    for k in (2, 3, 4, 5, 6, 7, 8, 9):

        decay = 0.25 * 1/k
        max_amount = exp(-decay) * k  # If last kmer was the same as current kmer

        sub_ptr = s_bytes
        for i in range(len(seq) - k):

            a = xxhasher(sub_ptr, k, 42)
            if last_visited.find(a) != last_visited.end():
                diff = i - last_visited[a]
                x = exp(-decay * diff)
                amount = (k * x) / max_amount

            else:
                amount = 0
            if i > k:
                tot_amount += amount
                total_seen += 1
            last_visited[a] = i
            sub_ptr += 1

    if total_seen == 0:
        return 0

    return tot_amount / total_seen
