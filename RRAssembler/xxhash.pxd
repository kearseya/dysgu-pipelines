from libc.stdint cimport uint64_t

cdef extern from "xxhash64.h" namespace "XXHash64" nogil:
    cdef uint64_t hash(void* input, uint64_t length, uint64_t seed) nogil
