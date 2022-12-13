from Graph cimport Graph
from libc.stdint cimport uint8_t, uint32_t, uint64_t

cdef extern from "cpp/KmerFinder.hpp":
    cdef cppclass KmerFinder:
        KmerFinder(Graph *, uint8_t, uint8_t) except +
        Graph *graph
        const uint8_t k
        const uint8_t max_variant_nodes
        uint64_t *found_kmers
        uint32_t *found_nodes
        uint64_t found_count

        void Find()
        void ReverseFoundKmers()
