from libc.stdint cimport uint8_t, uint32_t, uint64_t
from libcpp cimport bool

cdef extern from "cpp/node.hpp":
    struct node:
        uint32_t length
        uint64_t *sequences
        uint32_t sequences_len
        uint32_t *edges
        uint32_t *edges_in
        uint8_t edges_len
        uint8_t edges_in_len
        bool reference
