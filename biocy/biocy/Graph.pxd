from biocy.node cimport node
from libc.stdint cimport uint8_t, uint32_t

cdef extern from "cpp/Graph.hpp":
    cdef cppclass Graph:
        node *nodes
        uint32_t nodes_len
        char encoding[4]
        uint8_t encoding_map[256]

        @staticmethod
        Graph *FromFile(char *)
        @staticmethod
        Graph *FromGFAFile(char *)
        @staticmethod
        Graph *FromGFAFileEncoded(char *, char *)

        void Compress()

        void ToFile(char *)
