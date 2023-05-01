from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t
from libcpp.unordered_map cimport unordered_map

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
        uint32_t reference_index
        bool reference

cdef extern from "cpp/Graph.hpp":
    cdef cppclass Graph:
        Graph(char *) except +
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

        uint64_t HashMinKmer(char *, uint8_t)
        uint64_t HashMaxKmer(char *, uint8_t)
        uint64_t HashKmer(char *, uint8_t)
        char *DecodeKmer(uint64_t, uint8_t)

        void AddInEdges()

        uint32_t GetNextReferenceNodeID(uint32_t)

cdef extern from "cpp/KmerFinder.hpp":
    enum: FILTER_NODE_ID
    enum: FLAG_ALIGN_SIGNATURE_WINDOWS
    enum: FLAG_MINIMIZE_SIGNATURE_OVERLAP
    enum: FLAG_ONLY_SAVE_INITIAL_NODES
    enum: FLAG_SAVE_WINDOWS

    cdef cppclass VariantWindow:
        uint64_t *reference_kmers;
        uint64_t *variant_kmers;
        uint16_t reference_kmers_len;
        uint16_t variant_kmers_len;
        uint32_t max_frequency;

        void ReverseKmers(uint8_t k)

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

        void SetFilter(uint8_t, uint64_t)
        void RemoveFilter(uint8_t)
        void SetFlag(uint8_t, bool)
       
        unordered_map[uint64_t, uint32_t] CreateKmerFrequencyIndex()
        void SetKmerFrequencyIndex(unordered_map[uint64_t, uint32_t])
        bool HasKmerFrequencyIndex()

        KmerFinder *CreateWindowFinder()
        VariantWindow *FindVariantSignatures(uint32_t reference_node_id, uint32_t variant_node_id)
        VariantWindow *FindVariantSignaturesWithFinder(
                uint32_t reference_node_id, uint32_t variant_node_id, KmerFinder *kf)

