from libc.stdlib cimport malloc, free
from libc.string cimport strdup, strlen, memset
from libc.stdio cimport printf
from cpython cimport array
import numpy as np

cdef extern from "graph.h":
    struct node:
        unsigned int len
        unsigned long long *sequences
        unsigned short sequences_len
        unsigned long *edges
        unsigned char edges_len
        unsigned char reference

    struct graph:
        node *nodes
        unsigned long nodes_len

cdef extern from "kmer_finder.h":
    struct kmer_finder:
        unsigned char k
        unsigned char max_variant_nodes

        unsigned long long *found_kmers
        unsigned long *found_nodes
        unsigned long long found_len
        unsigned long long found_count

        unsigned long long kmer_buffer
        unsigned long long kmer_buffer_ext

        unsigned long *path_buffer
        unsigned char path_buffer_len

        unsigned long long kmer_mask
        unsigned char shift

        unsigned char variant_counter

    kmer_finder *init_kmer_finder(graph *graph, unsigned char k, unsigned char max_variant_nodes)
    void free_kmer_finder(kmer_finder *kf)
    void find_kmers(kmer_finder *kf)

#ctypedef node_t Node_t
#ctypedef graph_t Graph_t

#cdef class Node:
#    cdef Node_t data
#
#    def __cinit__(self, char *sequence, unsigned char edge_count):
#        self.data.sequence = <char *> malloc(len(sequence) * sizeof(char))
#        self.data.edge_count = edge_count
#        self.data.edges = <unsigned int *> malloc(edge_count * sizeof(char))
#
#    def __dealloc__(self):
#        free(self.data.sequence)
#        free(self.data.edges)

cdef class Graph:
    cdef graph data

    def __cinit__(self):
        self.data.nodes = NULL

    """
    @staticmethod
    def from_sequence_edges(sequence, edges):
        node = Node()
        node.sequence_len = len(sequence)
        node.num_subsequences = 0 if node.sequence_len == 0 else 1 + node.sequence_len // 31
        node.sequence = np.empty((node.num_subsequences,), dtype=np.ulonglong)
        for i in range(node.num_subsequences):
            segment_end = min((i + 1) * 31, node.sequence_len)
            node.sequence[i] = hash_max_kmer(sequence.encode('ASCII')[(i * 31):segment_end], segment_end - (i * 31))
        node.edges = np.array(edges, dtype=np.int32)
        return node

    @staticmethod
    def from_obgraph(node, sequence, edges):
        node = Node()
        node.sequence_len = sequence.shape[0]
        node.num_subsequences = 0 if node.sequence_len == 0 else 1 + node.sequence_len // 31
        node.sequence = np.empty((node.num_subsequences,), dtype=np.ulonglong)
        for i in range(node.num_subsequences):
            segment_end = min((i + 1) * 31, node.sequence_len)
            node.sequence[i] = pack_max_kmer(sequence, i * 31, segment_end)
        node.edges = np.array(edges, dtype=np.int32)
        return node
    """

    @staticmethod
    cdef void init_node(node *n, sequence, unsigned int sequence_len, edges, unsigned char edges_len, is_ascii):
        n.reference = 0
        n.len = sequence_len
        if sequence_len != 0:
            n.sequences_len = 1 + sequence_len // 31
        else:
            n.sequences_len = 0
        n.sequences = <unsigned long long *> malloc(n.sequences_len * sizeof(unsigned long long))
        for i in range(n.sequences_len):
            segment_end = min((i + 1) * 31, n.len)
            if is_ascii:
                n.sequences[i] = hash_max_kmer(sequence, i * 31, segment_end)
            else:
                n.sequences[i] = pack_max_kmer(sequence, i * 31, segment_end)
        n.edges_len = edges_len
        n.edges = <unsigned long *> malloc(edges_len * sizeof(unsigned long))
        for i in range(edges_len):
            n.edges[i] = edges[i]

    @staticmethod
    def from_obgraph(obg):
        g = Graph()
        cdef node *n
        cdef unsigned int node_count = len(obg.nodes)
        cdef unsigned int i
        g.data.nodes = <node *> malloc(node_count * sizeof(node))
        g.data.nodes_len = node_count
        for i in range(node_count):
            n = g.data.nodes + i
            Graph.init_node(n,
                            obg.sequences[i],
                            obg.sequences[i].shape[0],
                            obg.edges[i],
                            obg.edges[i].shape[0],
                            False)
        ref = obg.linear_ref_nodes()
        for i in ref:
            (g.data.nodes + i).reference = 1
    
        return g

    @staticmethod
    def from_sequence_edge_lists(sequences, edges, ref=None):
        """
        Args:
            sequences: ["ACT", "G", "A", "GT"]
            edges: [[1, 2], [3], [3], []]
        List index determines a node's ID, and edges refer to what IDs are a node is connected to.
        Assumes the first node is the start node
        """
        g = Graph()
        cdef node *n
        cdef unsigned int node_count = len(sequences)
        cdef unsigned int i
        g.data.nodes = <node *> malloc(node_count * sizeof(node))
        g.data.nodes_len = node_count
        for i in range(node_count):
            n = g.data.nodes + i
            Graph.init_node(n,
                            strdup(sequences[i].encode('ASCII')),
                            len(sequences[i]),
                            edges[i],
                            len(edges[i]),
                            True)
        if ref is not None:
            for i in ref:
                (g.data.nodes + i).reference = 1
        else:
            for i in range(node_count):
                (g.data.nodes + i).reference = 1
        return g


    def create_kmer_index(self, k, max_variant_nodes=4):
        cdef kmer_finder *kf = init_kmer_finder(&self.data, k, max_variant_nodes)
        find_kmers(kf)
        #kmers, nodes = GraphKmerFinder(self, k, max_variant_nodes).find()
        kmers = np.empty((kf.found_count,), dtype=np.ulonglong)
        nodes = np.empty((kf.found_count,), dtype=np.uint32)
        for i in range(kf.found_count):
            kmers[i] = kf.found_kmers[i]
            nodes[i] = kf.found_nodes[i]
        free_kmer_finder(kf)
        return kmers, nodes

def hash_kmer(arr, k):
    cdef unsigned long long hashed = 0
    for i in range(k):
        hashed |= (arr[i] & 6) << ((k - i - 1) << 1)
    return hashed >> 1

cdef unsigned long long hash_max_kmer(arr, unsigned int start, unsigned int end):
    cdef unsigned long long hashed = 0
    cdef unsigned long long val
    for i in range(start, end):
        val = arr[i]
        hashed |= (val & 6) << (61 - (i - start) * 2)
    return hashed

cdef unsigned long long pack_max_kmer(arr, unsigned int start, unsigned int end):
    cdef unsigned long long packed = 0
    cdef unsigned long long val
    for i in range(start, end):
        val = arr[i]
        packed |= (val << (62 - (i - start) * 2))
    return packed
