from libc.stdlib cimport malloc, free
from libc.string cimport strdup, strlen, memset, memcpy
from libc.stdio cimport printf
from cpython cimport array
import numpy as np
cimport numpy as cnp

cdef extern from "graph.h":
    struct node:
        unsigned int length
        unsigned long long *sequences
        unsigned short sequences_len
        unsigned long *edges
        unsigned char edges_len
        unsigned char reference

    struct graph:
        node *nodes
        unsigned long nodes_len

    void from_obgraph(graph *graph, unsigned int node_count,
                      unsigned char *sequences, unsigned int *sequence_lens,
                      unsigned int *edges, long long *edges_lens)
    int from_file(char *filepath, graph *graph)
    void to_file(char *filepath, graph *graph)
    int from_file_gfa_encoded(char *filepath, graph *graph, char *encoding, char flags)

    void set_encoding(graph *graph, char *encoding)

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
    void reverse_kmer_endian(kmer_finder *kf)

cdef class Graph:
    cdef graph data

    def __cinit__(self, encoding):
        self.data.nodes = NULL
        set_encoding(&(self.data), encoding)

    @staticmethod
    cdef void init_node(node *n, sequence, unsigned long sequence_len, edges, unsigned char edges_len, is_ascii):
        n.reference = 0
        n.length = sequence_len
        if sequence_len != 0:
            n.sequences_len = 1 + sequence_len // 32
        else:
            n.sequences_len = 0
        n.sequences = <unsigned long long *> malloc(n.sequences_len * sizeof(unsigned long long))
        for i in range(n.sequences_len):
            segment_end = min((i + 1) * 32, n.length)
            if is_ascii:
                n.sequences[i] = hash_max_kmer(sequence, i * 32, segment_end)
            else:
                n.sequences[i] = pack_max_kmer(sequence, i * 32, segment_end)
        n.edges_len = edges_len
        n.edges = <unsigned long *> malloc(edges_len * sizeof(unsigned long))
        for i in range(edges_len):
            n.edges[i] = edges[i]

    @staticmethod
    def from_obgraph(obg, encoding="ACGT"):
        if not Graph.is_valid_encoding(encoding):
            return None
        g = Graph(encoding.encode('ASCII'))
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
        """
        cdef cnp.ndarray[unsigned char, ndim=1, mode="c"] sequences = obg.sequences._data
        cdef cnp.ndarray[unsigned int, ndim=1, mode="c"] sequence_lens = obg.nodes
        cdef cnp.ndarray[unsigned int, ndim=1, mode="c"] edges = obg.edges._data
        cdef cnp.ndarray[long long, ndim=1, mode="c"] edges_len = obg.edges.shape.lengths.copy(order='C')
        from_obgraph(&(g.data), len(obg.nodes),
                     <unsigned char *> sequences.data,
                     <unsigned int *> sequence_lens.data,
                     <unsigned int *> edges.data,
                     <long long *> edges_len.data)
        cdef unsigned int i
        for i in obg.linear_ref_nodes():
            (g.data.nodes + i).reference = 1
        return g
        """

    @staticmethod
    def from_gfa(filepath, encoding="ACGT", optimize=True):
        cdef char flags = 0
        if optimize:
            flags += 1
        g = Graph(encoding.encode('ASCII'))
        cdef char *fpath = strdup(filepath.encode('ASCII'))
        from_file_gfa_encoded(fpath, &(g.data), encoding.encode('ASCII'), flags)
        free(fpath)
        return g

    @staticmethod
    def from_file(filepath):
        g = Graph("ACGT".encode('ASCII'))
        cdef char *fpath = strdup(filepath.encode('ASCII'))
        cdef int ret = from_file(fpath, &(g.data))
        free(fpath)
        if ret == 1:
            print("The specified file was of an invalid format.")
            raise
        return g

    def to_file(self, filepath):
        cdef char *fpath = strdup(filepath.encode('ASCII'))
        to_file(fpath, &(self.data))
        free(fpath)

    @staticmethod
    def from_sequence_edge_lists(sequences, edges, encoding="ACGT", ref=None):
        """
        Args:
            sequences: ["ACT", "G", "A", "GT"]
            edges: [[1, 2], [3], [3], []]
        List index determines a node's ID, and edges refer to what IDs are a node is connected to.
        Assumes the first node is the start node
        """
        if not Graph.is_valid_encoding(encoding):
            return None
        g = Graph(encoding.encode('ASCII'))
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

    @staticmethod
    def is_valid_encoding(encoding):
        if len(encoding) != 4:
            raise "Graph encoding must be a permutation of ACGT."
            return False
        encoding = encoding.upper()
        for i in "ACGT":
            if i not in encoding:
                raise "Graph encoding must be a permutation of ACGT."
                return False
        return True

    def create_kmer_index(self, k, max_variant_nodes=31, big_endian=True):
        if k < 1 or k > 31:
            raise "create_kmer_index: k must be between 1 and 31 inclusive"
        if max_variant_nodes <= 0 or max_variant_nodes > k:
            max_variant_nodes = k
        print("Finding kmers...")
        cdef kmer_finder *kf = init_kmer_finder(&(self.data), k, max_variant_nodes)
        find_kmers(kf)
        if not big_endian:
            reverse_kmer_endian(kf)
        print("Copying to numpy arrays")
        kmers = np.empty((kf.found_count,), dtype=np.ulonglong)
        nodes = np.empty((kf.found_count,), dtype=np.uint32)
        cdef cnp.ndarray[unsigned long long, ndim=1, mode="c"] c_kmers = kmers
        cdef cnp.ndarray[unsigned int, ndim=1, mode="c"] c_nodes = nodes
        memcpy(c_kmers.data, kf.found_kmers, sizeof(unsigned long long) * kf.found_count)
        memcpy(c_nodes.data, kf.found_nodes, sizeof(unsigned int) * kf.found_count)
        free_kmer_finder(kf)
        print("Done")
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
