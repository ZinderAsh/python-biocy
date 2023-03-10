# distutils: language = c++

from libc.stdlib cimport malloc, free
from libc.string cimport strdup, strlen, memset, memcpy
from libc.stdio cimport printf
from libc.stdint cimport uint8_t, uint32_t, uint64_t
from cpython cimport array
import numpy as np
cimport numpy as cnp
from npstructures import RaggedArray

cimport biocy.biocpp as cpp
from biocy.hashing cimport pack_max_kmer_with_offset, decode_kmer_by_map, fill_map_by_encoding, hash_min_kmer_by_encoding

cdef uint8_t get_node_base(cpp.node *n, uint32_t index):
    cdef uint32_t sequence_index = index // 32
    cdef uint8_t sub_index = index % 32

    return ((n.sequences[sequence_index]) >> ((31 - sub_index) * 2)) & 3

cdef class Graph:
    cdef cpp.Graph *data

    def __cinit__(self):
        self.data = NULL

    def __dealloc__(self):
        del self.data

    cdef void init_node(self, cpp.node *n, sequence, uint32_t sequence_len, edges, uint8_t edges_len, is_ascii):
        cdef char *ascii_seq
        cdef cnp.ndarray[char, ndim=1, mode="c"] numpy_seq
        cdef uint32_t i
        cdef uint64_t hashed_seq
        n.reference = 0
        n.reference_index = 0
        n.length = sequence_len
        if sequence_len != 0:
            n.sequences_len = 1 + sequence_len // 32
        else:
            n.sequences_len = 0
        n.sequences = <uint64_t *> malloc(n.sequences_len * sizeof(uint64_t))
        for i in range(n.sequences_len):
            segment_end = min((i + 1) * 32, n.length)
            if is_ascii:
                ascii_seq = strdup(sequence)
                n.sequences[i] = self.data.HashMaxKmer(ascii_seq + (i * 32), segment_end - (i * 32))
                free(ascii_seq)
            else:
                numpy_seq = sequence
                n.sequences[i] = pack_max_kmer_with_offset(numpy_seq.data, i * 32, segment_end - (i * 32))
        n.edges_len = edges_len
        n.edges = <uint32_t *> malloc(edges_len * sizeof(uint32_t))
        for i in range(edges_len):
            n.edges[i] = edges[i]
        n.edges_in_len = 0
        n.edges_in = NULL
    
    @staticmethod
    def from_obgraph(obg, encoding="ACGT"):
        if not has_obgraph:
            raise Exception("Graph.from_obgraph requires the 'obgraph' module to function.")
        if not Graph.is_valid_encoding(encoding):
            return None
        g = Graph()
        g.data = new cpp.Graph(encoding.encode('ASCII'))
        cdef cpp.node *n
        cdef uint32_t node_count = len(obg.nodes)
        cdef uint32_t edge_count
        cdef uint32_t edge_id
        cdef uint32_t i
        cdef uint32_t j
        cdef uint32_t index
        cdef uint8_t byte
        g.data.nodes = <cpp.node *> malloc(node_count * sizeof(cpp.node))
        g.data.nodes_len = node_count
        for i in range(node_count):
            n = g.data.nodes + i
            g.init_node(n,
                        obg.sequences[i],
                        obg.sequences[i].shape[0],
                        obg.edges[i],
                        obg.edges[i].shape[0],
                        False)
        g.data.AddInEdges()
        ref = obg.linear_ref_nodes_and_dummy_nodes_index
        for i, byte in enumerate(ref):
            (g.data.nodes + i).reference = byte

        return g

    def to_obgraph(self):
        if not has_obgraph:
            raise Exception("Graph.to_obgraph requires the 'obgraph' module to function.")
        cdef cpp.node *n
        cdef uint32_t i
        cdef uint32_t j
        cdef uint64_t reference_length = 0
        cdef uint32_t root_id
        cdef cpp.node *root

        node_count = self.data.nodes_len

        node_lengths = np.empty((node_count,), dtype=np.uint32)
        edge_lengths = np.empty((node_count,), dtype=np.uint8)
        chromosome_start_nodes = None
        for i in range(node_count):
            n = self.data.nodes + i
            node_lengths[i] = n.length
            edge_lengths[i] = n.edges_len

            if chromosome_start_nodes is None:
                if n.edges_len > 0 and n.edges_in_len == 0:
                    if not self.node_has_parent(i):
                        chromosome_start_nodes = [i]

        sequences_val = np.empty((sum(node_lengths),), dtype=np.uint8)
        edges_val = np.empty((sum(edge_lengths),), dtype=np.uint32)
        cdef uint64_t base_index = 0
        cdef uint64_t edge_index = 0;
        for i in range(node_count):
            n = self.data.nodes + i
            for j in range(n.length):
                sequences_val[base_index] = get_node_base(n, j)
                base_index += 1
            for j in range(n.edges_len):
                edges_val[edge_index] = n.edges[j]
                edge_index += 1

        sequences = RaggedArray(sequences_val, node_lengths)
        edges = RaggedArray(edges_val, edge_lengths)

        linear_ref_nodes_index = np.zeros((node_count,), dtype=np.uint8)
        linear_ref_nodes_and_dummy_nodes_index = np.zeros((node_count,), dtype=np.uint8)
        for i in range(node_count):
            n = self.data.nodes + i
            if n.reference:
                linear_ref_nodes_and_dummy_nodes_index[i] = 1
                if n.length > 0:
                    linear_ref_nodes_index[i] = 1

        node_to_ref_offset = None
        ref_to_node_offset = None
        if False and chromosome_start_nodes is not None:
            #node_to_ref_offset = np.zeros((node_count,), dtype=np.uint64)
            root_id = chromosome_start_nodes[0]
            root = self.data.nodes + root_id
            count = 0
            while root.edges_len > 0:
                for i in range(root.edges_len):
                    n = self.data.nodes + root.edges[i]
                    print(root_id, root.reference_index, "->", root.edges[i], n.reference_index)
                root_id = self.data.GetNextReferenceNodeID(root_id)
                root = self.data.nodes + root_id
                count += 1
            print(count)

        return ob.Graph(
            node_lengths,
            sequences,
            edges,
            node_to_ref_offset=node_to_ref_offset,
            ref_offset_to_node=ref_to_node_offset,
            chromosome_start_nodes=chromosome_start_nodes,
            linear_ref_nodes_index=linear_ref_nodes_index,
            linear_ref_nodes_and_dummy_nodes_index=linear_ref_nodes_and_dummy_nodes_index
        )

    cdef node_has_parent(self, node_id):
        cdef uint32_t i
        cdef uint32_t j
        for i in range(self.data.nodes_len):
            parent = self.data.nodes + i
            for j in range(parent.edges_len):
                if parent.edges[j] == node_id:
                    return True
        return False

    @staticmethod
    def from_gfa(filepath, encoding="ACGT", compress=True):
        cdef char flags = 0
        cdef char *fpath = strdup(filepath.encode('ASCII'))
        cdef cpp.Graph *cpp_graph = cpp.Graph.FromGFAFileEncoded(fpath, encoding.encode('ASCII'))
        if compress:
            cpp_graph.Compress()
        g = Graph()
        g.data = cpp_graph
        free(fpath)
        return g

    @staticmethod
    def from_file(filepath):
        cdef char *fpath = strdup(filepath.encode('ASCII'))
        cdef cpp.Graph *cpp_graph = cpp.Graph.FromFile(fpath)
        free(fpath)
        if cpp_graph == NULL:
            print("The specified file was of an invalid format.")
            raise
        g = Graph()
        g.data = cpp_graph
        return g

    def to_file(self, filepath):
        cdef char *fpath = strdup(filepath.encode('ASCII'))
        self.data.ToFile(fpath)
        free(fpath)

    def print_node_data(self, node_id):
        cdef uint32_t i = node_id
        cdef uint8_t j
        cdef cpp.node *n = self.data.nodes + i
        output = "Node ID: " + str(node_id)
        output += "\nLength: " + str(n.length)
        output += "\nEdges Out Len: " + str(n.edges_len)
        output += "\nEdges Out:"
        for j in range(n.edges_len):
            output += " " + str(n.edges[j])
        output += "\nEdges In Len: " + str(n.edges_in_len)
        output += "\nEdges In:"
        for j in range(n.edges_in_len):
            output += " " + str(n.edges_in[j])
        print(output)


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
        cdef cpp.node *n
        cdef unsigned int node_count = len(sequences)
        cdef unsigned int i
        cdef char *sequence
        g.data = new cpp.Graph(encoding.encode('ASCII'))
        g.data.nodes = <cpp.node *> malloc(node_count * sizeof(cpp.node))
        g.data.nodes_len = node_count
        for i in range(node_count):
            n = g.data.nodes + i
            sequence = strdup(sequences[i].encode('ASCII'))
            g.init_node(n,
                        sequence,
                        len(sequences[i]),
                        edges[i],
                        len(edges[i]),
                        True)
            free(sequence)
        g.data.AddInEdges()
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
            raise "Graph encoding must be a permutation of ACGT (length was not 4)."
        encoding = encoding.upper()
        for i in "ACGT":
            if i not in encoding:
                raise "Graph encoding must be a permutation of ACGT (did not include all characters)."
        return True

def hash_kmer(kmer, k, encoding="ACGT"):
    return hash_min_kmer_by_encoding(kmer.encode('ASCII'), k, encoding.encode('ASCII'))
