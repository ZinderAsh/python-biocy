from libc.stdlib cimport malloc, free
from libc.string cimport strdup, strlen, memset
from libc.stdio cimport printf
from cpython cimport array

cdef struct node_t:
    unsigned int id
    char *sequence
    unsigned char edge_count
    unsigned int *edges
    char visited
ctypedef node_t Node_t

cdef struct graph_t:
    Node_t *nodes
    int node_count
ctypedef graph_t Graph_t

cdef class Node:
    cdef Node_t data

    def __cinit__(self, char *sequence, unsigned char edge_count):
        self.data.sequence = <char *> malloc(len(sequence) * sizeof(char))
        self.data.edge_count = edge_count
        self.data.edges = <unsigned int *> malloc(edge_count * sizeof(char))

    def __dealloc__(self):
        free(self.data.sequence)
        free(self.data.edges)

cdef class Graph:
    cdef Graph_t data

    def __cinit__(self):
        self.data.nodes = NULL

    def root(self):
        return None #Node(self.data.start_node.sequence, self.data.start_node.edge_count)

    @staticmethod
    def from_sequence_edge_lists(sequences, edges):
        """
        Args:
            sequences: ["ACT", "G", "A", "GT"]
            edges: [[1, 2], [3], [3], []]
        List index determines a node's ID, and edges refer to what IDs are a node is connected to.
        Assumes the first node is the start node
        """
        cdef unsigned int node_count = len(edges)
        cdef Node_t *node
        cdef unsigned int edge_count
        cdef unsigned int i, j
        graph = Graph()
        graph.data.nodes = <Node_t *> malloc(node_count * sizeof(Node_t))
        graph.data.node_count = node_count
        for i in range(node_count):
            node = graph.data.nodes + i
            edge_count = len(edges[i])
            node.sequence = strdup(sequences[i].encode('ASCII'))
            node.edge_count = edge_count
            node.visited = 0
            if edge_count > 0:
                node.edges = <unsigned int *> malloc(edge_count * sizeof(unsigned int))
                for j in range(edge_count):
                    node.edges[j] = edges[i][j]
            else:
                node.edges = NULL

        return graph

    cdef unsigned int get_kmers(self, char *kmer, unsigned char kmer_len, int node_start_index, int node_index, unsigned char k, dict index):
        cdef Node_t *node = self.data.nodes + node_index
        cdef unsigned int i
        cdef unsigned int j
        cdef unsigned int seq_len = strlen(node.sequence)
        cdef unsigned char new_kmer_len = 0
        cdef char extending = (node_start_index != node_index)
        if not extending:
            node.visited = 1
        for i in range(seq_len):
            kmer[kmer_len] = node.sequence[i]
            if kmer_len >= k - 1:
                if not kmer + kmer_len - k + 1 in index:
                    index[kmer + kmer_len - k + 1] = []
                index[kmer + kmer_len - k + 1].append(node_start_index)
            if extending:
                if kmer_len >= k - 1:
                    return kmer_len + 1
                else:
                    kmer_len += 1
            elif kmer_len + 1 == k * 2: 
                for j in range(k * 2 - 2):
                    kmer[j] = kmer[j + 1]
            else:
                kmer_len += 1
        if extending:
            for i in range(node.edge_count):
                new_kmer_len = self.get_kmers(kmer, kmer_len, node_start_index, node.edges[i], k, index)
                while (new_kmer_len > kmer_len):
                    new_kmer_len -= 1
                    kmer[new_kmer_len] = 0
            return kmer_len
        for i in range(0 if kmer_len >= k else k - 1 - kmer_len, k):
            if kmer_len >= k - i:
                if i == k - 1:
                    kmer[0] = 0
                    kmer_len = 0
                elif i > 0:
                    for j in range(k - i - 1):
                        kmer[j] = kmer[j + 1]
                    kmer[k - i - 1] = 0
                    kmer_len = k - i - 1
                else:
                    for j in range(k - 1):
                        kmer[j] = kmer[j + kmer_len - k + 1]
                    for j in range(k - 1, kmer_len):
                        kmer[j] = 0
                    kmer_len = k - 1
            for j in range(node.edge_count):
                if kmer_len > 0 or not self.data.nodes[node.edges[j]].visited:
                    new_kmer_len = self.get_kmers(kmer, kmer_len, node.edges[j] if kmer_len == 0 else node_start_index, node.edges[j], k, index)
                    while (new_kmer_len > kmer_len):
                        new_kmer_len -= 1
                        kmer[new_kmer_len] = 0
        return kmer_len

    def create_kmer_index(self, int k):
        index = {}
        cdef int i
        if (self.data.nodes[0].visited):
            for i in range(self.data.node_count):
                self.data.nodes[i].visited = 0
        cdef char *kmer = <char *> malloc(sizeof(char) * (k * 2 + 1))
        memset(kmer, 0, k * 2 + 1)
        self.get_kmers(kmer, 0, 0, 0, k, index)
        free(kmer)
        return index

    @staticmethod
    def generate():
        return None

    def print_graph(self):
        cdef Node_t *node = self.data.nodes
        cdef int edge_count = node.edge_count
        cdef int i
        printf("Graph:\n0:%s\n", node.sequence)
        while edge_count > 0:
            for i in range(edge_count):
                printf("%u:%s ", node.edges[i], (self.data.nodes + node.edges[i]).sequence)
            printf("\n")
            node = self.data.nodes + node.edges[0]
            edge_count = node.edge_count
