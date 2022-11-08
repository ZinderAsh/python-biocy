import numpy as np

class Node:
    def __init__(self, sequence, edges):
        self.sequence = np.frombuffer(sequence.encode('ASCII'), dtype=np.dtype('B'))
        self.edges = np.array(edges, dtype=np.int32)
        self.reference = False

class GraphKmerFinder:

    def __init__(self, graph, k, max_variant_nodes):
        self.graph = graph
        self.k = k
        self.max_variant_nodes = max_variant_nodes
        self.kmers = None
        self.nodes = None
        self.kmer_buffer = np.empty(k * 2, dtype=np.dtype('B'))
        self.path_buffer = np.empty(k * 4, dtype=np.uint32)

    def find(self):
        self.kmers = []
        self.nodes = []
        for node_id in range(len(self.graph.nodes)):
            self.get_kmers(self.k, node_id)
        return self.kmers, self.nodes

    def get_kmers(self, k, node_id, variant_cnt=0, path_len=0, kmer_len=0, rec_kmer_len=0):
        # Update lengths and buffers
        self.path_buffer[path_len] = node_id
        path_len += 1
        node = self.graph.nodes[node_id]
        recursive = (path_len > 1)

        # If the node is a variant node, ensure max_variant_nodes is adhered to.
        if not node.reference:
            variant_cnt += 1
            if variant_cnt > self.max_variant_nodes:
                return

        # Add all bases to kmer_buffer, shifting and adding to result as necessary.
        for base in node.sequence:
            self.kmer_buffer[kmer_len] = base
            if kmer_len == self.k * 2 - 1:
                # kmer buffer full, left shift
                self.kmer_buffer[0:kmer_len] = self.kmer_buffer[1:kmer_len+1]
            else:
                kmer_len += 1
            if kmer_len >= self.k:
                # full kmer collected, add to result
                kmer = hash_kmer(self.kmer_buffer[(kmer_len-self.k):kmer_len], self.k)
                for i in range(path_len):
                    self.kmers.append(kmer)
                    self.nodes.append(self.path_buffer[i])
            # if recursive, stop once no more bases from starting node are involved.
            if recursive and kmer_len >= rec_kmer_len + self.k - 1:
                return

        # shorten kmer buffer to last k-1 bases before the first recursion.
        if not recursive:
            new_kmer_len = min(kmer_len, k-1)
            self.kmer_buffer[0:new_kmer_len] = self.kmer_buffer[(kmer_len-new_kmer_len):kmer_len]
            kmer_len = new_kmer_len

        # Recurse to all edges, remembering the kmer_len recursion started at.
        for edge in node.edges:
            rkl = rec_kmer_len if recursive else kmer_len
            self.get_kmers(k, edge, variant_cnt=variant_cnt, path_len=path_len,
                           kmer_len=kmer_len, rec_kmer_len=rkl)

class Graph:

    def __init__(self):
        self.nodes = []

    @staticmethod
    def from_sequence_edge_lists(sequences, edges, ref=None):
        """
        Args:
            sequences: ["ACT", "G", "A", "GT"]
            edges: [[1, 2], [3], [3], []]
        List index determines a node's ID, and edges refer to what IDs are a node is connected to.
        Assumes the first node is the start node
        """
        graph = Graph()
        for i in range(len(sequences)):
            node = Node(sequences[i], edges[i])
            graph.nodes.append(node)
        if ref is not None:
            for i in ref:
                graph.nodes[i].reference = True
        else:
            for node in graph.nodes:
                node.reference = True
        return graph


    def create_kmer_index(self, k, max_variant_nodes=4):
        kmers, nodes = GraphKmerFinder(self, k, max_variant_nodes).find()
        return kmers, nodes


def hash_kmer(arr, k):
    hashed = 0
    for i in range(k):
        hashed |= (arr[i] & 6) << ((k - i - 1) << 1)
    return hashed >> 1
