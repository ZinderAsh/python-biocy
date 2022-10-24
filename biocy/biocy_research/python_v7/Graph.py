import numpy as np
from graph_kmer_index import kmer_hash_to_sequence

class Node:
    def __init__(self, sequence, edges):
        self.sequence_len = len(sequence)
        self.num_segments = 0 if self.sequence_len == 0 else 1 + self.sequence_len // 31
        self.sequence = np.empty((self.num_segments,), dtype=np.ulonglong)
        for i in range(self.num_segments):
            segment_end = min((i + 1) * 31, self.sequence_len)
            self.sequence[i] = hash_max_kmer(sequence.encode('ASCII')[(i * 31):segment_end], segment_end - (i * 31))
        #print(kmer_hash_to_sequence(np.right_shift(self.sequence[0], np.uint64(64 - self.sequence_len * 2)), self.sequence_len))
        self.edges = np.array(edges, dtype=np.int32)
        self.reference = False

# This one is specifically for 0 < k < 16. There will be a separate class for 15 < k < 32 because it needs two buffer int64 variables.
class GraphKmerFinder:

    def __init__(self, graph, k, max_variant_nodes):
        self.graph = graph
        self.k = k
        self.max_variant_nodes = max_variant_nodes
        self.kmers = None
        self.nodes = None
        self.kmer_mask = np.uint64((1 << (k * 2)) - 1)
        self.kmer_buffer = 0
        self.path_buffer = np.empty(k * 4, dtype=np.uint32)

    def find(self):
        self.kmers = []
        self.nodes = []
        for node_id in range(len(self.graph.nodes)):
            print(node_id)
            if self.graph.nodes[node_id].sequence_len != 0:
                self.get_kmers(self.k, node_id)
        return self.kmers, self.nodes

    def get_kmers_recursive(self, k, node_id, variant_cnt, path_len, start_kmer_len, kmer_len):
        # Update lenghts and buffers
        self.path_buffer[path_len] = node_id
        path_len += 1
        node = self.graph.nodes[node_id]
 
        # If the node is a variant node, ensure max_variant_nodes is adhered to.       
        if not node.reference:
            variant_cnt += 1
            if variant_cnt > self.max_variant_nodes:
                return

        if node.num_segments == 1:
            # Clear and replace the end of the buffer with the node currently being visited.
            self.kmer_buffer = np.bitwise_and(self.kmer_buffer, np.left_shift(np.uint64(-1), np.uint64(64 - (kmer_len * 2))))
            self.kmer_buffer = np.bitwise_or(self.kmer_buffer, np.right_shift(node.sequence[0], np.uint64(kmer_len * 2)))
            new_kmer_len = min(kmer_len + node.sequence_len, start_kmer_len + k - 1)
            # Quickly make sure kmer_len catches up to avoid needless iterations in the following loop
            if kmer_len < k - 1:
                kmer_len = min(k - 1, new_kmer_len)
            # Iterate through and index the kmers
            while kmer_len < new_kmer_len:
                kmer_len += 1
                if kmer_len >= k:
                    # Index the kmer for all nodes in the path
                    for i in range(path_len):
                        self.nodes.append(self.path_buffer[i])
                        self.kmers.append(np.bitwise_and(np.right_shift(self.kmer_buffer, np.uint64(64 - (kmer_len * 2))), self.kmer_mask))
        else:
            pass
        
        # Stop if the node that first called the recursive function no longer has any of their bases included in results.
        if kmer_len >= start_kmer_len + k - 1:
            return

        # Continue recursion
        for edge in node.edges:
            self.get_kmers_recursive(k, edge, variant_cnt, path_len, start_kmer_len, kmer_len)

    def get_kmers(self, k, node_id):
        # Update lengths and buffers
        self.path_buffer[0] = node_id
        node = self.graph.nodes[node_id]

        # If the node is a variant node, ensure max_variant_nodes is adhered to.
        variant_cnt = 0 if node.reference else 1
        if variant_cnt > self.max_variant_nodes:
            return

        # Store the node's sequence in the buffer.
        self.kmer_buffer = node.sequence[0]
        kmer_len = min(k, node.sequence_len)
        # If at least k bases are stored, iterate through and index the kmers.
        if kmer_len == k:
            kmer_len -= 1
            while kmer_len < node.sequence_len:
                kmer_len += 1
                self.nodes.append(node_id)
                self.kmers.append(np.bitwise_and(np.right_shift(self.kmer_buffer, np.uint64(64 - (kmer_len * 2))), self.kmer_mask))
        if node.num_segments > 1:
            pass # TODO: Handle nodes with more than 31 bases.
        
        # Keep only the last k - 1 bases in the buffer.
        if kmer_len >= k:
            new_kmer_len = min(kmer_len, k - 1)
            self.kmer_buffer = np.left_shift(self.kmer_buffer, np.uint64((kmer_len - new_kmer_len) * 2))
            kmer_len = new_kmer_len

        # Recurse to all edges, remembering the kmer_len recursion started at.
        for edge in node.edges:
            self.get_kmers_recursive(k, edge, variant_cnt, 1, kmer_len, kmer_len)

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

def hash_max_kmer(arr, k):
    hashed = 0
    for i in range(k):
        hashed |= (arr[i] & 6) << (61 - i * 2)
    return hashed
