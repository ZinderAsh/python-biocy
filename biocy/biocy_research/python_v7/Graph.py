import numpy as np
from graph_kmer_index import kmer_hash_to_sequence

class Node:
    def __init__(self):
        self.reference = False

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


class GraphKmerFinder:

    def __init__(self, graph, k, max_variant_nodes):
        self.graph = graph
        self.k = k
        self.max_variant_nodes = max_variant_nodes
        self.kmers = None
        self.nodes = None
        self.kmer_mask = np.uint64((1 << (k * 2)) - 1)
        self.shift = np.uint64((32 - k) * 2)
        self.kmer_buffer = 0 # Main kmer buffer
        self.kmer_buffer_ext = 0 # Extended kmer buffer
        self.path_buffer = np.empty(k * 4, dtype=np.uint32)

    def find(self):
        self.kmers = []
        self.nodes = []
        for node_id in range(len(self.graph.nodes)):
            if self.graph.nodes[node_id].sequence_len != 0:
                self.get_kmers(self.k, node_id)
        return self.kmers, self.nodes

    def get_kmers_recursive(self, k, node_id, variant_cnt, path_len, kmer_len, kmer_ext_len):
        # Update lenghts and buffers
        self.path_buffer[path_len] = node_id
        path_len += 1
        node = self.graph.nodes[node_id]
 
        # If the node is a variant node, ensure max_variant_nodes is adhered to.       
        if not node.reference:
            variant_cnt += 1
            if variant_cnt > self.max_variant_nodes:
                return

        if node.num_subsequences != 0:
            if kmer_ext_len == 0:
                # First recursion, replace entire buffer.
                self.kmer_buffer_ext = node.sequence[0]
            else:
                # Clear and replace the end of the buffer with the node currently being visited.
                self.kmer_buffer_ext = np.bitwise_and(self.kmer_buffer_ext,
                                                      np.left_shift(np.uint64(-1),
                                                                    np.uint64(64 - (kmer_ext_len * 2))))
                self.kmer_buffer_ext = np.bitwise_or(self.kmer_buffer_ext,
                                                     np.right_shift(node.sequence[0],
                                                                    np.uint64(kmer_ext_len * 2)))
            sequence_len = min(k - 1, kmer_ext_len + node.sequence_len)
            # Check if there are enough bases to index new kmers
            if kmer_len + sequence_len >= k:
                if kmer_ext_len < k - kmer_len - 1:
                    kmer_ext_len = k - kmer_len - 1
                while kmer_ext_len < sequence_len:
                    kmer_ext_len += 1
                    # Combine the main buffer and extended buffer to form the final kmer.
                    kmer_hash = np.bitwise_or(np.left_shift(self.kmer_buffer,
                                                            np.uint64(kmer_ext_len * 2)),
                                              np.right_shift(self.kmer_buffer_ext,
                                                             np.uint64(64 - kmer_ext_len * 2)))
                    kmer_hash = np.bitwise_and(kmer_hash, self.kmer_mask)
                    # Add the kmer to the index for every node in the path.
                    for i in range(path_len):
                        self.nodes.append(self.path_buffer[i])
                        self.kmers.append(kmer_hash)
            else:
                kmer_ext_len = sequence_len
                
            # Stop if the node that first called the recursive function no longer has any of their bases included in results.
            if kmer_ext_len >= k - 1:
                return

        # Continue recursion
        for edge in node.edges:
            self.get_kmers_recursive(k, edge, variant_cnt, path_len, kmer_len, kmer_ext_len)

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
        sequence_len = node.sequence_len
        kmer_len = min(k, sequence_len)
        subsequence_idx = 1
        subsequence_pos = np.uint64(0)
        # If at least k bases are stored, iterate through and index the kmers.
        if kmer_len == k:
            kmer_len -= 1
            while kmer_len < sequence_len:
                kmer_len += 1
                self.nodes.append(node_id)
                self.kmers.append(np.bitwise_and(np.right_shift(self.kmer_buffer,
                                                                np.uint64(64 - (kmer_len * 2))),
                                                 self.kmer_mask))
                # Check if the buffer is full, and if so, shift values along as much as k allows
                if kmer_len == 31 and subsequence_idx < node.num_subsequences:
                    self.kmer_buffer = np.left_shift(self.kmer_buffer, self.shift)
                    subsequence = np.left_shift(node.sequence[subsequence_idx], subsequence_pos)
                    self.kmer_buffer = np.bitwise_or(self.kmer_buffer,
                                                     np.right_shift(subsequence,
                                                                    np.uint64(62 - self.shift)))
                    subsequence_pos = np.uint64(subsequence_pos + self.shift)
                    if subsequence_pos > 62:
                        subsequence_idx += 1
                        if subsequence_idx < node.num_subsequences:
                            subsequence_pos = np.uint64(subsequence_pos - 62)
                            self.kmer_buffer = np.bitwise_or(self.kmer_buffer,
                                                             np.right_shift(node.sequence[subsequence_idx],
                                                                            np.uint64(62 - subsequence_pos)))
                    elif subsequence_pos == 62:
                        subsequence_idx += 1
                        subsequence_pos = np.uint64(0)
                    kmer_len = k - 1
                    sequence_len -= 32 - k

        
        # Right-align the buffer for easier bit-wise operations with kmer_buffer_ext in recursion.
        self.kmer_buffer = np.right_shift(self.kmer_buffer, 
                                          np.uint64(64 - kmer_len * 2))
        # Fake the kmer being max k - 1 as no more bases than that matter during recursion.
        kmer_len = min(kmer_len, k - 1)

        # Recurse to all edges, remembering the kmer_len recursion started at.
        for edge in node.edges:
            self.get_kmers_recursive(k, edge, variant_cnt, 1, kmer_len, 0)

    def print_buffers(self, msg=None):
        if msg:
            print(msg)
        print(kmer_hash_to_sequence(self.kmer_buffer, 32),
              kmer_hash_to_sequence(self.kmer_buffer_ext, 32))

class Graph:

    def __init__(self):
        self.nodes = []

    @staticmethod
    def from_obgraph(obg):
        graph = Graph()
        for i in range(len(obg.nodes)):
            node = Node.from_obgraph(obg.nodes[i], obg.sequences[i], obg.edges[i])
            graph.nodes.append(node)
        ref = obg.linear_ref_nodes()
        for i in ref:
            graph.nodes[i].reference = True

        return graph

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
            node = Node.from_sequence_edges(sequences[i], edges[i])
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

def pack_max_kmer(arr, start, end):
    packed = 0
    for i in range(start, end):
        packed |= (arr[i] << (62 - (i - start) * 2))
    return packed
