import numpy as np

class Node:
    def __init__(self, sequence, edges):
        self.sequence = np.frombuffer(sequence.encode('ASCII'), dtype=np.dtype('B'))
        self.edges = np.array(edges, dtype=np.int32)

class Graph:

    def __init__(self):
        self.nodes = []

    @staticmethod
    def from_sequence_edge_lists(sequences, edges):
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
        return graph

    def get_kmers(self, k, kmers, nodes, node_id, kmer_buf, nodes_buf,
                  kmer_len=0, nodes_len=0, rec_kmer_len=0, recurse=False):
        nodes_buf[nodes_len] = node_id
        nodes_len += 1
        for base in self.nodes[node_id].sequence:
            kmer_buf[kmer_len] = base
            if kmer_len == k * 2 - 1: # kmer buffer full, left shift
                for i in range(kmer_len):
                    kmer_buf[i] = kmer_buf[i + 1]
            else:
                kmer_len += 1
            if kmer_len >= k: # full kmer collected, add to result
                kmer = hash_kmer(kmer_buf[(kmer_len-k):kmer_len], k)
                for i in range(nodes_len):
                    kmers.append(kmer)
                    nodes.append(nodes_buf[i])
            if recurse and kmer_len >= rec_kmer_len + k - 1:
                return # stop recursion
        if not recurse: # shorten kmer buffer to last k-1 bases before recursion
            new_kmer_len = min(kmer_len, k-1)
            for i in range(new_kmer_len):
                kmer_buf[i] = kmer_buf[kmer_len - new_kmer_len + i]
            kmer_len = new_kmer_len
        for edge in self.nodes[node_id].edges: # recurse
            rkl = rec_kmer_len if recurse else kmer_len
            self.get_kmers(k, kmers, nodes, edge, kmer_buf, nodes_buf,
                           kmer_len=kmer_len, nodes_len=nodes_len, rec_kmer_len=rkl, recurse=True)


    def create_kmer_index(self, k):
        kmers = []
        nodes = []
        kmer_buf = np.empty(k * 2, dtype=np.dtype('B'))
        nodes_buf = np.empty(k, dtype=np.uint32)
        for node_id in range(len(self.nodes)):
            self.get_kmers(k, kmers, nodes, node_id, kmer_buf, nodes_buf)
        return kmers, nodes


def hash_kmer(arr, k):
    hashed = 0
    for i in range(k):
        hashed |= (arr[i] & 6) << ((k - i - 1) << 1)
    return hashed >> 1
