class Node:
    def __init__(self, sequence, edges):
        self.sequence = sequence
        self.edges = edges

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

    def get_kmers(self, k, index, start_node_id, node_id, kmer_buf="", kmer_len=0, rec_kmer_len=0, recurse=False):
        for base in self.nodes[node_id].sequence:
            kmer_buf += base
            kmer_len += 1
            if kmer_len >= k:
                kmer = kmer_buf[(kmer_len-k):kmer_len].encode('ASCII')
                if kmer not in index:
                    index[kmer] = []
                index[kmer].append(start_node_id)
            if recurse and kmer_len >= rec_kmer_len + k - 1:
                break
        if not recurse:
            kmer_buf = kmer_buf[max([kmer_len - k + 1, 0]):kmer_len]
            kmer_len = min([kmer_len, k-1])
        if not recurse or kmer_len < rec_kmer_len + k - 1:
            for edge in self.nodes[node_id].edges:
                rkl = rec_kmer_len if recurse else kmer_len
                self.get_kmers(k, index, start_node_id, edge, kmer_buf=kmer_buf, kmer_len=kmer_len, rec_kmer_len=rkl, recurse=True)


    def create_kmer_index(self, k):
        index = {}
        for node_id in range(len(self.nodes)):
            self.get_kmers(k, index, node_id, node_id)
        return index

