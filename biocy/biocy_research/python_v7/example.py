from Graph import Graph
from graph_kmer_index import kmer_hash_to_sequence

def print_pairs(kmers, nodes, k):
    for i in range(len(kmers)):
        print(i, kmer_hash_to_sequence(kmers[i], k), kmers[i], nodes[i])

nodes = ["AGTA", "G", "CT", "A", "CTA", "G", "A", "T"]
edges = [[1, 2], [3], [3], [4], [5, 6], [7], [7], []]
reference = [0, 1, 3, 4, 5, 7]
graph = Graph.from_sequence_edge_lists(nodes, edges, ref=reference)
print("Finding 3-mers")
kmers, nodes = graph.create_kmer_index(3)
print_pairs(kmers, nodes, 3)
print("Finding 4-mers")
kmers, nodes = graph.create_kmer_index(4)
print_pairs(kmers, nodes, 4)
print("Finding 5-mers")
kmers, nodes = graph.create_kmer_index(5)
print_pairs(kmers, nodes, 5)
