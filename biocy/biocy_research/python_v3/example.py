from Graph import Graph

def print_pairs(kmers, nodes):
    for i in range(len(kmers)):
        print(kmers[i], nodes[i])

nodes = ["AGTA", "G", "CT", "A", "CTA", "G", "A", "T"]
edges = [[1, 2], [3], [3], [4], [5, 6], [7], [7], []]
graph = Graph.from_sequence_edge_lists(nodes, edges)
print("Finding 3-mers")
kmers, nodes = graph.create_kmer_index(3)
print_pairs(kmers, nodes)
print("Finding 4-mers")
kmers, nodes = graph.create_kmer_index(4)
print_pairs(kmers, nodes)
print("Finding 5-mers")
kmers, nodes = graph.create_kmer_index(5)
print_pairs(kmers, nodes)
