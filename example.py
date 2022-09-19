from biocy import Graph

nodes = ["AGTA", "G", "CT", "A", "CTA", "G", "A", "T"]
edges = [[1, 2], [3], [3], [4], [5, 6], [7], [7], []]
graph = Graph.from_sequence_edge_lists(nodes, edges)
graph.print_graph()
print("Finding 3-mers")
index = graph.create_kmer_index(3)
print(index)
print("Finding 4-mers")
index = graph.create_kmer_index(4)
print(index)
print("Finding 5-mers")
index = graph.create_kmer_index(5)
print(index)
