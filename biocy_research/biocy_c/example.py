from biocy import Graph
from obgraph import Graph as OBGraph
from graph_kmer_index.kmer_finder import DenseKmerFinder
from timeit import timeit

def print_pairs(kmers, nodes):
    for i in range(len(kmers)):
        print(kmers[i], nodes[i])

obgraph = OBGraph.from_file("../../data/chr21_v2.npz")
graph = Graph.from_obgraph(obgraph)
#graph = Graph.from_file("../tests/data/chr21.bcg")

k = 31

def time_dense_kmer_finder():
    finder = DenseKmerFinder(obgraph, k=k, max_variant_nodes=k)
    finder.find()

def time_new_kmer_finder():
    graph = Graph.from_obgraph(obgraph)
    graph.create_kmer_index(k, max_variant_nodes=k)

#graph.to_file("../tests/data/chr21.bcg")

#time_dense_kmer_finder()
#print(len(obgraph.nodes))
#time_new_kmer_finder()

#kmers, nodes = graph.create_kmer_index(k, max_variant_nodes=k)
#print(kmers[:10])
#print(kmers[-10:])

time_new_kmer_finder()

print("Finding kmers with kage")
#ret = 0;
#ret = timeit(lambda: time_new_kmer_finder(), number=1)
#print("Done in", ret)

exit(0)

nodes = ["AGTA", "G", "CT", "A", "CTA", "G", "A", "T"]
edges = [[1, 2], [3], [3], [4], [5, 6], [7], [7], []]
reference = [0, 1, 3, 4, 5, 7]
graph = Graph.from_sequence_edge_lists(nodes, edges, ref=reference)
print("Finding 3-mers")
kmers, nodes = graph.create_kmer_index(3)
print_pairs(kmers, nodes)
print("Finding 4-mers")
kmers, nodes = graph.create_kmer_index(4)
print_pairs(kmers, nodes)
print("Finding 5-mers")
kmers, nodes = graph.create_kmer_index(5)
print_pairs(kmers, nodes)
