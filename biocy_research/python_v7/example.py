from Graph import Graph, Node, hash_kmer
from obgraph import Graph as OBGraph
from graph_kmer_index import kmer_hash_to_sequence
from graph_kmer_index.kmer_finder import DenseKmerFinder

def print_pairs(kmers, nodes, k):
    for i in range(len(kmers)):
        print(i, kmer_hash_to_sequence(kmers[i], k), kmers[i], nodes[i])

if False:
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

print("Creating graph from npz file...")
obg = OBGraph.from_file("data/example_graph.npz")

graph = Graph.from_obgraph(obg)

#print(graph.nodes[1].sequence_len)
#print(len(obg.sequences[1]))
#for i in range(graph.nodes[1].num_subsequences):
#    print(i, kmer_hash_to_sequence(graph.nodes[1].sequence[i], 32),
#          obg.sequences[1][(i*31):min((i+1)*31, graph.nodes[1].sequence_len)])

#exit()
k = 6
print(f"Indexing {k}mers with new solution...")
kmers, nodes = graph.create_kmer_index(k, max_variant_nodes=1000)
#count = 0
#kmer = hash_kmer(b"cacgtc", 6)
#for i in kmers:
#    if i == kmer:
#        count += 1
#print(count)
#exit(0)
print(f"Indexing {k}mers with kage solution...")
finder = DenseKmerFinder(obg, k=k, max_variant_nodes=1000)
finder.find()
ob_kmers, ob_nodes = finder.get_found_kmers_and_nodes()

#i = 0
#while nodes[i] == 1 or ob_nodes[i] == 1:
#    i += 1
#    print(i, kmer_hash_to_sequence(kmers[i], k), kmer_hash_to_sequence(ob_kmers[i], k))
