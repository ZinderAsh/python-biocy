from biocy import Graph, KmerFinder, hash_kmer
from obgraph import Graph as OBGraph
from graph_kmer_index import kmer_hash_to_sequence
from graph_kmer_index.kmer_finder import DenseKmerFinder
from timeit import timeit
import numpy as np
import npstructures as nps

def print_pairs(kmers, nodes):
    for i in range(len(kmers)):
        print(kmers[i], nodes[i])

#obgraph = OBGraph.from_file("../data/obgraph_chr1.npz")
#print(len(obgraph.nodes))

#exit(0)

print("Loading OBGraph")
obgraph = OBGraph.from_file("../data/obgraph_chr1.npz")
print("Creating Biocy Graph")
graph = Graph.from_obgraph(obgraph)

print("Reading variant to nodes file")
variant_to_nodes = np.load("../data/variant_to_nodes_chr1.npz")
ref_nodes = variant_to_nodes['ref_nodes']
var_nodes = variant_to_nodes['var_nodes']

print("Initializing kmer finder")
kmer_finder = KmerFinder(graph, 31)

#hash_map = nps.Counter([4, 10, 8, 40, 2207261504872468052])
#hash_map.count([4, 4, 8, 8, 8, 10, 2207261504872468052, 2207261504872468052])
#hash_map.count([11, 11])
#kmer_finder.set_kmer_frequency_index(hash_map)

print(kmer_finder.find_kmers_spanning_node(ref_nodes[300], max_variant_nodes=2))

exit(0)

kmer_finder.create_frequency_index()

print("Finding identifying windows for", len(ref_nodes), "variants")
ref, var = kmer_finder.find_variant_signatures(ref_nodes, var_nodes, align_windows=False, minimize_overlaps=False)

exit(0)

ref_nodes = variant_to_nodes['ref_nodes']
var_nodes = variant_to_nodes['var_nodes']
print("Finding identifying windows for", len(ref_nodes), "variants")
ref, var = kmer_finder.find_variant_signatures(ref_nodes, var_nodes, align_windows=False, minimize_overlaps=False)
ref2, var2 = kmer_finder.find_variant_signatures(ref_nodes, var_nodes, align_windows=True, minimize_overlaps=False)

print("Found windows for", len(ref), "variants")

for i in range(len(ref)):
    if len(ref2[i]) > 100 or len(var2[i]) > 100:
        print(i, len(ref2[i]), len(var2[i]))

diff = 0
for i in range(len(ref)):
    if ref[i][0] != ref2[i][0] or var[i][0] != var2[i][0]:
        diff += 1
        #print(i)
        #print("----------------------------")
        #print("ref1", ref[i])
        #print("var1", var[i])
        #print("ref2", ref2[i])
        #print("var2", var2[i])

print(diff)

#for i in range(10):
#    print("----------------------------")
#    print("ref1", ref[i])
#    print("var1", var[i])
#    print("ref2", ref2[i])
#    print("var2", var2[i])


variant_to_nodes.close()

#graph = Graph.from_sequence_edge_lists(
#    ["ACTGA", "G", "A", "C", "G", "T", "G", "A", "C", "G", "T", "ACTGATACTACTAC"],
#    [[1, 2], [3, 4], [3, 4], [5, 6], [5, 6], [7, 8], [7, 8], [9, 10], [9, 10], [11], [11], []],
#    ref=[0, 1, 3, 5, 7, 9, 11]
#)

exit(0)

k = 5
kmer_finder = KmerFinder(graph, 5)
ref, var = kmer_finder.find_variant_signatures([1, 3, 5, 7, 9], [2, 4, 6, 8, 10], reverse_kmers=True)

print("Ref")
for i in ref:
    kmer_str = ""
    for j in i:
        kmer_str += " " + kmer_hash_to_sequence(j, 3)
    print(kmer_str)
print("")
print("Var")
for i in var:
    kmer_str = ""
    for j in i:
        kmer_str += " " + kmer_hash_to_sequence(j, 3)
    print(kmer_str)

exit(0)

obgraph = OBGraph.from_file("tests/data/example_graph.npz")
graph = Graph.from_obgraph(obgraph)
obgraph2 = graph.to_obgraph()

k = 12

finder = DenseKmerFinder(obgraph, k=k, max_variant_nodes=k-2)
finder.find()
kmers, nodes = finder.get_found_kmers_and_nodes()

finder2 = DenseKmerFinder(obgraph2, k=k, max_variant_nodes=k-2)
finder2.find()
ob_kmers, ob_nodes = finder.get_found_kmers_and_nodes()

print(len(kmers))
print(len(ob_kmers))

counts = []
ob_counts = []
kmers_to_count = []
for i in range(1, 5):
    index = round(len(kmers) * i / 5)
    counts.append(0)
    ob_counts.append(0)
    kmers_to_count.append(kmers[index])

for j, kmer in enumerate(kmers_to_count):
    print("Counting", kmer)
    for i in range(len(kmers)):
        if i % 1000000 == 0:
            print(i, "/", len(kmers))
        if kmers[i] == kmer:
            counts[j] += 1
        if ob_kmers[i] == kmer:
            ob_counts[j] += 1

for i, kmer in enumerate(kmers_to_count):
    print(kmer, counts[i], ob_counts[i])

exit(0)


k = 31

def time_dense_kmer_finder():
    finder = DenseKmerFinder(obgraph, k=k, max_variant_nodes=k)
    finder.find()

def time_new_kmer_finder():
    #graph = Graph.from_gfa("../../data/chr21.gfa", compress=True)
    graph = Graph.from_file("../../data/chr21.bcg")
    #kmers, nodes = graph.create_kmer_index(k, max_variant_nodes=k)
    graph.to_file("../../data/chr21.bcg")
    print(len(nodes))

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
