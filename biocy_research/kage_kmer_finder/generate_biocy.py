from biocy import Graph
from obgraph import Graph as OBGraph
from graph_kmer_index.kmer_finder import DenseKmerFinder
from timeit import timeit
import numpy as np

#obgraph = OBGraph.from_file("example_graph.npz")
#graph = Graph.from_obgraph(obgraph)
#graph = Graph.from_file("../../data/biocy_graph.bcg")

def dense_kmer_finder():
    graph = Graph.from_gfa("../../data/vg_graph.gfa", compress=False)
    k = 31
    kmers, nodes = graph.create_kmer_index(k, max_variant_nodes=k, big_endian=False)
    kmers.tofile("biocy_kmers.npy")
    nodes.tofile("biocy_nodes.npy")
    
    graph = Graph.from_gfa("../../data/vg_graph.gfa", compress=True)
    k = 31
    kmers, nodes = graph.create_kmer_index(k, max_variant_nodes=k, big_endian=False)
    kmers.tofile("biocy_kmers2.npy")
    nodes.tofile("biocy_nodes2.npy")

dense_kmer_finder()

