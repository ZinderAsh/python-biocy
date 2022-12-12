from biocy import Graph
from obgraph import Graph as OBGraph
from graph_kmer_index.kmer_finder import DenseKmerFinder
from timeit import timeit
import numpy as np

obgraph = OBGraph.from_file("example_graph.npz")
k = 31

def dense_kmer_finder():
    finder = DenseKmerFinder(obgraph, k=k, max_variant_nodes=k)
    finder.find() 
    ob_kmers, ob_nodes = finder.get_found_kmers_and_nodes()
    ob_kmers.tofile("ob_kmers.npy")
    ob_nodes.tofile("ob_nodes.npy")

dense_kmer_finder()

