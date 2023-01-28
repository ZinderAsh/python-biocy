from biocy import Graph
from obgraph import Graph as OBGraph
from graph_kmer_index.kmer_finder import DenseKmerFinder
from timeit import timeit
import numpy as np

obgraph = OBGraph.from_file("example_graph.npz")
k = 31

empty_nodes = 0

for i in obgraph.nodes[1:]:
    if i == 0:
        empty_nodes += 1

print("Total Nodes:", len(obgraph.nodes) - 1)
print("Empty Nodes:", empty_nodes)

print("nodes", obgraph.nodes.dtype, obgraph.nodes)
print("edges", obgraph.edges.dtype, obgraph.edges)
print("sequences", obgraph.sequences.dtype, obgraph.sequences)
print("node_to_ref_offset", obgraph.node_to_ref_offset.dtype, obgraph.node_to_ref_offset)
print("ref_offset_to_node", obgraph.ref_offset_to_node.dtype, obgraph.ref_offset_to_node)
print("_linear_ref_nodes_cache", obgraph._linear_ref_nodes_cache)
print("chromosome_start_nodes", "python list", obgraph.chromosome_start_nodes)
print("allele_frequencies", obgraph.allele_frequencies)
print("numeric_node_sequences", obgraph.numeric_node_sequences)
print("linear_ref_nodes_index", obgraph.linear_ref_nodes_index.dtype, obgraph.linear_ref_nodes_index)
print("linear_ref_nodes_and_dummy_nodes_index", obgraph.linear_ref_nodes_and_dummy_nodes_index.dtype, obgraph.linear_ref_nodes_and_dummy_nodes_index)

