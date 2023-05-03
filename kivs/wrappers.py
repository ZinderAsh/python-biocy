from npstructures import Counter
from . import KmerFinder, Graph
import numpy as np


def get_variant_signatures(obgraph, variant_to_nodes, counter, k=31):
    graph = Graph.from_obgraph(obgraph)

    ref_nodes = variant_to_nodes.ref_nodes
    var_nodes = variant_to_nodes.var_nodes

    # Counter is a HashTable, make a Counter instead
    counter = Counter(counter._keys.ravel().copy(), counter._values.ravel().copy())
    counter._values = counter._values.astype(np.int64)

    kmer_finder = KmerFinder(graph, k, reverse_kmers=True)
    kmer_finder.set_frequency_index(counter)

    results = kmer_finder.find_variant_signatures(
            ref_nodes, var_nodes, max_variant_nodes=200, align_windows=True, minimize_overlaps=True)

    return results
