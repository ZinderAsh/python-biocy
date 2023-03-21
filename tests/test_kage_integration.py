import sys
import obgraph
import numpy as np
from shared_memory_wrapper import from_file
from biocy import KmerFinder, Graph
from npstructures import Counter
from obgraph.variant_to_nodes import VariantToNodes


def small_test():
    g = obgraph.Graph.from_dicts({
            1: "ACTG",
            2: "A",
            3: "G",
            4: "CCCC"
        },
        edges={
            1: [2, 3],
            2: [4],
            3: [4]
        },
        linear_ref_nodes=[1, 2, 4],
        chromosome_start_nodes={"chr1": 1}
    )

    variant_to_nodes = obgraph.variant_to_nodes.VariantToNodes(np.array([2], dtype=np.uint32), np.array([3], dtype=np.uint32))
    counter = Counter([0, 1, 2], [1, 10, 1])
    n = 1000000
    keys = np.unique(np.random.randint(0, n*2, n, dtype=np.uint64))

    #counter = Counter(keys, np.random.randint(0, 5, len(keys)))
    counter = _get_counter("linear_kmers_counter.npz")
    #counter = from_file("testcounter5.npz")
    print(counter._keys.dtype, counter._values.dtype)
    print("N keys:", len(counter._keys))

    #graph = Graph.from_obgraph(g)
    graph = Graph.from_obgraph(obgraph.Graph.from_file("obgraph.npz"))

    kmer_finder = KmerFinder(graph, k=3, max_variant_nodes=4)
    kmer_finder.set_kmer_frequency_index(counter)
    results = kmer_finder.find_variant_signatures(variant_to_nodes.ref_nodes, variant_to_nodes.var_nodes)
    print(results)


def _get_counter(counter_fn):
    kmer_counter = from_file(counter_fn)
    counter = kmer_counter.counter
    # Counter is a HashTable, make a Counter instead
    counter = Counter(counter._keys.ravel(), counter._values.ravel())
    counter._values = counter._values.astype(np.int64)
    return counter


def main(obgraph_fn, kmer_counter_fn, variant_to_nodes_fn):
    k = 31
    variant_to_nodes = VariantToNodes.from_file(variant_to_nodes_fn)
    ref_nodes = variant_to_nodes.ref_nodes
    var_nodes = variant_to_nodes.var_nodes

    graph = Graph.from_obgraph(obgraph.Graph.from_file(obgraph_fn))
    counter = _get_counter(kmer_counter_fn)

    kmer_finder = KmerFinder(graph, k, max_variant_nodes=4)
    kmer_finder.set_kmer_frequency_index(counter)

    print("Finding signatures")
    results = kmer_finder.find_variant_signatures(ref_nodes, var_nodes)
    print(results)


if __name__ == "__main__":
    #small_test()
    main(*sys.argv[1:])
    # example:
    # python test_kage_integration.py obgraph.npz linear_kmers_counter.npz variant_to_nodes.npz
