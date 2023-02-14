from obgraph import Graph as OBGraph
from graph_kmer_index.kmer_finder import DenseKmerFinder
from graph_kmer_index import kmer_hash_to_sequence
from biocy import Graph, KmerFinder, hash_kmer
from os.path import exists
import pytest

def hash_all(arr):
    if len(arr) == 0:
        return []
    k = len(arr[0])
    return [hash_kmer(i, k) for i in arr]

def compare_kmer_node_lists(kmers_a, nodes_a, kmers_b, nodes_b, k, print_lists=False):
    assert len(kmers_a) == len(kmers_b)
    assert len(nodes_a) == len(nodes_b)
    assert len(kmers_a) == len(nodes_a)
    counts_a = {}
    for i in range(len(kmers_a)):
        if kmers_a[i] not in counts_a:
            counts_a[kmers_a[i]] = 0
        counts_a[kmers_a[i]] += 1
    counts_b = {}
    for i in range(len(kmers_b)):
        if kmers_b[i] not in counts_b:
            counts_b[kmers_b[i]] = 0
        counts_b[kmers_b[i]] += 1
    all_keys = set(counts_a.keys())
    all_keys.update(counts_b.keys())
    for i in all_keys:
        a = counts_a[i] if i in counts_a else 0
        b = counts_b[i] if i in counts_b else 0
        if a != b and a > 0:
            print(kmer_hash_to_sequence(i, k), a, b)
    # Print entire kmer lists
    if print_lists:
        for i in range(max(len(kmers_a), len(kmers_b))):
            if i < len(kmers_a) and i < len(kmers_b):
                print(f'{kmer_hash_to_sequence(kmers_a[i], k)}, {nodes_a[i]}\t{kmer_hash_to_sequence(kmers_b[i], k)}, {nodes_b[i]}')
            elif i < len(kmers_a):
                print(f'{kmer_hash_to_sequence(kmers_a[i], k)}, {nodes_a[i]}')
            elif i < len(kmers_b):
                print(f'{" " * k},  \t{kmer_hash_to_sequence(kmers_b[i], k)}, {nodes_b[i]}')
    for i in all_keys:
        assert counts_a[i] == counts_b[i]

# The last tests here are quite large. The graph they share looks like this:
#         A   G   G   C   T
# ACTGA - G - C - T - A - G - ACTGA
#
@pytest.mark.parametrize("nodes,edges,ref,k,max_var", [ 
            (["ACTGACTGACTG", "ACTAGTC"], [[1], []], [0, 1], 3, 4),
            (["ACTGACTGACTG", "ACTGCA"], [[1], []], [0, 1], 4, 4),
            (["ACTAG", "GAC", "TGA", "CTGAGT"], [[1], [2], [3], []], [0, 1, 2, 3], 3, 4),
            (["ACTAG", "G", "A", "GTACTCA"], [[1, 2], [3], [3], []], [0, 1, 3], 3, 4),
            (["ACTAGATTTAGGCTA" * 4, "G", "A", "GTACTCA" * 10], [[1, 2], [3], [3], []], [0, 1, 3], 3, 4),
            (["ACTAGATTTAGGCTA" * 4, "GTACTAA" * 12, "ATGACTA" * 12, "GTACTCA" * 10], [[1, 2], [3], [3], []], [0, 1, 3], 3, 4),
            (["AGTAGA", "G", "CT", "ACTA", "G", "A", "TCATA"], [[1, 2], [3], [3], [4, 5], [6], [6], []],
             [0, 1, 3, 4, 6], 3, 4),
            (["AGTAGA", "G", "CT", "ACTA", "G", "A", "TCATA"], [[1, 2], [3], [3], [4, 5], [6], [6], []],
             [0, 1, 3, 4, 6], 4, 4),
            (["AGTAGA", "G", "CT", "ACTA", "G", "A", "TCATA"], [[1, 2], [3], [3], [4, 5], [6], [6], []],
             [0, 1, 3, 4, 6], 5, 4),
            (["AGTAGAACTGACTTCAGGTACTTA" * 10] * 7, [[1, 2], [3], [3], [4, 5], [6], [6], []],
             [0, 1, 3, 4, 6], 5, 4),
            (["ACTGA", "G", "A", "C", "G", "T", "G", "A", "C", "G", "T", "ACTGA"],
             [[1, 2], [3, 4], [3, 4], [5, 6], [5, 6], [7, 8], [7, 8], [9, 10], [9, 10], [11], [11], []],
             [0, 1, 3, 5, 7, 9, 11], 4, 250),
            (["ACTGA", "G", "A", "C", "G", "T", "G", "A", "C", "G", "T", "ACTGA"],
             [[1, 2], [3, 4], [3, 4], [5, 6], [5, 6], [7, 8], [7, 8], [9, 10], [9, 10], [11], [11], []],
             [0, 1, 3, 5, 7, 9, 11], 6, 2),
            (["ACTGA", "G", "A", "C", "G", "T", "G", "A", "C", "G", "T", "ACTGA"],
             [[1, 2], [3, 4], [3, 4], [5, 6], [5, 6], [7, 8], [7, 8], [9, 10], [9, 10], [11], [11], []],
             [0, 1, 3, 5, 7, 9, 11], 6, 4),
            (["ACTGA", "G", "A", "C", "G", "T", "G", "A", "C", "G", "T", "ACTGA"],
             [[1, 2], [3, 4], [3, 4], [5, 6], [5, 6], [7, 8], [7, 8], [9, 10], [9, 10], [11], [11], []],
             [0, 1, 3, 5, 7, 9, 11], 6, 6)
            ])
def test_kmer_index_against_kage(nodes, edges, ref, k, max_var):
    graph = Graph.from_sequence_edge_lists(nodes, edges, ref=ref)
    kmer_finder = KmerFinder(graph, k, max_variant_nodes=max_var)
    res_kmers, res_nodes = kmer_finder.find(reverse_kmers=True)
    ob_node_sequences = {}
    ob_edges = {}
    ob_linear_ref_nodes = []
    for i in range(len(nodes)):
        ob_node_sequences[i + 1] = nodes[i]
        if len(edges[i]) > 0:
            ob_edges[i + 1] = []
            for edge in edges[i]:
                ob_edges[i + 1].append(edge + 1)
    ob_linear_ref_nodes = [i + 1 for i in ref]
    obgraph = OBGraph.from_dicts(
        node_sequences=ob_node_sequences,
        edges=ob_edges,
        linear_ref_nodes=ob_linear_ref_nodes
    )
    finder = DenseKmerFinder(obgraph, k=k, max_variant_nodes=max_var)
    finder.find()
    ob_kmers, ob_nodes = finder.get_found_kmers_and_nodes()
    compare_kmer_node_lists(res_kmers, res_nodes, ob_kmers, ob_nodes, k)

# The last test here is quite large. The graph looks like this (* means a node with an empty sequence):
#      *       A       *       C
#    /   \   /   \   /   \   /   \
# AC - G - * - * - G - T - * - A - GT
#
@pytest.mark.parametrize("nodes,edges,ref,k,max_var,expected_nodes,expected_kmers", [
            (["ACTG", "ACTA"], [[1], []], [0, 1], 3, 4,
             [0, 0, 0, 1, 0, 1, 1, 1],
             hash_all(['ACT', 'CTG', 'TGA', 'TGA', 'GAC', 'GAC', 'ACT', 'CTA'])),
            (["ACT", "", "A", "GCA"], [[1, 2], [3], [3], []], [0, 1], 3, 4,
             [0, 0, 2, 0, 2, 3, 2, 3, 3, 0, 1, 3, 0, 1, 3],
             hash_all(['ACT', 'CTA', 'CTA', 'TAG', 'TAG', 'TAG', 'AGC', 'AGC', 'GCA', 'CTG', 'CTG', 'CTG',
                       'TGC', 'TGC', 'TGC'])),
            (["AC", "G", "", "", "", "A", "G", "T", "", "", "A", "C", "GT"],
             [[1, 2], [3], [3], [4, 5], [6], [6], [7, 8], [9], [9], [10, 11], [12], [12], []],
             [0, 1, 3, 4, 6, 7, 9, 10, 12], 3, 10,
             [0, 1, 0, 2, 3, 5, 0, 2, 3, 4, 6,
              0, 1, 3, 5, 0, 1, 3, 4, 6,
              0, 2, 3, 5, 6, 0, 2, 3, 4, 6, 7,
              0, 2, 3, 4, 6, 8, 9, 10,
              0, 2, 3, 4, 6, 8, 9, 11,
              1, 3, 5, 6, 1, 3, 4, 6, 7,
              1, 3, 4, 6, 8, 9, 10,
              1, 3, 4, 6, 8, 9, 11,
              5, 6, 7, 5, 6, 8, 9, 10,
              5, 6, 8, 9, 11, 6, 7, 9, 10,
              6, 7, 9, 11, 6, 8, 9, 10, 12,
              6, 8, 9, 11, 12, 7, 9, 10, 12,
              7, 9, 11, 12, 10, 12, 11, 12],
             hash_all(['ACG', 'ACG', 'ACA', 'ACA', 'ACA', 'ACA', 'ACG', 'ACG', 'ACG', 'ACG', 'ACG',
                       'CGA', 'CGA', 'CGA', 'CGA', 'CGG', 'CGG', 'CGG', 'CGG', 'CGG',
                       'CAG', 'CAG', 'CAG', 'CAG', 'CAG', 'CGT', 'CGT', 'CGT', 'CGT', 'CGT', 'CGT',
                       'CGA', 'CGA', 'CGA', 'CGA', 'CGA', 'CGA', 'CGA', 'CGA',
                       'CGC', 'CGC', 'CGC', 'CGC', 'CGC', 'CGC', 'CGC', 'CGC',
                       'GAG', 'GAG', 'GAG', 'GAG', 'GGT', 'GGT', 'GGT', 'GGT', 'GGT',
                       'GGA', 'GGA', 'GGA', 'GGA', 'GGA', 'GGA', 'GGA',
                       'GGC', 'GGC', 'GGC', 'GGC', 'GGC', 'GGC', 'GGC',
                       'AGT', 'AGT', 'AGT', 'AGA', 'AGA', 'AGA', 'AGA', 'AGA',
                       'AGC', 'AGC', 'AGC', 'AGC', 'AGC', 'GTA', 'GTA', 'GTA', 'GTA',
                       'GTC', 'GTC', 'GTC', 'GTC', 'GAG', 'GAG', 'GAG', 'GAG', 'GAG',
                       'GCG', 'GCG', 'GCG', 'GCG', 'GCG', 'TAG', 'TAG', 'TAG', 'TAG',
                       'TCG', 'TCG', 'TCG', 'TCG', 'AGT', 'AGT', 'CGT', 'CGT']))
            ])
def test_kmer_empty_nodes(nodes, edges, ref, k, max_var, expected_nodes, expected_kmers):
    graph = Graph.from_sequence_edge_lists(nodes, edges, ref=ref)
    kmer_finder = KmerFinder(graph, k, max_variant_nodes=max_var)
    res_kmers, res_nodes = kmer_finder.find()
    compare_kmer_node_lists(res_kmers, res_nodes, expected_kmers, expected_nodes, k)

@pytest.mark.slow
@pytest.mark.parametrize("file,k,max_var", [
        ("tests/data/example_graph.npz", 4, 4),
        ("tests/data/example_graph.npz", 6, 6),
        ("tests/data/example_graph.npz", 8, 8),
        ("tests/data/example_graph.npz", 12, 12),
        ("tests/data/example_graph.npz", 16, 16),
        ("tests/data/example_graph.npz", 24, 24),
        ("tests/data/example_graph.npz", 31, 31)
    ])
def test_obgraph_against_kage_big_graph(file, k, max_var):
    obgraph = OBGraph.from_file(file)
    graph = Graph.from_obgraph(obgraph)
    kmer_finder = KmerFinder(graph, k, max_variant_nodes=max_var)
    res_kmers, res_nodes = kmer_finder.find()
    fname = f'tests/data/kage_results_{k}mer_{max_var}var.txt'
    ob_kmers, ob_nodes = [], []
    if exists(fname):
        f = open(fname, "r")
        lines = f.readlines()
        for l in lines:
            seg = l.split(" ")
            if len(seg) == 2:
                ob_kmers.append(int(seg[0]))
                ob_nodes.append(int(seg[1]))
        f.close()
    else:
        finder = DenseKmerFinder(obgraph, k=k, max_variant_nodes=max_var)
        finder.find()
        ob_kmers, ob_nodes = finder.get_found_kmers_and_nodes()
        f = open(fname, "w")
        for i in range(len(ob_kmers)):
            f.write(f'{str(ob_kmers[i])} {str(ob_nodes[i])}\n')
        f.close()
    compare_kmer_node_lists(res_kmers, res_nodes, ob_kmers, ob_nodes, k)

