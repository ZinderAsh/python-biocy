from obgraph import Graph as OBGraph
from graph_kmer_index.kmer_finder import DenseKmerFinder
from graph_kmer_index import kmer_hash_to_sequence
from Graph import Graph
import pytest

@pytest.mark.parametrize("nodes,edges,ref,k,max_var", [
            (["ACTGACTGACTG", "ACTAGTC"], [[1], []], [0, 1], 3, 4),
            (["ACTGACTGACTG", "ACTGCA"], [[1], []], [0, 1], 4, 4),
            (["ACTAG", "GAC", "TGA", "CTGAGT"], [[1], [2], [3], []], [0, 1, 2, 3], 3, 4),
            (["ACTAG", "G", "A", "GTACTCA"], [[1, 2], [3], [3], []], [0, 1, 3], 3, 4),
            (["AGTAGA", "G", "CT", "ACTA", "G", "A", "TCATA"], [[1, 2], [3], [3], [4, 5], [6], [6], []],
             [0, 1, 3, 4, 6], 3, 4),
            (["AGTAGA", "G", "CT", "ACTA", "G", "A", "TCATA"], [[1, 2], [3], [3], [4, 5], [6], [6], []],
             [0, 1, 3, 4, 6], 4, 4),
            (["AGTAGA", "G", "CT", "ACTA", "G", "A", "TCATA"], [[1, 2], [3], [3], [4, 5], [6], [6], []],
             [0, 1, 3, 4, 6], 5, 4),
            (["ACTGA", "G", "A", "C", "G", "T", "G", "A", "C", "G", "T", "ACTGA"],
             [[1, 2], [3, 4], [3, 4], [5, 6], [5, 6], [7, 8], [7, 8], [9, 10], [9, 10], [11], [11], []],
             [0, 1, 3, 5, 7, 9, 11], 4, 1000),
            (["ACTGA", "G", "A", "C", "G", "T", "G", "A", "C", "G", "T", "ACTGA"],
             [[1, 2], [3, 4], [3, 4], [5, 6], [5, 6], [7, 8], [7, 8], [9, 10], [9, 10], [11], [11], []],
             [0, 1, 3, 5, 7, 9, 11], 4, 2),
            (["ACTGA", "G", "A", "C", "G", "T", "G", "A", "C", "G", "T", "ACTGA"],
             [[1, 2], [3, 4], [3, 4], [5, 6], [5, 6], [7, 8], [7, 8], [9, 10], [9, 10], [11], [11], []],
             [0, 1, 3, 5, 7, 9, 11], 4, 4),
            (["ACTGA", "G", "A", "C", "G", "T", "G", "A", "C", "G", "T", "ACTGA"],
             [[1, 2], [3, 4], [3, 4], [5, 6], [5, 6], [7, 8], [7, 8], [9, 10], [9, 10], [11], [11], []],
             [0, 1, 3, 5, 7, 9, 11], 4, 6)
            ])
def test_kmer_index(nodes, edges, ref, k, max_var):
    graph = Graph.from_sequence_edge_lists(nodes, edges, ref=ref)
    res_kmers, res_nodes = graph.create_kmer_index(k, max_variant_nodes=max_var)
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
    assert len(res_kmers) == len(ob_kmers)
    assert len(res_nodes) == len(ob_nodes)
    assert len(res_kmers) == len(res_nodes)
    res_counts = {}
    for i in range(len(res_kmers)):
        if res_kmers[i] not in res_counts:
            res_counts[res_kmers[i]] = 0
        res_counts[res_kmers[i]] += 1
    ob_counts = {}
    for i in range(len(ob_kmers)):
        if ob_kmers[i] not in ob_counts:
            ob_counts[ob_kmers[i]] = 0
        ob_counts[ob_kmers[i]] += 1
    for i in res_counts:
        print(kmer_hash_to_sequence(i, k),
              res_counts[i] if i in res_counts else 0,
              ob_counts[i] if i in ob_counts else 0)
    for i in res_counts:
        assert res_counts[i] == ob_counts[i]
