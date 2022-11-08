from obgraph import Graph as OBGraph
from graph_kmer_index.kmer_finder import DenseKmerFinder
from graph_kmer_index import kmer_hash_to_sequence
from Graph import Graph
import pytest

@pytest.mark.parametrize("nodes,edges,k", [
            (["ACTGACTGACTG", "ACTAGTC"], [[1], []], 3),
            (["ACTGACTGACTG", "ACTGCA"], [[1], []], 4),
            (["ACTAG", "GAC", "TGA", "CTGAGT"], [[1], [2], [3], []], 3),
            (["ACTAG", "G", "A", "GTACTCA"], [[1, 2], [3], [3], []], 3),
            (["AGTAGA", "G", "CT", "ACTA", "G", "A", "TCATA"], [[1, 2], [3], [3], [4, 5], [6], [6], []], 3),
            (["AGTAGA", "G", "CT", "ACTA", "G", "A", "TCATA"], [[1, 2], [3], [3], [4, 5], [6], [6], []], 4),
            (["AGTAGA", "G", "CT", "ACTA", "G", "A", "TCATA"], [[1, 2], [3], [3], [4, 5], [6], [6], []], 5)])
def test_kmer_index(nodes, edges, k):
    graph = Graph.from_sequence_edge_lists(nodes, edges)
    res_kmers, res_nodes = graph.create_kmer_index(k)
    ob_node_sequences = {}
    ob_edges = {}
    ob_linear_ref_nodes = []
    for i in range(len(nodes)):
        ob_node_sequences[i + 1] = nodes[i]
        if len(edges[i]) > 0:
            ob_edges[i + 1] = []
            for edge in edges[i]:
                ob_edges[i + 1].append(edge + 1)
    ob_linear_ref_nodes.append(1)
    while ob_linear_ref_nodes[-1] in ob_edges:
        ob_linear_ref_nodes.append(ob_edges[ob_linear_ref_nodes[-1]][0])
    obgraph = OBGraph.from_dicts(
        node_sequences=ob_node_sequences,
        edges=ob_edges,
        linear_ref_nodes=ob_linear_ref_nodes
    )
    finder = DenseKmerFinder(obgraph, k=k, max_variant_nodes=100000)
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
