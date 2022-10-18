from obgraph import Graph as OBGraph
from graph_kmer_index.kmer_finder import DenseKmerFinder
from Graph import Graph
import pytest

@pytest.mark.parametrize("nodes,edges,k", [
            (["ACTGACTGACTG", "ACTAGTC"], [[1], []], 3),
            (["ACTGACTGACTG", "ACTGCA"], [[1], []], 4),
            (["ACTAG", "GAC", "TGA", "CTGAGT"], [[1], [2], [3], []], 3),
            (["ACTAG", "G", "A", "GTACTCA"], [[1, 2], [3], [3], []], 3),
            (["AGTAGA", "G", "CT", "A", "CTA", "G", "A", "TCATA"], [[1, 2], [3], [3], [4], [5, 6], [7], [7], []], 3),
            (["AGTAGA", "G", "CT", "A", "CTA", "G", "A", "TCATA"], [[1, 2], [3], [3], [4], [5, 6], [7], [7], []], 4),
            (["AGTAGA", "G", "CT", "A", "CTA", "G", "A", "TCATA"], [[1, 2], [3], [3], [4], [5, 6], [7], [7], []], 5)])
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
    finder = DenseKmerFinder(obgraph, k=k)
    finder.find()
    ob_kmers, ob_nodes = finder.get_found_kmers_and_nodes()
    assert len(res_kmers) == len(ob_kmers)
    assert len(res_nodes) == len(ob_nodes)
    assert len(res_kmers) == len(res_nodes)
    res_kmers = sorted(res_kmers)
    res_nodes = sorted(res_nodes)
    ob_kmers = sorted(ob_kmers)
    ob_nodes = sorted(ob_nodes)
    for i in range(len(res_kmers)):
        assert res_kmers[i] == ob_kmers[i]
        assert res_nodes[i] + 1 == ob_nodes[i]
    
