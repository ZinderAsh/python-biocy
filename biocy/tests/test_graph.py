import biocy as bc
import pytest

@pytest.mark.parametrize("nodes,edges,k,expected", [
            (["ACTGACTGACTG"], [[]], 3,
             {b'ACT': [0, 0, 0], b'CTG': [0, 0, 0], b'TGA': [0, 0], b'GAC': [0, 0]}),
            (["ACTGACTGACTG"], [[]], 4,
             {b'ACTG': [0, 0, 0], b'CTGA': [0, 0], b'TGAC': [0, 0], b'GACT': [0, 0]}),
            (["ACT", "GAC", "TGA", "CTG"], [[1], [2], [3], []], 3,
             {b'ACT': [0, 1, 2], b'CTG': [0, 1, 3], b'TGA': [0, 2], b'GAC': [1, 2]}),
            (["ACT", "G", "A", "GTAC"], [[1, 2], [3], [3], []], 4,
             {b'ACTG': [0], b'CTGG': [0], b'TGGT': [0], b'GGTA': [1], b'ACTA': [0], b'CTAG': [0], b'TAGT': [0], b'AGTA': [2], b'GTAC': [3]}),
            (["AGTA", "G", "CT", "A", "CTA", "G", "A", "T"], [[1, 2], [3], [3], [4], [5, 6], [7], [7], []], 3,
             {b'AGT': [0, 4], b'GTA': [0], b'TAG': [0, 4], b'TAC': [0, 2], b'AGA': [0], b'ACT': [0, 3], b'GAC': [1], b'CTA': [4, 2], b'TAA': [4], b'AAT': [4]}),
            (["AGTA", "G", "CT", "A", "CTA", "G", "A", "T"], [[1, 2], [3], [3], [4], [5, 6], [7], [7], []], 4,
             {b'AGTA': [0], b'GTAG': [0], b'GTAC': [0], b'TAGA': [0], b'TACT': [0, 2], b'AGAC': [0], b'ACTA': [0, 3], b'GACT': [1], b'CTAG': [4], b'CTAA': [4], b'TAGT': [4], b'TAAT': [4], b'CTAC': [2]}),
            (["AGTA", "G", "CT", "A", "CTA", "G", "A", "T"], [[1, 2], [3], [3], [4], [5, 6], [7], [7], []], 5,
             {b'AGTAG': [0], b'AGTAC': [0], b'GTAGA': [0], b'GTACT': [0], b'TAGAC': [0], b'TACTA': [0, 2], b'AGACT': [0], b'ACTAC': [0], b'GACTA': [1], b'ACTAG': [3], b'ACTAA': [3], b'CTAGT': [4], b'CTAAT': [4], b'CTACT': [2]})])
def test_kmer_index(nodes, edges, k, expected):
    graph = bc.Graph.from_sequence_edge_lists(nodes, edges)
    index = graph.create_kmer_index(k)
    assert len(index) == len(expected)
    for i in index:
        assert len(index[i]) == len(expected[i])
        for a, b in zip(index[i], expected[i]):
            assert a == b
