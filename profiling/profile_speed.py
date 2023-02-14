import biocy as bc

def main():
    nodes = ["AGTA", "G", "CT", "A", "CTA", "G", "A", "T"]
    edges = [[1, 2], [3], [3], [4], [5, 6], [7], [7], []]
    graph = bc.Graph.from_sequence_edge_lists(nodes, edges)
    index = graph.create_kmer_index(5) 
