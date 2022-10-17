from Graph_with_stats import Graph
import random

BASES = ["A", "C", "T", "G"]

def rand_sequence(seq_len=1):
    seq = ""
    for _ in range(seq_len):
        seq += BASES[int(random.random() * 4 // 1)]
    return seq

def create_matrix_graph(width, height, seq_len=1):
    nodes = []
    edges = []
    for _ in range(width * height):
        nodes.append(rand_sequence(seq_len))
        edges.append([])
    for w in range(width - 1):
        for h in range(height):
            for h2 in range(height):
                edges[w * height + h].append((w + 1) * height + h2)
    return Graph.from_sequence_edge_lists(nodes, edges)

def measure_k(f=None):
    for i in [8, 12, 16, 20]:
        graph = create_matrix_graph(i, i, seq_len=1)
        for k in range(2, 6 if i == 20 else 7):
            graph.create_kmer_index(k)
            print(f'N = {len(graph.nodes)}, E = {graph.total_edges}, k = {k}')
            print(f'  calls: {graph.get_kmers_initial_calls + graph.get_kmers_recursive_calls}')
            total_calls = graph.get_kmers_initial_calls + graph.get_kmers_recursive_calls
            if f: f.write(f'{len(graph.nodes)},{graph.total_edges},{k},{total_calls}\n')

def measure_n_e(f=None):
    k = 4
    for i in range(2, 32):
        graph = create_matrix_graph(i, i, seq_len=1)
        graph.create_kmer_index(k)
        print(f'N = {len(graph.nodes)}, E = {graph.total_edges}, k = {k}')
        print(f'  calls: {graph.get_kmers_initial_calls + graph.get_kmers_recursive_calls}')
        total_calls = graph.get_kmers_initial_calls + graph.get_kmers_recursive_calls
        if f: f.write(f'{len(graph.nodes)},{graph.total_edges},{k},{total_calls}\n')

f = open("measures.csv", "w")
f.write("N,E,k,calls\n")
measure_n_e(f)
measure_k(f)
#f.close()

def measure_n_e_2(f=None):
    k = 3
    for w, h in [(64, 1), (32, 2), (16, 4), (8, 8), (4, 16), (2, 32)]:
        graph = create_matrix_graph(w, h, seq_len=1)
        graph.create_kmer_index(k)
        print(f'{w}x{h}: N = {len(graph.nodes)}, E = {graph.total_edges}, k = {k}')
        print(f'  calls: {graph.get_kmers_initial_calls + graph.get_kmers_recursive_calls}')
        total_calls = graph.get_kmers_initial_calls + graph.get_kmers_recursive_calls
        if f: f.write(f'{len(graph.nodes)},{graph.total_edges},{k},{total_calls}\n')

measure_n_e_2(f)
f.close()
