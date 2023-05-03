---
layout: default
title: KIVS - Graph
permalink: /docs/Graph/
---

# class Graph

`from kivs import Graph`

Representation of a directional genome graph of reference nodes and variant nodes.

## Method Summary

### Static Methods

| **Return Type** | **Method and Description** |
|---|---|
| Graph | `from_file(filepath: str)`<br>Load a Graph from a file at *filepath*.<br>The file extension for BioCy graphs is ".bcg". |
| Graph | `from_obgraph(graph: obgraph.Graph)`<br>Requires the obgraph module.<br>Construct a Graph from an obgraph. |
| Graph | `from_gfa(filepath: str)`<br>Construct a Graph from a GFA file at *filepath*. |
| Graph | `from_sequence_edge_lists(sequences, edges, encoding="ACGT", ref=None)`<br>Construct a Graph using the given a list of sequences and edges.<br>Parameters:<br>- sequences: List of strings of bases, such as "ACCGTA".<br>- edges: 2D list of node indexes in the format `edges[from_node] = [to_node1, to_node2, ...]`<br>- *[encoding]*: String containing A, C, T and G. These bases will be encoded as 0-3 in order.<br>- *[ref]*: List of nodes that are reference nodes (not variant nodes). If not specified, all nodes are treated as reference nodes. |
| Boolean | `ìs_valid_encoding(encoding)`<br>Returns whether an encoding is valid.<br>For an encoding to be valid, it must be a permutation of "ACGT". |

### Instance Methods

| **Return Type** | **Method and Description** |
|---|---|
| None | `to_file(filepath: str)`<br>Saves a Graph to a file at *filepath*.<br>The file extension for BioCy graphs is ".bcg". |
| obgraph.Graph | `to_obgraph()`<br>Requires the obgraph module.<br>Construct an obgraph from this BioCy graph. |
| None | `print_node_data(node_id: int)`<br>Prints various information about the node with id *node_id*. |
