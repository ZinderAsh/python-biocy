---
layout: default
title: Biocy - KmerFinder
permalink: /docs/KmerFinder/
---

# class KmerFinder

`from biocy import KmerFinder`

Used to analyze and retrieve various information about k-mers in a Graph.

## Method Summary

### Constructors

| **Method and Description** |
|---|
| `KmerFinder(graph: Graph, k: int, max_variant_nodes=255)`<br>Initialize a KmerFinder ready for further function calls.<br>Parameters:<br>- graph: The biocy Graph to analyze.<br>- k: integer between 1 and 31 inclusive. The length of the kmers to analyze.<br>- *[max\_variant\_nodes]*: Ignore paths that cross through at least *max_variant_nodes* variant nodes. |

### Instance Methods

| **Return Type** | **Method and Description** |
|---|---|
| np.array(np.uint64),<br>np.array(np.uint32) | `find(reverse_kmers=False)`<br>Returns two NumPy arrays of equal length.<br>Explores the graph and returns all k-mers (hashed) in the first array and the nodes they were found at in the second array. If a k-mer spans more than one node, the k-mer will have one entry in the arrays for each node it spans.<br>If `reverse_kmers` is set to `True`, all k-mer hashes will be reversed before returning. |
| RaggedArray(np.uint64),<br>RaggedArray(np.uint64) | `find_variant_signatures(ref_node_ids, var_node_ids, reverse_kmers=False)`<br>Returns two RaggedArrays of hashed k-mers with equal length.<br>For each pair of reference node and variant node provided, finds a window of k length spanning both nodes with the rarest k-mers in the Graph, and returns that k-mer. If there are multiple possible paths through the graph for the found window, all possible k-mers for that window are returned.<br>Parameters:<br>- ref\_node\_ids: A list of all reference nodes to find windows for.<br>- var\_node\_ids: A list of all variant nodes to find windows for. Must be same length as ref\_node\_ids.<br>- *[reverse\_kmers]*: If True, reversed all k-mer hashes before returning. |
| None | `set_kmer_frequency_index(npstructures.Counter index)`<br>Set a specific kmer frequency index for use when finding variant signatures. |
