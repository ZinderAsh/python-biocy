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
| `KmerFinder(graph: Graph, k: int, reverse_kmers=False)`<br>Initialize a KmerFinder ready for further function calls.<br>Parameters:<br>- graph: The biocy Graph to analyze.<br>- k: integer between 1 and 31 inclusive. The length of the kmers to analyze.<br>- *[reverse\_kmers]*: If set to `True`, all k-mer hashes will be reversed before returning. |

### Instance Methods

| **Return Type** | **Method and Description** |
|---|---|
| np.array(np.uint64),<br>np.array(np.uint32) | `find(max_variant_nodes=255, include_spanning_nodes=False, stdout=False)`<br>Returns two NumPy arrays of equal length.<br>Explores the graph and returns all k-mers (hashed) in the first array and the nodes they were found at in the second array.<br>Parameters:<br>- *[include_spanning_nodes]*: If set to True and a k-mer spans more than one node, the k-mer will have one entry in the arrays for each node it spans.<br>- *[max\_variant\_nodes]*: Ignore paths that cross through at least *max_variant_nodes* variant nodes.<br>- *[stdout]*: If set to True, results will be printed to stdout instead of returning as arrays. |
| np.array(np.uint64) | `find_kmers_spanning_node(node_id: int, max_variant_nodes=255, include_spanning_nodes=False, stdout=False)`<br>Explores the graph and returns all k-mers (hashed) that include the specified node.<br>Parameters:<br>- node\_id: The node to find k-mers for.<br>- *[include_spanning_nodes]*: If set to True and a k-mer spans more than one node, the k-mer will have one entry in the arrays for each node it spans.<br>- *[max\_variant\_nodes]*: Ignore paths that cross through at least *max_variant_nodes* variant nodes.<br>- *[stdout]*: If set to True, results will be printed to stdout instead of returning as arrays. |
| RaggedArray(np.uint64),<br>RaggedArray(np.uint64) | `find_variant_signatures(ref_node_ids, var_node_ids, max_variant_nodes=255, minimize_overlaps=False, align_windows=False)`<br>Returns two RaggedArrays of hashed k-mers with equal length.<br>For each pair of reference node and variant node provided, finds a window of k length spanning both nodes with the rarest k-mers in the Graph, and returns that k-mer. If there are multiple possible paths through the graph for the found window, all possible k-mers for that window are returned.<br>Parameters:<br>- ref\_node\_ids: A list of all reference nodes to find windows for.<br>- var\_node\_ids: A list of all variant nodes to find windows for. Must be same length as ref\_node\_ids.<br>- *[max\_variant\_nodes]*: Ignore paths that cross through at least *max_variant_nodes* variant nodes.<br>- *[minimize\_overlaps]*: If True, the variant signature will attempt to avoid having k-mers that were reference signature candidates, and vice versa.<br>- *[align\_windows]*: If True, the variant signature and reference signature will be aligned, such that they span the same area of the genome. |
| dict | `create_frequency_index(max_variant_nodes=255, set_index=True)`<br>Creates and returns a dictionary where the keys are hashed k-mers and the values are their frequency in the graph.<br>Parameters:<br>- *[max\_variant\_nodes]*: Ignore paths that cross through at least *max_variant_nodes* variant nodes.<br>- *[set\_index]*: If True, the created index will also be set as this KmerFinder's index. |
| None | `set_frequency_index(index)`<br>Set a specific kmer frequency index for use when finding variant signatures. Only accepts Python dictionaries and npstructures.Counter objects. |
