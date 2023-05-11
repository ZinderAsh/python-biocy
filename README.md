# KIVS: Graph **K**-mer **I**ndexer and **V**ariant **S**ignature Finder

Master's project about discovering and indexing k-mers from genome graphs.

## Description

This project aims to create a high-performance cython-based python module for creating genome graphs and discovering/indexing the k-mers within the graph.

## Getting Started

### Dependencies

* Python
* Cython
* PyTest (for running tests)
* npstructures

#### Optional Dependencies

* obgraph

### Setup

* Clone the repository: `git clone https://github.com/ZinderAsh/masters-project-genotyping.git`
* cd into the directory: `cd masters-project-genotyping`
* Run `pip install .` inside the folder

### Usage

* Run `pytest` to perform the tests ensuring that everything is working correctly.
* To use in a separate Python program (example using obgraph):
```python
from kivs import Graph, KmerFinder
from obgraph import Graph as OBGraph

# Read graph from obgraph file
obgraph = OBGraph.from_file("directory/obgraph.npz")
graph = Graph.from_obgraph(obgraph)

# Write graph to biocy graph file
graph.to_file("directory/bcgraph.kivs")

# Read graph from biocy graph file
graph = Graph.from_file("directory/bcgraph.kivs")

# Find all 5-mers that span a maximum of 5 variants
k = 5
kmer_finder = KmerFinder(graph, k)
kmers, nodes = kmer_finder.find(max_variant_nodes=5)

# Find identifying windows for variants
# reverse_kmers specify if the kmers should be mirrored before returning
k = 5
kmer_finder = KmerFinder(graph, k, reverse_kmers=True)
ref, var = kmer_finder.find_variant_signatures(
    [1, 3, 5, 7, 9], [2, 4, 6, 8, 10],
    align_windows=True, minimize_overlaps=True)

```

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
