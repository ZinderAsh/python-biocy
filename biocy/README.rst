=====
BioCy
=====

Master's project about discovering and indexing k-mers from genome graphs.

## Description

This project aims to create a high-performance cython-based python module for creating genome graphs and discovering/indexing the k-mers within the graph.

## Getting Started

### Dependencies

* Python
* Cython
* PyTest (for running tests)

### Setup

* Clone the repository: `git clone https://github.com/ZinderAsh/masters-project-genotyping.git`
* cd into the directory: `cd masters-project-genotyping`
* Run `make` inside the folder

### Usage

* Run `python example.py` for an example output.
* Run `pytest` to perform the tests ensuring that everything is working correctly.
* To use in a separate Python program:
```python
from Graph import Graph

graph = Graph.from_sequence_edge_lists(
	["ACT", "TGA"], # List of all nodes to be in the graph
	[[1], []]) # Indexes of nodes each node is connected to

graph.print_graph() # The graph can be visualized

index = graph.create_kmer_index(3) # Takes k for k-mer length as parameter

print(index[b'ACT']) # Index dictionary uses bytestrings as keys
```

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
