# distutils: language = c++

cimport biocy.biocpp as cpp
cimport biocy.hashing as hashing

try:
    import obgraph as ob
    has_obgraph = True
except ImportError:
    has_obgraph = False

include "Graph.pyx"
include "KmerFinder.pyx"

def test_func():
    return "Hello"
