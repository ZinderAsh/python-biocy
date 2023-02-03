# distutils: language = c++

cimport biocy.biocpp as cpp

include "Graph.pyx"
include "KmerFinder.pyx"

def test_func():
    return "Hello"
