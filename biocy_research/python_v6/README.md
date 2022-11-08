Added support for empty nodes.
Accomplished this simply by expanding the node path buffer from length k to length k * 4.
This is assuming a worst-case where there are no more than 3 empty nodes for each base in the kmer.
Given that k will only support values up to 31, scaling this buffer further should be no issue if necessary.
Also added a thorough test for this functionality.
