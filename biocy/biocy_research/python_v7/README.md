Large rewrite that stores the graph node sequences as 2-bit encoded hashes instead of bytestrings of ASCII characters.
Supports nodes with infinitely long sequences and empty nodes.
This implementation allows for very fast kmer buffer management with left/right shift and bitwise and/or.
It also saves the additional work previously done where every base of a kmer had to be hashed when added to the index.
The current implementation uses only two 64-bit buffer variables. One of which is named `kmer_buffer_ext` and is used exclusively for recursive calls.
While it would be possible to create an even faster implementation with a single buffer, this solution would only support up to a max of 15-mers before the recursive calls don't have enough buffer space. As large kmers are a priority for this program's use-case, and small kmers are already very fast due to low complexity, implementing this separate solution for 1-mers through 15-mers is not a priority.
