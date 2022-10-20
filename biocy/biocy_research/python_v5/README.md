Moved finding kmers to a separate class to...
- Keep most recursive parameters as class variables, adding less to the stack.
- Better readability

Renamed `kmer_buf` and `nodes_buf` to `kmer_buffer` and `path_buffer` for clarity. `nodes_len` has also been changed to `path_len`.
Removed the recurse parameter in favor of a recursive variable as `path_len > 1` as it amounts to the same thing.
Rewrote array shifting of the `kmer_buffer` to use numpy array slicing methods (still in-place).
Added spacing and better comments for readability.
