# distutils: language = c++

cdef class KmerFinder:
    cdef Graph graph
    cdef int k
    cdef int max_variant_nodes

    def __cinit__(self, Graph graph, int k, int max_variant_nodes=31):
        if k < 1 or k > 31:
            raise "KmerFinder: k must be between 1 and 31 inclusive"
        self.graph = graph
        self.k = k
        self.max_variant_nodes = max_variant_nodes

    def find(self, reverse_kmers=False):
        print("Finding kmers...")
        cdef cpp.KmerFinder *kf = new cpp.KmerFinder(self.graph.data, self.k, self.max_variant_nodes)
        kf.Find()
        if reverse_kmers:
            kf.ReverseFoundKmers()
        print("Copying to numpy arrays")
        kmers = np.empty((kf.found_count,), dtype=np.ulonglong)
        nodes = np.empty((kf.found_count,), dtype=np.uint32)
        cdef cnp.ndarray[unsigned long long, ndim=1, mode="c"] c_kmers = kmers
        cdef cnp.ndarray[unsigned int, ndim=1, mode="c"] c_nodes = nodes
        memcpy(c_kmers.data, kf.found_kmers, sizeof(unsigned long long) * kf.found_count)
        memcpy(c_nodes.data, kf.found_nodes, sizeof(unsigned int) * kf.found_count)
        del kf
        print("Done")
        return kmers, nodes
