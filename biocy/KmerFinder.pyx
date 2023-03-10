# distutils: language = c++

from libcpp.unordered_map cimport unordered_map
from npstructures import RaggedArray, Counter

cdef class KmerFinder:
    cdef Graph graph
    cdef int k
    cdef int max_variant_nodes
    cdef public object _kmer_frequency_index
    cdef public object _kmer_frequency_index_starts
    cdef public object _kmer_frequency_index_lengths

    def __cinit__(self, Graph graph, int k, int max_variant_nodes=255):
        if k < 1 or k > 31:
            raise "KmerFinder: k must be between 1 and 31 inclusive"
        self.graph = graph
        self.k = k
        self.max_variant_nodes = max_variant_nodes
        self._kmer_frequency_index = None
        self._kmer_frequency_index_starts = None
        self._kmer_frequency_index_lengths = None

    def find(self, reverse_kmers=False):
        print("Finding kmers...")
        cdef cpp.KmerFinder *kf = new cpp.KmerFinder(self.graph.data, self.k, self.max_variant_nodes)
        kf.Find()
        if reverse_kmers:
            kf.ReverseFoundKmers()
        print("Copying to numpy arrays")
        kmers = np.empty((kf.found_count,), dtype=np.uint64)
        nodes = np.empty((kf.found_count,), dtype=np.uint32)
        cdef cnp.ndarray[unsigned long long, ndim=1, mode="c"] c_kmers = kmers
        cdef cnp.ndarray[unsigned int, ndim=1, mode="c"] c_nodes = nodes
        memcpy(c_kmers.data, kf.found_kmers, sizeof(unsigned long long) * kf.found_count)
        memcpy(c_nodes.data, kf.found_nodes, sizeof(unsigned int) * kf.found_count)
        del kf
        print("Done")
        return kmers, nodes

    def find_variant_signatures(self, reference_node_ids, variant_node_ids, reverse_kmers=False):
        if len(reference_node_ids) != len(variant_node_ids):
            raise "find_identifying_windows_for_variants: reference_node_ids and variant_node_ids must have the same length."
        cdef cpp.KmerFinder *kf = new cpp.KmerFinder(self.graph.data, self.k, self.max_variant_nodes)
        if self._kmer_frequency_index is not None:
            self.set_kmer_finder_frequency_index(kf)
        else:
            print("Finding kmers...")
            kf.SetFlag(cpp.FLAG_ONLY_SAVE_INITIAL_NODES, True)
            kf.Find()
            print("Creating frequency index...")
            kf.SetKmerFrequencyIndex(kf.CreateKmerFrequencyIndex())
            kf.SetFlag(cpp.FLAG_ONLY_SAVE_INITIAL_NODES, False)
        
        print("Finding windows...")
        cdef cpp.VariantWindow *window
        cdef uint32_t i
        cdef uint32_t j
        cdef uint32_t pair_count = len(reference_node_ids)
        cdef uint32_t reference_kmer_index = 0
        cdef uint32_t variant_kmer_index = 0
        cdef uint32_t additional_reference_kmers = 0
        cdef uint32_t additional_variant_kmers = 0
        reference_kmers = np.empty((pair_count,), dtype=np.uint64)
        variant_kmers = np.empty((pair_count,), dtype=np.uint64)
        reference_kmer_lens = np.empty((pair_count,), dtype=np.uint32)
        variant_kmer_lens = np.empty((pair_count,), dtype=np.uint32)
        
        cdef cpp.KmerFinder *window_finder = kf.CreateWindowFinder()

        for i in range(pair_count):
            window = kf.FindVariantSignaturesWithFinder(reference_node_ids[i], variant_node_ids[i], window_finder)
            
            if reverse_kmers:
                window.ReverseKmers(kf.k)

            reference_kmer_lens[i] = window.reference_kmers_len
            variant_kmer_lens[i] = window.variant_kmers_len

            additional_reference_kmers += window.reference_kmers_len - 1
            additional_variant_kmers += window.variant_kmers_len - 1

            if reference_kmer_index + window.reference_kmers_len > len(reference_kmers):
                reference_kmers.resize((len(reference_kmers) + additional_reference_kmers,))
                additional_reference_kmers = 0
            if variant_kmer_index + window.variant_kmers_len > len(variant_kmers):
                variant_kmers.resize((len(variant_kmers) + additional_variant_kmers,))
                additional_variant_kmers = 0

            for j in range(window.reference_kmers_len):
                reference_kmers[reference_kmer_index] = window.reference_kmers[j]
                reference_kmer_index += 1
            for j in range(window.variant_kmers_len):
                variant_kmers[variant_kmer_index] = window.variant_kmers[j]
                variant_kmer_index += 1
            
            del window

        del kf
        
        reference_kmers = RaggedArray(reference_kmers, shape=reference_kmer_lens, dtype=np.uint64)
        variant_kmers = RaggedArray(variant_kmers, shape=variant_kmer_lens, dtype=np.uint64)
        

        return reference_kmers, variant_kmers

    cdef set_kmer_finder_frequency_index(self, cpp.KmerFinder *kf):
        self._kmer_frequency_index_starts = np.ascontiguousarray(self._kmer_frequency_index._values._shape.starts)
        self._kmer_frequency_index_lengths = np.ascontiguousarray(self._kmer_frequency_index._values.lengths)
        cdef cnp.ndarray[long, ndim=1, mode="c"] c_keys = self._kmer_frequency_index._keys.ravel()
        cdef cnp.ndarray[long, ndim=1, mode="c"] c_values = self._kmer_frequency_index._values.ravel()
        cdef cnp.ndarray[long, ndim=1, mode="c"] c_starts = self._kmer_frequency_index_starts
        cdef cnp.ndarray[long, ndim=1, mode="c"] c_lengths = self._kmer_frequency_index_lengths
        cdef cpp.NPStructuresHashTable *hash_table = new cpp.NPStructuresHashTable()
        hash_table.mod = self._kmer_frequency_index._mod
        hash_table.keys = c_keys.data
        hash_table.values = c_values.data
        hash_table.starts = c_starts.data
        hash_table.lengths = c_lengths.data;
        hash_table.sizeof_key_dtype = self._kmer_frequency_index._keys.dtype.itemsize
        hash_table.sizeof_value_dtype = self._kmer_frequency_index._values.dtype.itemsize
        hash_table.sizeof_start_dtype = self._kmer_frequency_index_starts.dtype.itemsize
        hash_table.sizeof_length_dtype = self._kmer_frequency_index_lengths.dtype.itemsize
        if kf.nps_frequency_index != NULL:
            del kf.nps_frequency_index
        kf.nps_frequency_index = hash_table

    def set_kmer_frequency_index(self, frequency_index):
        if isinstance(frequency_index, Counter):
            self._kmer_frequency_index = frequency_index
        else:
            raise "A frequency index is required to be a Counter from npstructures."

