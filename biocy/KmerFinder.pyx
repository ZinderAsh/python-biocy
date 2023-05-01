# distutils: language = c++

from libcpp cimport bool
from libcpp.unordered_map cimport unordered_map
from cython.operator import dereference, postincrement
from npstructures import RaggedArray, Counter

cdef class KmerFinder:
    cdef Graph graph
    cdef int k
    cdef bool reverse_kmers
    cdef public object _kmer_frequency_index

    def __cinit__(self, Graph graph, int k, bool reverse_kmers=False):
        if k < 1 or k > 31:
            raise "KmerFinder: k must be between 1 and 31 inclusive"
        self.graph = graph
        self.k = k
        self.reverse_kmers = reverse_kmers
        self._kmer_frequency_index = None

    def find(self, bool include_spanning_nodes=False, int max_variant_nodes=255, bool stdout=False):
        #print("Finding kmers...")
        cdef cpp.KmerFinder *kf = new cpp.KmerFinder(self.graph.data, self.k, max_variant_nodes)
        kf.SetFlag(cpp.FLAG_ONLY_SAVE_INITIAL_NODES, not include_spanning_nodes)
        kf.SetFlag(cpp.FLAG_TO_STDOUT, stdout)
        kf.Find()
        if self.reverse_kmers:
            kf.ReverseFoundKmers()
        #print("Copying to numpy arrays")
        if stdout:
            del kf
            return None
        kmers = np.empty((kf.found_count,), dtype=np.uint64)
        nodes = np.empty((kf.found_count,), dtype=np.uint32)
        cdef cnp.ndarray[unsigned long long, ndim=1, mode="c"] c_kmers = kmers
        cdef cnp.ndarray[unsigned int, ndim=1, mode="c"] c_nodes = nodes
        memcpy(c_kmers.data, kf.found_kmers, sizeof(unsigned long long) * kf.found_count)
        memcpy(c_nodes.data, kf.found_nodes, sizeof(unsigned int) * kf.found_count)
        del kf
        #print("Done")
        return kmers, nodes

    def find_kmers_spanning_node(self, int node_id, int max_variant_nodes=255, bool stdout=False):
        #print("Finding kmers...")
        cdef cpp.KmerFinder *kf = new cpp.KmerFinder(self.graph.data, self.k, max_variant_nodes)
        kf.SetFlag(cpp.FLAG_TO_STDOUT, stdout)
        kf.FindKmersSpanningNode(node_id)
        if self.reverse_kmers:
            kf.ReverseFoundKmers()
        #print("Copying to numpy arrays")
        if stdout:
            del kf
            return None
        kmers = np.empty((kf.found_count,), dtype=np.uint64)
        nodes = np.empty((kf.found_count,), dtype=np.uint32)
        cdef cnp.ndarray[unsigned long long, ndim=1, mode="c"] c_kmers = kmers
        cdef cnp.ndarray[unsigned int, ndim=1, mode="c"] c_nodes = nodes
        memcpy(c_kmers.data, kf.found_kmers, sizeof(unsigned long long) * kf.found_count)
        memcpy(c_nodes.data, kf.found_nodes, sizeof(unsigned int) * kf.found_count)
        del kf
        #print("Done")
        return kmers, nodes

    def find_variant_signatures(self, reference_node_ids, variant_node_ids,
                                int max_variant_nodes=255, bool minimize_overlaps=False, bool align_windows=False):
        if len(reference_node_ids) != len(variant_node_ids):
            raise "find_identifying_windows_for_variants: reference_node_ids and variant_node_ids must have the same length."
        cdef cpp.KmerFinder *kf = new cpp.KmerFinder(self.graph.data, self.k, max_variant_nodes)
        if self._kmer_frequency_index is not None:
            print("Creating frequency index from Counter...")
            self.set_kmer_finder_frequency_index(kf)
        else:
            self.create_frequency_index(max_variant_nodes=max_variant_nodes)
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
        
        kf.SetFlag(cpp.FLAG_MINIMIZE_SIGNATURE_OVERLAP, minimize_overlaps)
        kf.SetFlag(cpp.FLAG_ALIGN_SIGNATURE_WINDOWS, align_windows)
        cdef cpp.KmerFinder *window_finder = kf.CreateWindowFinder()

        for i in range(pair_count):
            window = kf.FindVariantSignaturesWithFinder(reference_node_ids[i], variant_node_ids[i], window_finder)
            
            if self.reverse_kmers:
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

        del window_finder
        del kf
        
        reference_kmers = RaggedArray(reference_kmers, shape=reference_kmer_lens, dtype=np.uint64)
        variant_kmers = RaggedArray(variant_kmers, shape=variant_kmer_lens, dtype=np.uint64)
        

        return reference_kmers, variant_kmers 

    def create_frequency_index(self, max_variant_nodes=255, set_index=True):
        cdef cpp.KmerFinder *kf = new cpp.KmerFinder(self.graph.data, self.k, max_variant_nodes)
        cdef unordered_map[uint64_t, uint32_t] frequency_index

        print("Finding kmers...")
        kf.SetFlag(cpp.FLAG_ONLY_SAVE_INITIAL_NODES, True)
        kf.Find()
        if self.reverse_kmers:
            kf.ReverseFoundKmers()

        print("Creating frequency index...")
        frequency_index = kf.CreateKmerFrequencyIndex()
        del kf

        py_frequency_index = {}
        cdef unordered_map[uint64_t, uint32_t].iterator it = frequency_index.begin()
        while it != frequency_index.end():
            py_frequency_index[dereference(it).first] = dereference(it).second
            postincrement(it)

        if set_index:
            self._kmer_frequency_index = py_frequency_index

        return py_frequency_index
        

    def set_frequency_index(self, frequency_index):
        if isinstance(frequency_index, Counter):
            self._kmer_frequency_index = frequency_index
        elif isinstance(frequency_index, dict):
            self._kmer_frequency_index = frequency_index
        else:
            raise "A frequency index is required to be a Python dictionary or a Counter from npstructures."

    cdef set_kmer_finder_frequency_index(self, cpp.KmerFinder *kf):
        if isinstance(self._kmer_frequency_index, Counter):
            self.set_kmer_finder_frequency_index_from_counter(kf)
        elif isinstance(self._kmer_frequency_index, dict):
            self.set_kmer_finder_frequency_index_from_dict(kf)

    cdef set_kmer_finder_frequency_index_from_counter(self, cpp.KmerFinder *kf):
        cdef unordered_map[uint64_t, uint32_t] frequency_index;
        
        kmer_frequency_index_keys = np.ascontiguousarray(self._kmer_frequency_index._keys.ravel()).astype(np.uint64)
        kmer_frequency_index_values = np.ascontiguousarray(self._kmer_frequency_index._values.ravel()).astype(np.uint32)

        cdef cnp.ndarray[uint64_t, ndim=1, mode="c"] c_keys = kmer_frequency_index_keys
        cdef cnp.ndarray[uint32_t, ndim=1, mode="c"] c_values = kmer_frequency_index_values

        cdef uint64_t i
        cdef uint64_t key_count = len(kmer_frequency_index_keys)
        cdef uint64_t key
        for i in range(key_count):
            key = c_keys[i]
            if self.reverse_kmers:
                key = hashing.reverse_kmer(key, self.k)
            frequency_index[key] = c_values[i]
        
        kf.SetKmerFrequencyIndex(frequency_index)

    cdef set_kmer_finder_frequency_index_from_dict(self, cpp.KmerFinder *kf):
        cdef unordered_map[uint64_t, uint32_t] frequency_index;
        
        cdef uint64_t key
        cdef uint32_t value
        for key, value in self._kmer_frequency_index.items():
            if self.reverse_kmers:
                key = hashing.reverse_kmer(key, self.k)
            frequency_index[key] = value
        
        kf.SetKmerFrequencyIndex(frequency_index)
