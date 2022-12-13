#ifndef KMER_FINDER_H
#define KMER_FINDER_H

#include "Graph.hpp"
#include <stdint.h>

class KmerFinder {
public:
	Graph *graph;
	const uint8_t k;
	const uint8_t max_variant_nodes;

	uint64_t *found_kmers;
	uint32_t *found_nodes;
	uint64_t found_count;

private:
	uint64_t found_len;
	uint64_t kmer_buffer;
	uint64_t kmer_buffer_ext;
	uint32_t *path_buffer;
	uint16_t path_buffer_len;
	uint64_t kmer_mask;
	uint8_t kmer_buffer_shift;
	uint8_t variant_counter;

public:
	KmerFinder(Graph *graph, uint8_t k, uint8_t max_variant_nodes);
	~KmerFinder() {
		free(found_kmers);
		free(found_nodes);
		free(path_buffer);
	}

	void Find();
	void ReverseFoundKmers();

private:
	void FindKmersFromNode(uint32_t node_id);
	void FindKmersExtendedByEdge(uint32_t node_id, uint8_t kmer_len, uint8_t kmer_ext_len);
	void AddNodeFoundKmer(uint32_t node_id, uint64_t kmer);
};

#endif
