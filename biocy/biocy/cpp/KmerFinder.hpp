#ifndef KMER_FINDER_H
#define KMER_FINDER_H

#include "Graph.hpp"
#include <stdint.h>
#include <unordered_map>

#define FILTER_NODE_ID 1

class KmerFinder {
public:
	Graph *graph;
	const uint8_t k;
	const uint8_t max_variant_nodes;

	bool save_sequence_start_positions;
	bool save_sequence_kmer_positions;

	uint64_t *found_kmers;
	uint32_t *found_nodes;
	uint32_t *found_node_sequence_start_positions;
	uint16_t *found_node_sequence_kmer_positions;
	uint64_t found_count;

private:
	uint64_t found_len;
	uint64_t kmer_buffer;
	uint64_t kmer_buffer_ext;
	uint32_t *path_buffer;
	uint32_t start_position;
	uint16_t *kmer_position_buffer;
	uint16_t path_buffer_len;
	uint64_t kmer_mask;
	uint8_t kmer_buffer_shift;
	uint8_t variant_counter;
	uint8_t filters;
	uint32_t filter_node_id;

public:
	KmerFinder(Graph *graph, uint8_t k, uint8_t max_variant_nodes);
	~KmerFinder() {
		if (found_kmers) free(found_kmers);
		if (found_nodes) free(found_nodes);
		if (found_node_sequence_start_positions) free(found_node_sequence_start_positions);
		if (found_node_sequence_kmer_positions) free(found_node_sequence_kmer_positions);
		free(path_buffer);
		free(kmer_position_buffer);
	}
	
	void Reset();
	void InitializeFoundArrays();
	void Find();
	void FindKmersForVariant(uint32_t reference_node_id, uint32_t variant_node_id);
	void FindKmersSpanningNode(uint32_t center_node_id);
	void ReverseFoundKmers();
	std::unordered_map<uint64_t, uint32_t> CreateKmerIndex();
	
	void SetFilter(uint8_t filter, uint64_t value) {
		filters |= filter;
		switch (filter) {
			case FILTER_NODE_ID:
				filter_node_id = value;
				break;
			default:
				break;
		}
	}
	void RemoveFilter(uint8_t filter) {
		filters &= ~filter;
	}

private:
	uint64_t FindKmersFromNode(uint32_t node_id);
	uint64_t FindKmersExtendedByEdge(uint32_t node_id, uint8_t kmer_len, uint8_t kmer_ext_len);
	bool AddNodeFoundKmer(uint32_t node_id, uint64_t kmer, uint32_t start_position, uint16_t node_kmer_position);
};

#endif
