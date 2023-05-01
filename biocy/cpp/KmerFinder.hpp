#ifndef KMER_FINDER_H
#define KMER_FINDER_H

#include "Graph.hpp"
#include "node.hpp"
#include <stdint.h>
#include <stdio.h>
#include <unordered_map>
#include <vector>

#define FILTER_NODE_ID 1

#define FLAG_TO_STDOUT 1
#define FLAG_ALIGN_SIGNATURE_WINDOWS 1 << 4
#define FLAG_MINIMIZE_SIGNATURE_OVERLAP 1 << 5
#define FLAG_ONLY_SAVE_INITIAL_NODES 1 << 6
#define FLAG_SAVE_WINDOWS 1 << 7

class VariantWindow;

struct kmer_window {
	uint32_t node_id;
	uint64_t *kmers;
	uint16_t length;
	uint32_t start_position;
	uint16_t kmer_position;
	uint32_t max_frequency;
};

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
	struct kmer_window *found_windows;
	uint32_t found_window_count;

private:
	uint64_t found_len;
	uint32_t found_window_len;
	uint64_t kmer_buffer;
	uint64_t kmer_buffer_ext;
	uint32_t *path_buffer;
	uint32_t start_position;
	uint16_t *kmer_position_buffer;
	uint16_t path_buffer_len;
	uint64_t kmer_mask;
	uint8_t kmer_buffer_shift;
	uint8_t variant_counter;

	uint8_t flags;

	uint8_t filters;
	uint32_t filter_node_id;

	std::unordered_map<uint64_t, uint32_t> kmer_frequency_index;

public:
	KmerFinder(Graph *graph, uint8_t k, uint8_t max_variant_nodes);
	~KmerFinder() {
		if (found_kmers) free(found_kmers);
		if (found_nodes) free(found_nodes);
		if (found_node_sequence_start_positions) free(found_node_sequence_start_positions);
		if (found_node_sequence_kmer_positions) free(found_node_sequence_kmer_positions);
		if (found_windows) {
			for (uint32_t i = 0; i < found_window_count; i++) {
				free((found_windows + i)->kmers);
			}
			free(found_windows);
		};
		free(path_buffer);
		free(kmer_position_buffer);
	}
	
	void Reset();
	void InitializeFoundArrays();
	void Find();
	void FindKmersForVariant(uint32_t reference_node_id, uint32_t variant_node_id);
	void FindKmersSpanningNode(uint32_t center_node_id);
	KmerFinder *CreateWindowFinder();
	std::vector<VariantWindow *> FindWindowsForVariant(uint32_t reference_node_id, uint32_t variant_node_id);
	std::vector<VariantWindow *> FindWindowsForVariantWithFinder(uint32_t reference_node_id, uint32_t variant_node_id, KmerFinder *kf);
	VariantWindow *FindVariantSignatures(uint32_t reference_node_id, uint32_t variant_node_id);
	VariantWindow *FindVariantSignaturesWithFinder(uint32_t reference_node_id, uint32_t variant_node_id, KmerFinder *kf);
	VariantWindow *FindAlignedSignaturesWithFinder(uint32_t reference_node_id, uint32_t variant_node_id, KmerFinder *kf);
	VariantWindow *FindUnalignedSignaturesWithFinder(uint32_t reference_node_id, uint32_t variant_node_id, KmerFinder *kf);
	uint32_t GetWindowOverlap(std::vector<VariantWindow *> *windows, uint32_t window_index);
	void ReverseFoundKmers();
	std::unordered_map<uint64_t, uint32_t> CreateKmerFrequencyIndex();
	uint64_t GetKmerFrequency(uint64_t kmer);

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
	void SetFlag(uint8_t flag, bool value) {
		flags = value ? (flags | flag) : (flags & ~flag);
	}
	void SetKmerFrequencyIndex(std::unordered_map<uint64_t, uint32_t> index) {
		kmer_frequency_index = index;
	}
	bool HasKmerFrequencyIndex() {
		return !kmer_frequency_index.empty();
	}

private:
	uint64_t FindKmersFromNode(uint32_t node_id);
	uint64_t FindKmersExtendedByEdge(uint32_t node_id, uint8_t kmer_len, uint8_t kmer_ext_len);
	bool AddNodeFoundKmer(uint32_t node_id, uint64_t kmer, uint32_t start_position, uint16_t node_kmer_position);
	bool AddFoundWindowKmer(uint32_t node_id, uint64_t kmer, uint32_t start_position, uint16_t node_kmer_positions);
};

class VariantWindow {
public:
	uint64_t *reference_kmers;
	uint64_t *variant_kmers;
	uint16_t reference_kmers_len;
	uint16_t variant_kmers_len;
	uint32_t max_frequency;
	
	VariantWindow(struct kmer_window *ref, struct kmer_window *var) {
		reference_kmers_len = ref->length;
		variant_kmers_len = var->length;
		reference_kmers = (uint64_t *) malloc(sizeof(uint64_t) * reference_kmers_len);
		variant_kmers = (uint64_t *) malloc(sizeof(uint64_t) * variant_kmers_len);
		memcpy(reference_kmers, ref->kmers, sizeof(uint64_t) * reference_kmers_len);
		memcpy(variant_kmers, var->kmers, sizeof(uint64_t) * variant_kmers_len);
		max_frequency = ref->max_frequency + var->max_frequency;
	}
	~VariantWindow() {
		free(reference_kmers);
		free(variant_kmers);
	}

	void ReverseKmers(uint8_t k) {	
		for (uint64_t i = 0; i < reference_kmers_len; i++)
			reference_kmers[i] = reverse_kmer(reference_kmers[i], k);
		for (uint64_t i = 0; i < variant_kmers_len; i++)
			variant_kmers[i] = reverse_kmer(variant_kmers[i], k);
	}

	void Print() {
		printf("Ref:");
		for (uint8_t i = 0; i < reference_kmers_len; i++) printf(" %lu", reference_kmers[i]);
		printf("\nVar:");
		for (uint8_t i = 0; i < variant_kmers_len; i++) printf(" %lu", variant_kmers[i]);
		printf("\nFreq: %u\n", max_frequency);
	}

	void Print(KmerFinder *kf) {
		printf("Ref:");
		for (uint8_t i = 0; i < reference_kmers_len; i++) {
			char *kmer = kf->graph->DecodeKmer(reference_kmers[i], kf->k);
			printf(" %s", kmer);
			free(kmer);
		}
		printf("\nVar:");
		for (uint8_t i = 0; i < variant_kmers_len; i++) {
			char *kmer = kf->graph->DecodeKmer(variant_kmers[i], kf->k);
			printf(" %s", kmer);
			free(kmer);
		}
		printf("\nFreq: %u\n", max_frequency);
	}
};

#endif
