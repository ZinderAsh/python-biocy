#ifndef KMER_FINDER_H
#define KMER_FINDER_H

#include "graph.h"

struct kmer_finder {
	struct graph *graph;
	
	unsigned char k;
	unsigned char max_variant_nodes;
	
	unsigned long long *found_kmers;
	unsigned long *found_nodes;
	unsigned long long found_len;
	unsigned long long found_count;

	unsigned long long kmer_buffer;
	unsigned long long kmer_buffer_ext;

	unsigned long *path_buffer;
	unsigned char path_buffer_len;

	unsigned long long kmer_mask;
	unsigned char shift;

	unsigned char variant_counter;
};

struct kmer_finder *init_kmer_finder(struct graph *graph, unsigned char k, unsigned char max_variant_nodes);
void free_kmer_finder(struct kmer_finder *kf);
void find_kmers(struct kmer_finder *kf);
void reverse_kmer_endian(struct kmer_finder *kf);

#endif
