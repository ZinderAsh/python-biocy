#include "kmer_finder.h"
#include <stdlib.h>

void add_found(struct kmer_finder *kf, unsigned long node_id, unsigned long long kmer);
void get_kmers(struct kmer_finder *kf, unsigned long node_id);
void get_kmers_recursive(struct kmer_finder *kf, unsigned long node_id, unsigned char kmer_len, unsigned char kmer_ext_len);
unsigned long long full_mask = -1L;

struct kmer_finder *init_kmer_finder(struct graph *graph, unsigned char k, unsigned char max_variant_nodes) {
	struct kmer_finder *kf = (struct kmer_finder *) malloc(sizeof(struct kmer_finder));
	kf->graph = graph;
	kf->k = k;
	kf->max_variant_nodes = max_variant_nodes;
	kf->kmer_mask = (1L << (k * 2)) - 1;
	kf->shift = (32 - k) * 2;
	kf->found_kmers = NULL;
	kf->found_nodes = NULL;
	kf->path_buffer = (unsigned long *) malloc(k * 4 * sizeof(unsigned long));
	return kf;
}

void free_kmer_finder(struct kmer_finder *kf) {
	if (kf->found_kmers) free(kf->found_kmers);
	if (kf->found_nodes) free(kf->found_nodes);
	free(kf->path_buffer);
	free(kf);
}

void find_kmers(struct kmer_finder *kf) {
	kf->path_buffer_len = 1;
	kf->variant_counter = 0;
	if (kf->found_kmers) free(kf->found_kmers);
	if (kf->found_nodes) free(kf->found_nodes);

	// Start found arrays with node number of slots (they are resized automatically when necessary)
	// This may be changed to be the length of the reference genome.
	kf->found_len = kf->graph->nodes_len * 1000;
	kf->found_kmers = malloc(kf->found_len * sizeof(unsigned long long));
	kf->found_nodes = malloc(kf->found_len * sizeof(unsigned long));
	kf->found_count = 0;

	for (unsigned int i = 0; i < kf->graph->nodes_len; i++) {
		if (kf->graph->nodes[i].len != 0) get_kmers(kf, i);
	}
}

void add_found(struct kmer_finder *kf, unsigned long node_id, unsigned long long kmer) {
	// Resize found arrays if they are full
	if (kf->found_count == kf->found_len) {
		kf->found_len *= 2;
		kf->found_kmers = realloc(kf->found_kmers, kf->found_len * sizeof(unsigned long long));
		kf->found_nodes = realloc(kf->found_nodes, kf->found_len * sizeof(unsigned long));
	}
	kf->found_kmers[kf->found_count] = kmer;
	kf->found_nodes[kf->found_count] = node_id;
	kf->found_count++;
}

void get_kmers(struct kmer_finder *kf, unsigned long node_id) {
	// Update lengths of buffers
	kf->path_buffer[0] = node_id;
	struct node *node = kf->graph->nodes + node_id;

	// Ensure max variant nodes
	if (!node->reference) {
		if (kf->variant_counter >= kf->max_variant_nodes) return;
		kf->variant_counter++;
	}

	// Store the node's sequence in the buffer
	kf->kmer_buffer = node->sequences[0];
	unsigned int node_len = node->len;
	unsigned char kmer_len = (kf->k < node_len) ? kf->k : node_len;
	unsigned short sequence_idx = 1;
	unsigned char sequence_pos = 0;

	// If at least k bases are stored, iterate and index kmers
	if (kmer_len == kf->k) {
		kmer_len--;
		while (kmer_len < node_len) {
			kmer_len++;
			add_found(kf, node_id, (kf->kmer_buffer >> (64 - kmer_len * 2)) & kf->kmer_mask);
			
			// If the buffer is full, shift values as much as possible and fill with new values
			if (kmer_len == 31 && sequence_idx < node->sequences_len) {
				kf->kmer_buffer <<= kf->shift;
				kf->kmer_buffer |= ((node->sequences[sequence_idx] << sequence_pos) >> (62 - kf->shift));
				sequence_pos += kf->shift;
				
				if (sequence_pos > 62) {
					sequence_idx++;
					if (sequence_idx < node->sequences_len) {
						sequence_pos -= 62;
						kf->kmer_buffer |= (node->sequences[sequence_idx] >> (62 - sequence_pos));
					}
				} else if (sequence_pos == 62) {
					sequence_idx++;
					sequence_pos = 0;
				}

				kmer_len = kf->k - 1;
				node_len -= 32 - kf->k;
			}
		}
	}

	// Right-align buffer for easier bit-wise operations with kmer_buffer_ext in recursion
	kf->kmer_buffer >>= 64 - kmer_len * 2;

	// Fake the kmer being max k - 1 long as no more bases are relevant during recursion
	kmer_len = (kmer_len < kf->k - 1) ? kmer_len : (kf->k - 1);

	// Visit all edges, remembering the current kmer buffer length
	for (unsigned char i = 0; i < node->edges_len; i++) {
		get_kmers_recursive(kf, node->edges[i], kmer_len, 0);
	}

	// Count down variant nodes when done with this node
	if (!node->reference) kf->variant_counter--;
}

void get_kmers_recursive(struct kmer_finder *kf, unsigned long node_id, unsigned char kmer_len, unsigned char kmer_ext_len) {
	// Update lengths and buffers
	kf->path_buffer[kf->path_buffer_len] = node_id;
	kf->path_buffer_len++;
	struct node *node = kf->graph->nodes + node_id;

	// Ensure max_variant_nodes is adhered to
	if (!node->reference) {
		if (kf->variant_counter >= kf->max_variant_nodes) return;
		kf->variant_counter++;
	}

	if (node->sequences_len != 0) {
		if (kmer_ext_len == 0) {
			// First recursion, replace entire buffer
			kf->kmer_buffer_ext = node->sequences[0];
		} else {
			// Update the buffer with the node currently being visited
			kf->kmer_buffer_ext &= (full_mask << (64 - kmer_ext_len * 2));
			kf->kmer_buffer_ext |= (node->sequences[0] >> (kmer_ext_len * 2));
		}
		unsigned long long kmer_hash;
		unsigned int node_len = kmer_ext_len + node->len;
		if (node_len > kf->k - 1) node_len = kf->k - 1;
		// Check if there are enough bases to index new kmers
		if (kmer_len + node_len >= kf->k) {
			if (kmer_ext_len < kf->k - kmer_len - 1) kmer_ext_len = kf->k - kmer_len - 1;
			while (kmer_ext_len < node_len) {
				kmer_ext_len++;
				// Combine the main buffer and extended buffer to form the final kmer
				kmer_hash = ((kf->kmer_buffer << kmer_ext_len * 2) |
						(kf->kmer_buffer_ext >> (64 - kmer_ext_len * 2))) & kf->kmer_mask;
				// Add the kmer to the found array for every node in the path
				for (int i = 0; i < kf->path_buffer_len; i++) {
					add_found(kf, kf->path_buffer[i], kmer_hash);
				}
			}
		} else {
			kmer_ext_len = node_len;
		}
	}

	if (kmer_ext_len < kf->k - 1) {
		for (unsigned char i = 0; i < node->edges_len; i++) {
			get_kmers_recursive(kf, node->edges[i], kmer_len, kmer_ext_len);
		}
	}

	kf->path_buffer_len--;
	if (!node->reference) kf->variant_counter--;
}





