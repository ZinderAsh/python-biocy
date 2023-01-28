#include "KmerFinder.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <stack>
#include <vector>
#include <algorithm>
#include <cmath>

#include "node.hpp"

void add_found(struct kmer_finder *kf, unsigned long node_id, unsigned long long kmer);
void get_kmers(struct kmer_finder *kf, unsigned long node_id);
void get_kmers_recursive(struct kmer_finder *kf, unsigned long node_id, unsigned char kmer_len, unsigned char kmer_ext_len);
unsigned long long full_mask = -1L;

KmerFinder::KmerFinder(Graph *graph, uint8_t k, uint8_t max_variant_nodes) : k(k), max_variant_nodes(max_variant_nodes) {
	this->graph = graph;
	kmer_mask = (1L << (k * 2)) - 1;
	kmer_buffer_shift = (33 - k) * 2;
	
	found_kmers = NULL;
	found_nodes = NULL;
	found_node_sequence_start_positions = NULL;
	found_node_sequence_kmer_positions = NULL;
	save_sequence_start_positions = false;
	save_sequence_kmer_positions = false;
	
	path_buffer = (uint32_t *) malloc(sizeof(uint32_t) * k * 4);
	kmer_position_buffer = (uint16_t *) malloc(sizeof(uint16_t) * k * 4);
	
	filters = 0;
	filter_node_id = 0;
}

void KmerFinder::Reset() {
	path_buffer_len = 1;
	variant_counter = 0;
	if (found_kmers) free(found_kmers);
	if (found_nodes) free(found_nodes);
	if (found_node_sequence_start_positions) free(found_node_sequence_start_positions);
	if (found_node_sequence_kmer_positions) free(found_node_sequence_kmer_positions);
	found_count = 0;
	found_len = 0;
}

void KmerFinder::InitializeFoundArrays() {	
	// Start found arrays with node number of slots (they are resized automatically when necessary)
	// This may be changed to be the length of the reference genome.
	found_len = graph->nodes_len * 1;
	found_kmers = (uint64_t *) malloc(found_len * sizeof(uint64_t));
	found_nodes = (uint32_t *) malloc(found_len * sizeof(uint32_t));
	if (save_sequence_start_positions) {
		found_node_sequence_start_positions = (uint32_t *) malloc(found_len * sizeof(uint32_t));	
	}
	if (save_sequence_kmer_positions) {
		found_node_sequence_kmer_positions = (uint16_t *) malloc(found_len * sizeof(uint16_t));	
	}
}

void KmerFinder::FindKmersForVariant(uint32_t reference_node_id, uint32_t variant_node_id) {
	Reset();

	bool old_save_sequence_start_positions = save_sequence_start_positions;
	bool old_save_sequence_kmer_positions = save_sequence_kmer_positions;
	save_sequence_start_positions = true;
	save_sequence_kmer_positions = true;

	InitializeFoundArrays();

	FindKmersSpanningNode(reference_node_id);
	FindKmersSpanningNode(variant_node_id);

	save_sequence_start_positions = old_save_sequence_start_positions;
	save_sequence_kmer_positions = old_save_sequence_kmer_positions;
}

void KmerFinder::FindKmersSpanningNode(uint32_t center_node_id) {
	SetFilter(FILTER_NODE_ID, center_node_id);

	std::vector<uint32_t> visited;
	std::stack<uint32_t> stack;
	stack.push(center_node_id);

	while (!stack.empty()) {
		uint32_t node_id = stack.top();
		stack.pop();
		visited.push_back(node_id);
		uint64_t local_found_count = FindKmersFromNode(node_id);

		if (local_found_count == 0) continue;

		struct node *node = graph->Get(node_id);
		for (uint8_t i = 0; i < node->edges_in_len; i++) {
			uint32_t edge_id = node->edges_in[i];
			if (std::find(visited.begin(), visited.end(), edge_id) == visited.end()) {
				stack.push(edge_id);
			}
		}
	}

	RemoveFilter(FILTER_NODE_ID);
}

void KmerFinder::Find() {
	Reset();

	InitializeFoundArrays();

	for (uint32_t i = 0; i < graph->nodes_len; i++) {
		if (graph->nodes[i].length != 0) FindKmersFromNode(i);
		if (0 && i % 100000 == 0) {
			printf("Progress: %d / %d nodes\n", i, graph->nodes_len);
		}
	}
}

bool KmerFinder::AddNodeFoundKmer(uint32_t node_id, uint64_t kmer, uint32_t start_position, uint16_t node_kmer_position) {
	// Check filters
	if (filters & FILTER_NODE_ID && node_id != filter_node_id) return false;
	// Resize found arrays if they are full
	if (found_count == found_len) {
		found_len *= 2;
		// printf("Attempting to allocate %lld slots for results\n", found_len);
		found_kmers = (uint64_t *) realloc(found_kmers, found_len * sizeof(uint64_t));
		found_nodes = (uint32_t *) realloc(found_nodes, found_len * sizeof(uint32_t));
		if (save_sequence_start_positions) {
			found_node_sequence_start_positions =
				(uint32_t *) realloc(found_node_sequence_start_positions, found_len * sizeof(uint32_t));	
		}
		if (save_sequence_kmer_positions) {
			found_node_sequence_kmer_positions =
				(uint16_t *) realloc(found_node_sequence_kmer_positions, found_len * sizeof(uint16_t));	
		}
		if (found_kmers == NULL || found_nodes == NULL) {
			printf("Failed to reallocate result arrays\n");
			exit(1);
		}
	}
	found_kmers[found_count] = kmer;
	found_nodes[found_count] = node_id;
	if (save_sequence_start_positions) found_node_sequence_start_positions[found_count] = start_position;
	if (save_sequence_kmer_positions) found_node_sequence_kmer_positions[found_count] = node_kmer_position;
	found_count++;
	return true;
}

void KmerFinder::ReverseFoundKmers() {
	for (uint64_t i = 0; i < found_count; i++) {
		uint64_t reverse = 0;
		uint64_t kmer = found_kmers[i];
		for (uint8_t j = 0; j < k; j++) {
			reverse |= ((kmer >> ((k - j - 1) * 2)) & 3L) << (j * 2);
		}
		found_kmers[i] = reverse;
	}
}

uint64_t KmerFinder::FindKmersFromNode(uint32_t node_id) {
	struct node *node = graph->nodes + node_id;

	uint64_t local_found_count = 0;

	// Ensure max variant nodes
	if (!node->reference) {
		if (variant_counter >= max_variant_nodes) return local_found_count;
		variant_counter++;
	}

	// Update lengths of buffers
	path_buffer[0] = node_id;
	kmer_position_buffer[0] = 0;

	// Store the node's sequence in the buffer
	kmer_buffer = node->sequences[0];
	uint32_t node_len = node->length;
	uint8_t kmer_len = (k < node_len) ? k : node_len;
	uint32_t sequence_idx = 1;
	uint8_t sequence_pos = 0;
	start_position = 0;

	// If at least k bases are stored, iterate and index kmers
	if (kmer_len == k) {
		kmer_len--;
		while (kmer_len < node_len) {
			kmer_len++;
			local_found_count += AddNodeFoundKmer(
					node_id,
					(kmer_buffer >> (64 - kmer_len * 2)) & kmer_mask,
					start_position,
					0);
			start_position++;
			
			// If the buffer is full, shift values as much as possible and fill with new values
			if (kmer_len == 32 && sequence_idx < node->sequences_len) {
				kmer_buffer <<= kmer_buffer_shift;
				kmer_buffer |= ((node->sequences[sequence_idx] << sequence_pos) >> (64 - kmer_buffer_shift));
				sequence_pos += kmer_buffer_shift;
				
				if (sequence_pos > 64) {
					sequence_idx++;
					if (sequence_idx < node->sequences_len) {
						sequence_pos -= 64;
						kmer_buffer |= (node->sequences[sequence_idx] >> (64 - sequence_pos));
					}
				} else if (sequence_pos == 64) {
					sequence_idx++;
					sequence_pos = 0;
				}

				kmer_len = k - 1;
				node_len -= 33 - k;
			}
		}
	}

	// Right-align buffer for easier bit-wise operations with kmer_buffer_ext in recursion
	kmer_buffer >>= 64 - kmer_len * 2;

	// Fake the kmer being max k - 1 long as no more bases are relevant during recursion
	kmer_len = (kmer_len < k - 1) ? kmer_len : (k - 1);

	// Visit all edges, remembering the current kmer buffer length
	for (uint8_t i = 0; i < node->edges_len; i++) {
		local_found_count += FindKmersExtendedByEdge(node->edges[i], kmer_len, 0);
	}

	// Count down variant nodes when done with this node
	if (!node->reference) variant_counter--;

	return local_found_count;
}

uint64_t KmerFinder::FindKmersExtendedByEdge(uint32_t node_id, uint8_t kmer_len, uint8_t kmer_ext_len) {
	struct node *node = graph->nodes + node_id;

	uint64_t local_found_count = 0;

	// Ensure max_variant_nodes is adhered to
	if (!node->reference) {
		if (variant_counter >= max_variant_nodes) return local_found_count;
		variant_counter++;
	}

	// Update lengths and buffers
	path_buffer[path_buffer_len] = node_id;
	kmer_position_buffer[path_buffer_len] = kmer_len + kmer_ext_len;
	path_buffer_len++;

	if (node->sequences_len != 0) {
		if (kmer_ext_len == 0) {
			// First recursion, replace entire buffer
			kmer_buffer_ext = node->sequences[0];
		} else {
			// Update the buffer with the node currently being visited
			kmer_buffer_ext &= (full_mask << (64 - kmer_ext_len * 2));
			kmer_buffer_ext |= (node->sequences[0] >> (kmer_ext_len * 2));
		}
		uint64_t kmer_hash;
		uint32_t node_len = kmer_ext_len + node->length;
		if (node_len > k - 1) node_len = k - 1;
		// Check if there are enough bases to index new kmers
		if (kmer_len + node_len >= k) {
			if (kmer_ext_len < k - kmer_len - 1) kmer_ext_len = k - kmer_len - 1;
			while (kmer_ext_len < node_len) {
				kmer_ext_len++;
				// Combine the main buffer and extended buffer to form the final kmer
				kmer_hash = ((kmer_buffer << kmer_ext_len * 2) |
						(kmer_buffer_ext >> (64 - kmer_ext_len * 2))) & kmer_mask;
				// Add the kmer to the found array for every node in the path
				for (uint16_t i = 0; i < path_buffer_len; i++) {
					uint16_t kmer_position = kmer_position_buffer[i];
					if (i > 0 && kmer_ext_len > k - kmer_len) {
						kmer_position -= kmer_ext_len + kmer_len - k;
					}
					uint32_t start_pos = 0;
					if (i == 0) {
						start_pos = start_position;
						if (kmer_ext_len > k - kmer_len) {
							start_pos += kmer_ext_len + kmer_len - k;
						}
					}
					local_found_count += AddNodeFoundKmer(
							path_buffer[i],
							kmer_hash,
							start_pos,
							kmer_position);
				}
			}
		} else {
			kmer_ext_len = node_len;
		}
	}

	if (kmer_ext_len < k - 1) {
		for (uint8_t i = 0; i < node->edges_len; i++) {
			local_found_count += FindKmersExtendedByEdge(node->edges[i], kmer_len, kmer_ext_len);
		}
	}

	path_buffer_len--;
	if (!node->reference) variant_counter--;

	return local_found_count;
}

std::unordered_map<uint64_t, uint32_t> KmerFinder::CreateKmerIndex() {
	std::unordered_map<uint64_t, uint32_t> index;
	int64_t reserved_buckets = found_count / (128 / k);
	//printf("Reserved: %ld\n", reserved_buckets);
	index.reserve(reserved_buckets);

	for (uint64_t i = 0; i < found_count; i++) {
		index[found_kmers[i]] += 1;
	}

	return index;
}

uint64_t unique_kmer_list(uint64_t *kmers, uint64_t length) {
	uint64_t *unique_kmers = (uint64_t *) malloc(sizeof(uint64_t) * length / 8);
	uint64_t unique_kmers_len = 0;
	for (uint64_t i = 0; i < length; i++) {
		uint64_t kmer = kmers[i];
		bool found = false;
		for (uint64_t j = 0; j < unique_kmers_len; j++) {
			if (unique_kmers[j] == kmer) {
				found = true;
				break;
			}
		}
		if (!found) {
			unique_kmers[unique_kmers_len++] = kmer;
		}
	}
	free(unique_kmers);
	return unique_kmers_len;
}

/*
int main(int argc, char** argv) {
	if (argc != 2) {
		printf("This program requires one argument: A filename for an npz file.\n");
	}
	Graph *graph = Graph::FromGFAFile(argv[1]);
	printf("%d nodes\n", graph->nodes_len);	

	graph->Compress();
	printf("%d nodes\n", graph->nodes_len);

	KmerFinder *kf = new KmerFinder(graph, 31, 31);
	kf->Find();
	printf("Done. %ld kmers found\n", kf->found_count);
	
	auto kmer_index = kf->CreateKmerIndex();
	printf("Kmer Count Size: %ld\n", kmer_index.size());
	printf("Kmer Count Bucket Count: %ld\n", kmer_index.bucket_count());
	
	delete kf;

	//uint32_t empty_node_count = graph->AddEmptyNodes();
	//printf("Added %u empty nodes.\n", empty_node_count);

	//empty_node_count = graph->AddEmptyNodes();
	//printf("Added another %u empty nodes.\n", empty_node_count);
	
	delete graph;
	return 0;
	graph->ToFile(argv[2]);
	
	kf = new KmerFinder(graph, 31, 31);
	kf->Find();
	printf("Done. %ld kmers found\n", kf->found_count);
	delete kf;
	
	delete graph;
	
	return 0;
}
*/
