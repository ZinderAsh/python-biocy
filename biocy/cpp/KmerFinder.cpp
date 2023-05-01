#include "KmerFinder.hpp"

#include <stdlib.h>
#include <stdio.h>
#include <stack>
#include <algorithm>
#include <cmath>

uint64_t full_mask = -1L;

KmerFinder::KmerFinder(Graph *graph, uint8_t k, uint8_t max_variant_nodes) : k(k), max_variant_nodes(max_variant_nodes) {
	this->graph = graph;
	kmer_mask = (1L << (k * 2)) - 1;
	kmer_buffer_shift = (33 - k) * 2;
	
	found_kmers = NULL;
	found_nodes = NULL;
	found_node_sequence_start_positions = NULL;
	found_node_sequence_kmer_positions = NULL;
	found_windows = NULL;
	save_sequence_start_positions = false;
	save_sequence_kmer_positions = false;
	
	path_buffer = (uint32_t *) malloc(sizeof(uint32_t) * k * 4);
	kmer_position_buffer = (uint16_t *) malloc(sizeof(uint16_t) * k * 4);
	
	flags = 0;
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
	if (found_windows) {
		for (uint32_t i = 0; i < found_window_count; i++) {
			free((found_windows + i)->kmers);
		}
		free(found_windows);
	};
	found_count = 0;
	found_len = 0;
}

void KmerFinder::InitializeFoundArrays() {	
	// Start found arrays with node number of slots (they are resized automatically when necessary)
	// This may be changed to be the length of the reference genome.
	if (flags & FLAG_TO_STDOUT) return;
	if (flags & FLAG_SAVE_WINDOWS) {
		found_window_len = k * 2;
		found_window_count = 0;
		found_windows = (struct kmer_window *) malloc(sizeof(struct kmer_window) * found_window_len);
	} else {
		found_len = graph->nodes_len * 1;
		found_count = 0;
		found_kmers = (uint64_t *) malloc(found_len * sizeof(uint64_t));
		found_nodes = (uint32_t *) malloc(found_len * sizeof(uint32_t));
		if (save_sequence_start_positions) {
			found_node_sequence_start_positions = (uint32_t *) malloc(found_len * sizeof(uint32_t));
		}
		if (save_sequence_kmer_positions) {
			found_node_sequence_kmer_positions = (uint16_t *) malloc(found_len * sizeof(uint16_t));
		}
	}
}

uint32_t KmerFinder::GetWindowOverlap(std::vector<VariantWindow *> *windows, uint32_t window_index) {
	uint32_t overlap = 0;
	VariantWindow *window = (*windows)[window_index];
	VariantWindow *other_window;
	uint64_t kmer;
	
	for (uint32_t i = 0; i < window->reference_kmers_len; i++) {
		kmer = window->reference_kmers[i];
		for (uint32_t w = 0; w < windows->size(); w++) {
			if (w == window_index) continue;
			other_window = (*windows)[w];
			for (uint32_t j = 0; j < other_window->variant_kmers_len; j++) {
				if (kmer == other_window->variant_kmers[j])
					overlap++;
			}
		}
	}
	
	for (uint32_t i = 0; i < window->variant_kmers_len; i++) {
		kmer = window->variant_kmers[i];
		for (uint32_t w = 0; w < windows->size(); w++) {
			if (w == window_index) continue;
			other_window = (*windows)[w];
			for (uint32_t j = 0; j < other_window->reference_kmers_len; j++) {
				if (kmer == other_window->reference_kmers[j])
					overlap++;
			}
		}
	}

	return overlap;
}

VariantWindow *KmerFinder::FindAlignedSignaturesWithFinder(uint32_t reference_node_id, uint32_t variant_node_id, KmerFinder *kf) {
	bool minimize_overlap = (flags & FLAG_MINIMIZE_SIGNATURE_OVERLAP);
	auto windows = FindWindowsForVariantWithFinder(reference_node_id, variant_node_id, kf);
	uint32_t min_frequency = windows[0]->max_frequency;
	uint32_t min_overlap = (minimize_overlap ? GetWindowOverlap(&windows, 0) : 0);
	uint32_t min_window_index = 0;
	for (uint32_t i = 1; i < windows.size(); i++) {
		uint32_t overlap = (minimize_overlap ? GetWindowOverlap(&windows, i) : 0);
		if (overlap < min_overlap) {
			min_overlap = overlap;
			min_frequency = windows[i]->max_frequency;
			min_window_index = i;
		} else if (overlap == min_overlap && windows[i]->max_frequency < min_frequency) {
			min_overlap = overlap;
			min_frequency = windows[i]->max_frequency;
			min_window_index = i;
		}
	}
	for (uint32_t i = 0; i < windows.size(); i++) {
		if (i != min_window_index) delete windows[i];
	}
	return windows[min_window_index];
}

VariantWindow *KmerFinder::FindUnalignedSignaturesWithFinder(uint32_t reference_node_id, uint32_t variant_node_id, KmerFinder *kf) {
	kf->FindKmersForVariant(reference_node_id, variant_node_id);

	struct kmer_window *found_window_end = (kf->found_windows + kf->found_window_count);
	struct kmer_window *variant_window_start = kf->found_windows;
	while (variant_window_start < found_window_end) {
		if (variant_window_start->node_id == variant_node_id) break;
		variant_window_start++;
	}

	if (variant_window_start == kf->found_windows) {
		printf("FATAL: Failed to find windows for the variant node.\n");
		return NULL;
	}

	struct kmer_window *ref_window, *var_window;
	struct kmer_window *best_ref_window = kf->found_windows;
	struct kmer_window *best_var_window = variant_window_start;
	uint32_t min_ref_overlap = 99999, min_var_overlap = 99999;

	for (ref_window = kf->found_windows; ref_window < variant_window_start; ref_window++) {
		uint32_t overlap = 0;
		if (flags & FLAG_MINIMIZE_SIGNATURE_OVERLAP) {
			for (var_window = variant_window_start; var_window < found_window_end; var_window++) {
				for (uint32_t i = 0; i < ref_window->length; i++) {
					for (uint32_t j = 0; j < var_window->length; j++) {
						if (ref_window->kmers[i] == var_window->kmers[j]) overlap++;
					}
				}
			}
		}
		if (overlap < min_ref_overlap) {
			best_ref_window = ref_window;
			min_ref_overlap = overlap;
		} else if (overlap == min_ref_overlap && ref_window->max_frequency < best_ref_window->max_frequency) {
			best_ref_window = ref_window;
		}
	}

	for (var_window = variant_window_start; var_window < found_window_end; var_window++) {
		uint32_t overlap = 0;
		if (flags & FLAG_MINIMIZE_SIGNATURE_OVERLAP) {
			for (ref_window = kf->found_windows; ref_window < variant_window_start; ref_window++) {
				for (uint32_t i = 0; i < var_window->length; i++) {
					for (uint32_t j = 0; j < ref_window->length; j++) {
						if (var_window->kmers[i] == ref_window->kmers[j]) overlap++;
					}
				}
			}
		}
		if (overlap < min_var_overlap) {
			best_var_window = var_window;
			min_var_overlap = overlap;
		} else if (overlap == min_var_overlap && var_window->max_frequency < best_var_window->max_frequency) {
			best_var_window = var_window;
		}
	}

	return new VariantWindow(best_ref_window, best_var_window);
}

VariantWindow *KmerFinder::FindVariantSignaturesWithFinder(uint32_t reference_node_id, uint32_t variant_node_id, KmerFinder *kf) {
	if (flags & FLAG_ALIGN_SIGNATURE_WINDOWS) {
		return FindAlignedSignaturesWithFinder(reference_node_id, variant_node_id, kf);
	} else {
		return FindUnalignedSignaturesWithFinder(reference_node_id, variant_node_id, kf);
	}
}

VariantWindow *KmerFinder::FindVariantSignatures(uint32_t reference_node_id, uint32_t variant_node_id) {
	KmerFinder *kf = CreateWindowFinder();

	VariantWindow *result = FindVariantSignaturesWithFinder(reference_node_id, variant_node_id, kf);

	delete kf;

	return result;
}

KmerFinder *KmerFinder::CreateWindowFinder() {
	KmerFinder *kf = new KmerFinder(graph, k, max_variant_nodes);
	if (kmer_frequency_index.empty()) {
		kmer_frequency_index = CreateKmerFrequencyIndex();
	}
	kf->SetKmerFrequencyIndex(kmer_frequency_index);
	kf->SetFlag(FLAG_SAVE_WINDOWS, true);
	return kf;
}

std::vector<VariantWindow *> KmerFinder::FindWindowsForVariantWithFinder(uint32_t reference_node_id, uint32_t variant_node_id, KmerFinder *kf) {
	kf->FindKmersForVariant(reference_node_id, variant_node_id);
	
	std::vector<VariantWindow *> windows;

	struct kmer_window *found_window_end = (kf->found_windows + kf->found_window_count);
	struct kmer_window *variant_window_start = kf->found_windows;
	while (variant_window_start < found_window_end) {
		if (variant_window_start->node_id == variant_node_id) break;
		variant_window_start++;
	}

	if (variant_window_start == kf->found_windows) {
		printf("FATAL: Failed to find windows for the variant node.\n");
		return windows;
	}

	uint32_t ref_length = graph->nodes[reference_node_id].length;
	uint32_t var_length = graph->nodes[variant_node_id].length;

	struct kmer_window *ref_window, *var_window;
	bool save_window;

	for (ref_window = kf->found_windows; ref_window < variant_window_start; ref_window++) {
		for (var_window = variant_window_start; var_window < found_window_end; var_window++) {
			save_window = false;
			int16_t ref_left_bases = ref_window->kmer_position - ref_window->start_position;
			int16_t var_left_bases = var_window->kmer_position - var_window->start_position;
			if (ref_left_bases == var_left_bases) {
				save_window = true;
			}
			if (!save_window) {
				int16_t ref_right_bases = k + ref_window->start_position - ref_length - ref_window->kmer_position;
				int16_t var_right_bases = k + var_window->start_position - var_length - var_window->kmer_position;
				if (ref_right_bases == var_right_bases) {
					save_window = true;
				}
			}
			if (save_window) {
				windows.push_back(new VariantWindow(ref_window, var_window));
			}
		}
	}

	return windows;
}

std::vector<VariantWindow *> KmerFinder::FindWindowsForVariant(uint32_t reference_node_id, uint32_t variant_node_id) {
	KmerFinder *kf = CreateWindowFinder();

	std::vector<VariantWindow *> result = FindWindowsForVariantWithFinder(reference_node_id, variant_node_id, kf);

	delete kf;

	return result;
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
	bool found_any = false;
	uint32_t visited_count = 0;

	while (!stack.empty()) {
		visited_count++;
		uint32_t node_id = stack.top();
		stack.pop();
		visited.push_back(node_id);
		struct node *node = graph->Get(node_id);

		if (node->length > 0) {
			uint64_t local_found_count = FindKmersFromNode(node_id);
			if (local_found_count == 0) {
				if (found_any) continue;
			} else
				found_any = true;
		}

		if (!found_any && visited_count > k * k) {
			printf("FATAL: Failed to find any kmers spanning node %u.\n", center_node_id);
			break;
		}

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

uint64_t KmerFinder::GetKmerFrequency(uint64_t kmer) {
	return kmer_frequency_index[kmer];
}

bool KmerFinder::AddFoundWindowKmer(uint32_t node_id, uint64_t kmer, uint32_t start_position, uint16_t node_kmer_position) {
	uint64_t kmer_frequency = GetKmerFrequency(kmer);
	for (uint64_t i = 0; i < found_window_count; i++) {
		struct kmer_window *window = found_windows + i;
		if (window->node_id == node_id &&
		    window->start_position == start_position &&
		    window->kmer_position == node_kmer_position) {
			window->length++;
			window->kmers = (uint64_t *) realloc(window->kmers, sizeof(uint64_t) * window->length);
			window->kmers[window->length - 1] = kmer;
			if (kmer_frequency > window->max_frequency) {
				window->max_frequency = kmer_frequency;
			}
			return true;
		}
	}
	if (found_window_count == found_window_len) {
		found_window_len = found_window_len * 3 / 2;
		found_windows = (struct kmer_window *) realloc(found_windows, sizeof(struct kmer_window) * found_window_len);
	}
	struct kmer_window *window = found_windows + found_window_count;
	window->node_id = node_id;
	window->length = 1;
	window->kmers = (uint64_t *) malloc(sizeof(uint64_t) * window->length);
	window->kmers[0] = kmer;
	window->start_position = start_position;
	window->kmer_position = node_kmer_position;
	window->max_frequency = kmer_frequency;
	
	found_window_count++;

	return true;
}

bool KmerFinder::AddNodeFoundKmer(uint32_t node_id, uint64_t kmer, uint32_t start_position, uint16_t node_kmer_position) {
	// Check filters
	if (filters & FILTER_NODE_ID && node_id != filter_node_id) return false;
	if (flags & FLAG_TO_STDOUT) {
		printf("%lu\t%u:%u:%u\n", kmer, node_id, start_position, node_kmer_position);
		return true;
	}
	if (flags & FLAG_SAVE_WINDOWS) {
		return AddFoundWindowKmer(node_id, kmer, start_position, node_kmer_position);
	}
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
		found_kmers[i] = reverse_kmer(found_kmers[i], k);
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
				uint16_t nodes_to_save = path_buffer_len;
				if (flags & FLAG_ONLY_SAVE_INITIAL_NODES) nodes_to_save = 1;
				// Add the kmer to the found array for every node in the path
				for (uint16_t i = 0; i < nodes_to_save; i++) {
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

std::unordered_map<uint64_t, uint32_t> KmerFinder::CreateKmerFrequencyIndex() {
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
