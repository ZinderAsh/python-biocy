#include "Graph.hpp"

#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <queue>
#include <vector>
#include <bits/stdc++.h>

#include "GFA.hpp"
#include "VCF.hpp"
#include "FASTA.hpp"

#define LINE_BUF_LEN 1024
#define DEFAULT_ENCODING "ACGT"

struct queue_node {
	uint32_t id;
	uint32_t depth;

	queue_node(uint32_t id, uint32_t depth)
		: id(id), depth(depth) {}
};

Graph *Graph::FromGFAFile(char *filepath) {
	return FromGFAFileEncoded(filepath, DEFAULT_ENCODING);
}

Graph *Graph::FromGFAFileEncoded(char *filepath, const char *encoding) {
	GFA *gfa = GFA::ReadFile(filepath, encoding);
	
	Graph *graph = new Graph(encoding);
	
	graph->nodes_len = gfa->node_count;
	graph->nodes = (struct node *) malloc(sizeof(struct node) * gfa->node_count);
	memset(graph->nodes, 0, sizeof(struct node) * gfa->node_count);

	for (uint32_t index = 0; index < graph->nodes_len; index++) {
		struct node *node = (graph->nodes + index);
		node->length = gfa->sequence_lengths[index];
		if (node->length > 0) {
			node->sequences = (uint64_t *) malloc(sizeof(uint64_t));
			node->sequences[0] = gfa->sequences[index];
			node->sequences_len = 1;
		} else {
			node->sequences = 0;
			node->sequences_len = 0;
		}
		node->edges_len = gfa->edges_out_lengths[index];
		node->edges = (uint32_t *) malloc(sizeof(uint32_t) * node->edges_len);
		for (uint8_t i = 0; i < node->edges_len; i++) {
			node->edges[i] = gfa->edges_out[index][i];
		}
		node->edges_in_len = gfa->edges_in_lengths[index];
		node->edges_in = (uint32_t *) malloc(sizeof(uint32_t) * node->edges_in_len);
		for (uint8_t i = 0; i < node->edges_in_len; i++) {
			node->edges_in[i] = gfa->edges_in[index][i];
		}
		node->reference_index = gfa->reference_indices[index];
		node->reference = gfa->reference_nodes[index];
	}

	delete gfa;
	
	return graph;
}

Graph *Graph::FromFastaVCF(char *fasta_filepath, char *vcf_filepath, int16_t chromosome) {
	return FromFastaVCFEncoded(fasta_filepath, vcf_filepath, chromosome, DEFAULT_ENCODING);
}

Graph *Graph::FromFastaVCFEncoded(char *fasta_filepath, char *vcf_filepath, int16_t chromosome, const char *encoding) {
	VCF *vcf = VCF::ReadFile(vcf_filepath, chromosome);
	FASTA *fasta = FASTA::ReadFile(fasta_filepath);
	fasta->GoToChromosome(chromosome);

	Graph *graph = new Graph(encoding);

	/*
	uint64_t i = 99025;
	for (i = 0; i < 100000; i++) {
		if (strlen(vcf->references[i]) != strlen(vcf->variants[i]) && strstr(vcf->variants[i], ",")) {
			printf("%lu %d,%lu: %s / %s\n", i, vcf->chromosomes[i], vcf->positions[i], vcf->references[i], vcf->variants[i]);
		}
	}*/
	/*
	for (uint64_t i = 0; i < vcf->length; i++) {
		uint8_t ref_len = strlen(vcf->references[i]);
		for (uint64_t j = 0; j < vcf->length; j++) {
			if (i == j) continue;
			if (vcf->positions[i] > vcf->positions[j]) continue;
			uint64_t position_diff = vcf->positions[j] - vcf->positions[i];
			if (position_diff < ref_len) {
				printf("1: %lu %d,%lu: %s / %s\n", i,
						vcf->chromosomes[i], vcf->positions[i],
						vcf->references[i], vcf->variants[i]);
				printf("2: %lu %d,%lu: %s / %s\n", j,
						vcf->chromosomes[j], vcf->positions[j],
						vcf->references[j], vcf->variants[j]);
			}
		}
	}
	*/

	uint64_t *sorted_variant_indices = (uint64_t *) malloc(sizeof(uint64_t) * vcf->length);
	for (uint64_t i = 0; i < vcf->length; i++) {
		sorted_variant_indices[i] = i;
	}
	std::sort(sorted_variant_indices, sorted_variant_indices + vcf->length,
			[&vcf](const uint64_t a, const uint64_t b) -> bool {
				return vcf->positions[a] < vcf->positions[b];
			});

	uint64_t reference_pos = 0;
	uint32_t graph_previous_reference_id = 0;
	uint64_t variant_idx = 0;
	uint64_t reference_index = 0;

	uint32_t *variant_to_reference_node_id = (uint32_t *) malloc(sizeof(uint32_t) * vcf->length);
	uint32_t *variant_to_variant_node_id = (uint32_t *) malloc(sizeof(uint32_t) * vcf->length);
	memset(variant_to_reference_node_id, 0, sizeof(uint32_t) * vcf->length);
	memset(variant_to_variant_node_id, 0, sizeof(uint32_t) * vcf->length);

	uint32_t previous_variant_ids[128];
	uint8_t previous_variant_ids_len = 0;
	uint32_t next_variant_ids[128];
	uint8_t next_variant_ids_len = 0;

	uint32_t variants_added = 0;
	uint32_t variants_skipped_overlap = 0;
	while (variant_idx < vcf->length) {
		uint32_t real_variant_idx = sorted_variant_indices[variant_idx];
		uint64_t variant_pos = vcf->positions[real_variant_idx];
		variant_idx++;
		if (variant_idx % 100000 == 0 || variant_idx + 1 == vcf->length) printf("%lu / %lu variants processed\n", variant_idx, vcf->length);
		if (variant_pos < reference_pos) {
			//printf("Variant overlap\n");
			variants_skipped_overlap++;
		} else {
			variants_added++;
			// Read bases and add a reference node leading up to the variant
			uint64_t to_read = variant_pos - reference_pos;
			if (to_read > 0) {
				char *sequence = fasta->ReadNext(to_read);
				if (sequence == NULL) {
					printf("No more bases in FASTA (1)\n");
					break;
				}
				uint32_t new_reference_id = graph->AddNode(sequence);
				(graph->nodes + new_reference_id)->reference = true;
				(graph->nodes + new_reference_id)->reference_index = reference_index++;
				if (graph->nodes_len > 1) {
					graph->AddEdge(graph_previous_reference_id, new_reference_id);
				}
				for (uint8_t i = 0; i < previous_variant_ids_len; i++) {
					graph->AddEdge(previous_variant_ids[i], new_reference_id);
				}
				previous_variant_ids_len = 0;
				graph_previous_reference_id = new_reference_id;
				reference_pos += to_read;
			}
			
			// Read the reference bases of the variant and add the reference node
			to_read = strlen(vcf->references[real_variant_idx]);
			uint32_t variant_reference_id;
			if (to_read > 0) {
				char *sequence = fasta->ReadNext(to_read);
				if (sequence == NULL) {
					printf("No more bases in FASTA (2)\n");
					break;
				}
				for (uint64_t i = 0; i < to_read; i++) {
					if ((vcf->references[real_variant_idx][i] | 0x20) != (sequence[i] | 0x20)) {
						printf("Reference sequence mismatch! %s != %s\n", vcf->references[real_variant_idx], sequence);
						break;
					}
				}
				variant_reference_id = graph->AddNode(sequence);
			} else { // Empty node
				variant_reference_id = graph->AppendEmptyNode();
			}
			(graph->nodes + variant_reference_id)->reference = true;
			(graph->nodes + variant_reference_id)->reference_index = reference_index++;
			graph->AddEdge(graph_previous_reference_id, variant_reference_id);
			for (uint8_t i = 0; i < previous_variant_ids_len; i++) {
				graph->AddEdge(previous_variant_ids[i], variant_reference_id);
			}

			// Create variant nodes
			char *variant_string = strdup(vcf->variants[real_variant_idx]);
			char *tok = strtok(variant_string, ",");
			while (tok != NULL) {
				uint32_t tok_len = strlen(tok);
				uint32_t variant_node_id;
				if (tok_len > 0) {
					variant_node_id = graph->AddNode(tok);
				} else {
					variant_node_id = graph->AppendEmptyNode();
				}
				graph->AddEdge(graph_previous_reference_id, variant_node_id);
				for (uint8_t i = 0; i < previous_variant_ids_len; i++) {
					graph->AddEdge(previous_variant_ids[i], variant_node_id);
				}
				next_variant_ids[next_variant_ids_len++] = variant_node_id;
				tok = strtok(NULL, ",");
			}
			free(variant_string);

			memcpy(previous_variant_ids, next_variant_ids, sizeof(uint32_t) * next_variant_ids_len);
			previous_variant_ids_len = next_variant_ids_len;
			next_variant_ids_len = 0;

			graph_previous_reference_id = variant_reference_id;
			reference_pos += to_read;
		}
	}

	uint64_t remaining = 0;
	while (true) {
		char *a = fasta->ReadNext(128);
		if (a == NULL) break;
		remaining += strlen(a);
	}
	
	if (remaining > 0) {
		fasta->GoToChromosome(chromosome);
		uint64_t to_read = reference_pos;
		while (to_read > 128) {
			fasta->ReadNext(128);
			to_read -= 128;
		}
		fasta->ReadNext(to_read);
		char *sequence = fasta->ReadNext(remaining);
		uint32_t final_reference_id = graph->AddNode(sequence);
		if (graph->nodes_len > 1) {
			graph->AddEdge(graph_previous_reference_id, final_reference_id);
		}
		for (uint8_t i = 0; i < previous_variant_ids_len; i++) {
			graph->AddEdge(previous_variant_ids[i], final_reference_id);
		}
	}

	printf("Graph has %u nodes\n", graph->nodes_len);
	printf("Variants in graph: %u\n", variants_added);
	printf("Variants skipped due to overlap: %u\n", variants_skipped_overlap);

	free(sorted_variant_indices);
	free(variant_to_reference_node_id);
	free(variant_to_variant_node_id);

	delete vcf;
	delete fasta;

	return graph;
}

uint32_t *Graph::CreateCompressedIDMap(uint32_t *return_compressed_node_count) {
	bool visited[nodes_len];
	memset(visited, false, sizeof(bool) * nodes_len);

	uint32_t *id_map = (uint32_t *) malloc(sizeof(uint32_t) * nodes_len);
	memset(id_map, 0, sizeof(uint32_t) * nodes_len);

	uint32_t compressed_node_count = 0;
	for (uint32_t node_id = 0; node_id < nodes_len; node_id++) {
		if (visited[node_id]) continue;
		visited[node_id] = true;
		id_map[node_id] = compressed_node_count;
		uint32_t edge_id = node_id;
		struct node *edge = (nodes + edge_id);
		while (edge->edges_len == 1) {
			edge_id = edge->edges[0];
			edge = (nodes + edge_id);
			if (visited[edge_id]) break;
			if (edge->edges_in_len > 1) break;
			id_map[edge_id] = compressed_node_count;
			visited[edge_id] = true;
		}
		compressed_node_count++;
	}
	
	std::cout << "Compressed Node Count: " << compressed_node_count << std::endl;
	(*return_compressed_node_count) = compressed_node_count;
	return id_map;
}

void Graph::Compress() {
	for (uint32_t node_id = 0; node_id < nodes_len; node_id++) {
		if ((nodes + node_id)->sequences_len > 1) {
			std::cout << "This graph is already compressed." << std::endl;
			return;
		}
	}

	uint32_t compressed_node_count;
	uint32_t *id_map = CreateCompressedIDMap(&compressed_node_count);

	if (compressed_node_count == nodes_len) {
		std::cout << "This graph has no nodes that can be compressed." << std::endl;
		free(id_map);
		return;
	}

	printf("Optimizing from %u nodes to %u nodes\n", nodes_len, compressed_node_count);

	struct node *compressed_nodes = (struct node *) malloc(sizeof(struct node) * compressed_node_count);
	memset(compressed_nodes, 0, sizeof(struct node) * compressed_node_count);

	bool visited[nodes_len];
	memset(visited, false, sizeof(bool) * nodes_len);

	compressed_node_count = 0;
	for (uint32_t node_id = 0; node_id < nodes_len; node_id++) {
		if (visited[node_id]) continue;
		visited[node_id] = true;
		struct node *node = (nodes + node_id);
		struct node *compressed_node = (compressed_nodes + compressed_node_count);
		uint32_t node_length = node->length;
		uint32_t edge_id = node_id;
		struct node *edge = node;
		while (edge->edges_len == 1) {
			edge_id = edge->edges[0];
			edge = (nodes + edge_id);
			if (visited[edge_id] /*|| !(edge->reference)*/) break;
			if (edge->edges_in_len > 1) break;
			node_length += edge->length;
		}
		if (node_length == 0) {
			compressed_node->length = 0;
			compressed_node->sequences = NULL;
			compressed_node->sequences_len = 0;
		} else {
			compressed_node->length = node_length;
			compressed_node->sequences_len = (31 + node_length) / 32;
			compressed_node->sequences = (uint64_t *) malloc(sizeof(uint64_t) * compressed_node->sequences_len);
			memset(compressed_node->sequences, 0, sizeof(uint64_t) * compressed_node->sequences_len);
			compressed_node->sequences[0] = node->sequences[0];
			node_length = node->length;
			edge_id = node_id;
			edge = node;
			while (edge->edges_len == 1) {
				uint32_t temp_edge_id = edge->edges[0];
				struct node *temp_edge = (nodes + temp_edge_id);
				if (visited[temp_edge_id] /*|| !(temp_edge->reference)*/) break;
				if (temp_edge->edges_in_len > 1) break;
				edge_id = temp_edge_id;
				edge = temp_edge;
				if (node_length % 32 == 0) {
					compressed_node->sequences[node_length / 32] = edge->sequences[0];
					node_length += edge->length;
				} else {
					compressed_node->sequences[node_length / 32] |= (edge->sequences[0] >> ((node_length % 32) * 2));
					uint32_t new_length = node_length + edge->length;
					if (node_length / 32 < new_length / 32 && new_length % 32 != 0)
						compressed_node->sequences[node_length / 32 + 1] = (edge->sequences[0] << ((31 - (node_length % 32)) * 2));
					node_length = new_length;
				}
				visited[edge_id] = true;
			}
		}
		compressed_node->edges_len = edge->edges_len;
		if (compressed_node->edges_len == 0) {
			compressed_node->edges = NULL;
		} else {
			compressed_node->edges = (uint32_t *) malloc(sizeof(uint32_t) * edge->edges_len);
			for (uint8_t i = 0; i < compressed_node->edges_len; i++) {
				compressed_node->edges[i] = id_map[edge->edges[i]];
			}
		}
		compressed_node->edges_in_len = node->edges_in_len;
		if (compressed_node->edges_in_len == 0) {
			compressed_node->edges_in = NULL;
		} else {
			compressed_node->edges_in = (uint32_t *) malloc(sizeof(uint32_t) * node->edges_in_len);
			for (uint8_t i = 0; i < node->edges_in_len; i++) {
				compressed_node->edges_in[i] = id_map[node->edges_in[i]];
			}
		}
		compressed_node->reference_index = node->reference_index;
		compressed_node->reference = node->reference;

		compressed_node_count++;
	}

	for (uint32_t node_id = 0; node_id < nodes_len; node_id++) {
		struct node *node = (nodes + node_id);
		if (node->edges) free(node->edges);
		if (node->edges_in) free(node->edges_in);
		if (node->sequences) free(node->sequences);
	}
	free(nodes);
	nodes = compressed_nodes;
	nodes_len = compressed_node_count;

	uint32_t new_reference_index = 1;
	uint32_t first_ref_node_id = GetReferenceNodeID(0);
	struct node *ref_node = (nodes + first_ref_node_id);

	while (true) {
		uint32_t next_id = 0;
		int64_t next_reference_index = -1;
		for (uint8_t i = 0; i < ref_node->edges_len; i++) {
			struct node *edge = (nodes + (ref_node->edges[i]));
			if (edge->reference && (next_reference_index == -1 || edge->reference_index < next_reference_index)) {
				next_id = ref_node->edges[i];
				next_reference_index = edge->reference_index;
			}
		}
		if (next_reference_index == -1) break;
		ref_node = (nodes + next_id);
		ref_node->reference_index = new_reference_index++;
	}

	free(id_map);
}

uint32_t Graph::AddNode(const char *sequence) {
	nodes_len++;
	nodes = (struct node *) realloc(nodes, sizeof(struct node) * nodes_len);
	struct node *new_node = (nodes + nodes_len - 1);
	new_node->length = strlen(sequence);
	new_node->sequences_len = (new_node->length + 31) / 32;
	if (new_node->sequences_len > 0) {
		new_node->sequences = (uint64_t *) malloc(sizeof(uint64_t) * new_node->sequences_len);
		for (uint32_t i = 0; i < new_node->sequences_len; i++) {
			uint8_t hash_len = new_node->length - i * 32;
			if (hash_len > 32) hash_len = 32;
			const char *offset = sequence + (i * 32);
			uint64_t hash = hash_max_kmer_by_map(offset, hash_len, encoding_map);
			new_node->sequences[i] = hash;
		}
	} else {
		new_node->sequences = NULL;
	}
	new_node->edges = NULL;
	new_node->edges_in = NULL;
	new_node->edges_len = 0;
	new_node->edges_in_len = 0;
	new_node->reference_index = 0;
	new_node->reference = false;

	return nodes_len - 1;
}

void Graph::AddEdge(uint32_t from_node_id, uint32_t to_node_id) {
	struct node *from_node = (nodes + from_node_id);
	struct node *to_node = (nodes + to_node_id);

	from_node->edges_len++;
	from_node->edges = (uint32_t *) realloc(from_node->edges, sizeof(uint32_t) * from_node->edges_len);
	from_node->edges[from_node->edges_len - 1] = to_node_id;

	to_node->edges_in_len++;
	to_node->edges_in = (uint32_t *) realloc(to_node->edges_in, sizeof(uint32_t) * to_node->edges_in_len);
	to_node->edges_in[to_node->edges_in_len - 1] = from_node_id;
}

uint32_t Graph::GetRequiredEmptyNodesFromNode(uint32_t from_node_id) {
	struct node *from_node = (nodes + from_node_id);
	
	uint16_t nodes_required = 0;

	for (uint8_t i = 0; i < from_node->edges_len; i++) {
		struct node *edge = (nodes + from_node->edges[i]);
		if (edge->edges_in_len > 1) {
			//uint16_t nodes_required = GetRequiredEmptyNodesBetween(from_node_id, from_node->edges[i]);
			nodes_required++;
		}
	}
	
	return nodes_required;
}

uint32_t Graph::GetRootNodeID() {
	for (uint32_t i = 0; i < nodes_len; i++) {
		if ((nodes + i)->edges_in_len == 0)
			return i;
	}
	std::cout << "FATAL: Did not find a root node for the graph." << std::endl;
	return 0;
}

uint32_t Graph::GetReferenceNodeID(uint32_t reference_index) {
	for (uint32_t i = 0; i < nodes_len; i++) {
		if ((nodes + i)->reference_index == 0)
			return i;
	}
	std::cout << "FATAL: Did not find a reference node with index " << reference_index << " for the graph." << std::endl;
	return 0;
}

uint32_t Graph::GetLastNodeID() {
	for (uint32_t i = 0; i < nodes_len; i++) {
		if ((nodes + i)->edges_len == 0)
			return i;
	}
	std::cout << "FATAL: Did not find an end node for the graph." << std::endl;
	return 0;
}

uint32_t Graph::GetNextReferenceNodeID(uint32_t previous_id) {
	struct node *node = (nodes + previous_id);
	if (!(node->reference)) {
		std::cout << "FATAL: Can't get next reference node of non-reference node." << std::endl;
		return 0;
	}
	uint32_t next_reference_index = node->reference_index + 1;
	for (uint8_t i = 0; i < node->edges_len; i++) {
		struct node *edge = (nodes + node->edges[i]);
		if (edge->reference && edge->reference_index == next_reference_index) {
			return node->edges[i];
		}
	}
	std::cout << "FATAL: Did not find a next reference node." << std::endl;
	return 0;
}

uint32_t Graph::GetRequiredEmptyNodeCount() {
	uint32_t min_node_depth[nodes_len];
	uint32_t max_node_depth[nodes_len];
	for (uint32_t i = 0; i < nodes_len; i++) {
		min_node_depth[i] = nodes_len;
		max_node_depth[i] = 0;
	}

	uint32_t root_id = GetRootNodeID();
	min_node_depth[root_id] = 0;
	max_node_depth[root_id] = 0;

	uint32_t stack[nodes_len];
	uint32_t stack_len = 1;
	stack[0] = root_id;
	
	while (stack_len > 0) {
		uint32_t node_id = stack[--stack_len];
		struct node *node = (nodes + node_id);
		uint32_t min_depth = min_node_depth[node_id] + 1;
		uint32_t max_depth = max_node_depth[node_id] + 1;
		if (max_depth > nodes_len) {
			std::cout << "ERROR: This graph has a cycle. Cannot add empty nodes." << std::endl;
			return 0;
		}
		for (uint8_t i = 0; i < node->edges_len; i++) {
			uint32_t edge_id = node->edges[i];
			bool updated = false;
			if (min_depth < min_node_depth[edge_id]) {
				min_node_depth[edge_id] = min_depth;
				updated = true;
			}
			if (max_depth > max_node_depth[edge_id]) {
				max_node_depth[edge_id] = max_depth;
				updated = true;
			}
			if (updated) stack[stack_len++] = edge_id;
		}
	}

	uint32_t last_node = GetLastNodeID();
	return max_node_depth[last_node] - min_node_depth[last_node];
}

bool Graph::NodeHasEdge(uint32_t node_id, uint32_t edge_id) {
	struct node *node = (nodes + node_id);

	for (uint8_t edge_idx = 0; edge_idx < node->edges_len; edge_idx++) {
		if (node->edges[edge_idx] == edge_id)
			return true;
	}
	return false;
}

void Graph::InitializeEmptyNode(uint32_t node_id) {
	struct node *node = (nodes + node_id);
	node->length = 0;
	node->sequences_len = 0;
	node->sequences = NULL;
	node->edges = NULL;
	node->edges_len = 0;
	node->edges_in = NULL;
	node->edges_in_len = 0;
	node->reference_index = 0;
	node->reference = false;
}

uint32_t Graph::CreateEmptyNodes() {
	uint32_t node_depth[nodes_len];
	uint32_t last_update_parent[nodes_len];
	bool in_queue[nodes_len];
	memset(node_depth, 0, sizeof(uint32_t) * nodes_len);
	memset(last_update_parent, 0, sizeof(uint32_t) * nodes_len);
	memset(in_queue, false, sizeof(bool) * nodes_len);

	uint32_t root_id = GetRootNodeID();
	uint32_t nodes_created = 0;

	uint32_t more_than_1 = 0;

	auto cmp = [&](uint32_t l, uint32_t r) { return node_depth[l] > node_depth[r]; };
	std::priority_queue<uint32_t, std::vector<uint32_t>, decltype(cmp)> minq(cmp);
	std::vector<uint32_t> tmpq;

	minq.push(root_id);
	in_queue[root_id] = true;
	
	while (!minq.empty()) {
		uint32_t node_id = minq.top();
		minq.pop();
		in_queue[node_id] = false;
		uint32_t depth = node_depth[node_id];
		//std::cout << "Node: " << node_id << ", Depth: " << depth << std::endl;
		struct node *node = (nodes + node_id);
		if (depth > nodes_len) {
			std::cout << "ERROR: This graph has a cycle. Cannot add empty nodes." << std::endl;
			return 0;
		}
		for (uint8_t i = 0; i < node->edges_len; i++) {
			uint32_t edge_id = node->edges[i];
			bool updated = false;
			uint32_t edge_depth = depth + 1;
			if (node_depth[edge_id] == 0) {
				// Depth has not been set yet
				updated = true;
			} else if (node_depth[edge_id] != edge_depth) {
				int32_t diff = node_depth[edge_id] - edge_depth;
				if (diff > 1 || diff < -1) more_than_1++;
				if (node_depth[edge_id] < edge_depth) {
					// Depth mismatch, need empty node between
					struct node *edge = (nodes + edge_id);
					bool found_empty_node = false;
					for (uint8_t j = 0; j < edge->edges_in_len; j++) {
						uint32_t in_edge_id = edge->edges_in[j];
						struct node *in_edge = (nodes + in_edge_id);
						if (in_edge->length == 0 && node_depth[in_edge_id] == depth + 1) {
							found_empty_node = true;
							break;
						}
					}
					if (!found_empty_node) {
						
					}
					
					updated = true;
				} else {
					std::cout << "Test" << std::endl;
				}
			}
			if (updated) {
				last_update_parent[edge_id] = node_id;
				if (in_queue[edge_id]) {
					while (true) {
						uint32_t tmp_node = minq.top();
						minq.pop();
						tmpq.push_back(tmp_node);
						if (tmp_node == edge_id) break;
					}
					node_depth[edge_id] = edge_depth;
					while (tmpq.size() > 0) {
						minq.push(tmpq[tmpq.size() - 1]);
						tmpq.pop_back();
					}
				} else {
					node_depth[edge_id] = edge_depth;
					minq.push(edge_id);
					in_queue[edge_id] = true;
				}
			}
		}
	}

	std::cout << more_than_1 << " cases of depth diff > 1" << std::endl;

	return nodes_created;
}

/*
bool Graph::ConnectToEmptyNode(uint32_t node_id, uint32_t edge_id) {
	struct node *node = (nodes + node_id);
	struct node *edge = (nodes + edge_id);
	
	uint32_t empty_node_id = 0;
	struct node *empty_node = NULL;

	// Check for an existing empty node
	FOR (uint8_t edge_idx = 0; edge_idx < edge->in_edges_len; edge_idx++) {
		if ((nodes + (edge->edges_in[edge_idx]))->length == 0) {
			empty_node_id = edge->edges_in[edge_idx];
			empty_node = (nodes + empty_node_id);
			empty_node->edges = (uint32_t *) realloc(edges, sizeof(uint32_t) * (empty_node->edges_len++));
			break;
		}
	}
}

uint32_t Graph::AddEmptyNodeBetween(uint32_t node_id, uint32_t edge_id) {
	struct node *node = (nodes + node_id);
	struct node *edge = (nodes + edge_id);
	
	uint32_t empty_node_id = 0;
	struct node *empty_node = NULL;

	// Check for an existing empty node
	for (uint8_t edge_idx = 0; edge_idx < edge->in_edges_len; edge_idx++) {
		if ((nodes + (edge->edges_in[edge_idx]))->length == 0) {
			empty_node_id = edge->edges_in[edge_idx];
			empty_node = (nodes + empty_node_id);
			empty_node->edges = (uint32_t *) realloc(edges, sizeof(uint32_t) * (empty_node->edges_len++));
			break;
		}
	}

	// Create a new empty node
	empty_node_id = nodes_len++;
	empty_node = (nodes + empty_node_id);

	InitializeEmptyNode(empty_node_id);

	empty_node->edges[0] = edge_id;
	empty_node->edges_in[0] = node_id;

	for (uint8_t edge_idx = 0; edge_idx < node->edges_len; edge_idx++) {
		if (node->edges[edge_idx] == edge_id) {
			node->edges[edge_idx] = empty_node_idx;
			break;
		}
	}
	for (uint8_t edge_idx = 0; edge_idx < edge->edges_in_len; edge_idx++) {
		if (node->edges_in[edge_idx] == node_id) {
			node->edges_in[edge_idx] = empty_node_idx;
			break;
		}
	}
}

uint32_t Graph::AddEmptyNodesForNode(uint32_t node_id) {
	struct node *node = (nodes + node_id);

	if (node->edges_len <= 1) return;

	uint32_t nodes_added = 0;

	for (uint8_t edge_idx = 0; edge_idx < node->edges_len; edge_idx++) {
		struct node *edge = (nodes + node->edges[edge_idx]);
		
		if (edge->edges_in_len <= 1) continue;

		for (uint8_t in_edge_idx = 0; in_edge_idx < edge->edges_in_len; in_edge_idx++) {
			if (NodeHasEdge(node_id, edge->edges_in[in_edge_idx])) {
				nodes_added += AddEmptyNodeBewteen(node_id, node->edges[edge_idx]);
			}
		}
	}

	return nodes_added;
}
*/

uint32_t Graph::AppendEmptyNode() {
	uint32_t new_node_id = nodes_len++;
	nodes = (struct node *) realloc(nodes, sizeof(struct node) * nodes_len);
	InitializeEmptyNode(new_node_id);
	return new_node_id;
}

uint32_t Graph::AddEmptyNodes() {
	/*
	uint32_t empty_node_count = GetRequiredEmptyNodeCount();

	if (empty_node_count == 0) {
		std::cout << "There are no empty nodes to add to the graph." << std::endl;
		return 0;
	}

	printf("Adding %u empty nodes\n", empty_node_count);

	nodes = (struct node *) realloc(nodes, sizeof(struct node) * (nodes_len + empty_node_count));
	
	uint32_t real_added = CreateEmptyNodes();

	printf("Actually added %u nodes\n", real_added);
	*/

	uint32_t empty_node_count = 0;

	for (uint32_t node_id = 0; node_id < nodes_len; node_id++) {
		struct node *node = (nodes + node_id);
		if (!(node->reference)) continue;
		for (uint8_t i = 0; i < node->edges_len; i++) {
			uint32_t edge_id = node->edges[i];
			struct node *edge = (nodes + edge_id);
			if (edge->reference && edge->reference_index == node->reference_index + 2) {
				uint32_t empty_node_id = AppendEmptyNode();
				node = (nodes + node_id);
				empty_node_count++;
				MoveEdgesToIntermediateNode(node_id, edge_id, empty_node_id);
			}
		}
	}

	return empty_node_count;
}

void Graph::MoveEdgesToIntermediateNode(uint32_t from_node_id, uint32_t to_node_id, uint32_t mid_node_id) {
	struct node *from_node = (nodes + from_node_id);
	struct node *to_node = (nodes + to_node_id);
	struct node *mid_node = (nodes + mid_node_id);
	
	int16_t from_node_edge_index = -1;
	int16_t to_node_edge_index = -1;

	for (uint8_t i = 0; i < from_node->edges_len; i++) {
		if (from_node->edges[i] == to_node_id) {
			from_node_edge_index = i;
			break;
		}
	}

	if (from_node_edge_index == -1) {
		std::cout << "FATAL: Failed to find edge to the specified node." << std::endl;
		return;
	}

	for (uint8_t i = 0; i < to_node->edges_in_len; i++) {
		if (to_node->edges_in[i] == from_node_id) {
			to_node_edge_index = i;
			break;
		}
	}

	if (to_node_edge_index == -1) {
		std::cout << "FATAL: Failed to find edge from the specified node." << std::endl;
		return;
	}

	from_node->edges[from_node_edge_index] = mid_node_id;
	to_node->edges_in[to_node_edge_index] = mid_node_id;

	mid_node->edges_len++;
	mid_node->edges_in_len++;

	if (mid_node->edges == NULL) {
		mid_node->edges = (uint32_t *) malloc(sizeof(uint32_t) * mid_node->edges_len);
	} else {
		mid_node->edges = (uint32_t *) realloc(mid_node->edges, sizeof(uint32_t) * mid_node->edges_len);
	}
	mid_node->edges[mid_node->edges_len - 1] = to_node_id;

	if (mid_node->edges_in == NULL) {
		mid_node->edges_in = (uint32_t *) malloc(sizeof(uint32_t) * mid_node->edges_in_len);
	} else {
		mid_node->edges_in = (uint32_t *) realloc(mid_node->edges_in, sizeof(uint32_t) * mid_node->edges_in_len);
	}
	mid_node->edges_in[mid_node->edges_in_len - 1] = from_node_id;
}

Graph *Graph::FromFile(char *filepath) {
	FILE *f = fopen(filepath, "rb");
	int format_code_len = strlen(BCG_FORMAT_CODE);

	for (int i = 0; i < format_code_len; i++) {
		if (fgetc(f) != BCG_FORMAT_CODE[i]) {
			return NULL;
		}
	}

	uint8_t version_number = fgetc(f);

	char encoding[5];
	encoding[4] = '\0';

	for (uint8_t i = 0; i < 4; i++) {
		encoding[i] = fgetc(f);
	}
	
	Graph *graph = new Graph(encoding);

	fread(&(graph->nodes_len), sizeof(uint32_t), 1, f);
	graph->nodes = (struct node *) malloc(sizeof(struct node) * graph->nodes_len);
	memset(graph->nodes, 0, sizeof(struct node) * graph->nodes_len);

	for (uint32_t i = 0; i < graph->nodes_len; i++) {
		struct node *n = (graph->nodes + i);
		fread(&(n->length), sizeof(uint32_t), 1, f);
		fread(&(n->sequences_len), sizeof(uint32_t), 1, f);
		n->sequences = (uint64_t *) malloc(sizeof(uint64_t) * n->sequences_len);
		memset(n->sequences, 0, sizeof(uint64_t) * n->sequences_len);
		fread(n->sequences, sizeof(uint64_t), n->sequences_len, f);
		uint8_t edges_len = 0;
		fread(&edges_len, sizeof(uint8_t), 1, f);
		n->edges_len = edges_len;
		n->edges = (uint32_t *) malloc(sizeof(uint32_t) * edges_len);
		memset(n->edges, 0, sizeof(uint32_t) * edges_len);
		fread(n->edges, sizeof(uint32_t), edges_len, f);
		uint8_t reference = 0;
		fread(&reference, sizeof(uint8_t), 1, f);
		n->reference = reference;
	}

	fclose(f);

	return graph;
}

void Graph::ToFile(char *filepath) {
	FILE *f = fopen(filepath, "wb");

	// Indicator for file format
	fwrite(BCG_FORMAT_CODE, sizeof(char), strlen(BCG_FORMAT_CODE), f);

	// Format version number
	fputc(1, f);

	// Encoding
	fwrite(encoding, sizeof(char), 4, f);

	// Write graph info (currently only number of nodes)
	fwrite(&(nodes_len), sizeof(uint32_t), 1, f);

	// Write all nodes
	for (uint32_t i = 0; i < nodes_len; i++) {
		struct node *n = (nodes + i);
		fwrite(&(n->length), sizeof(uint32_t), 1, f);
		fwrite(&(n->sequences_len), sizeof(uint32_t), 1, f);
		fwrite(n->sequences, sizeof(uint64_t), n->sequences_len, f);
		uint8_t edges_len = n->edges_len;
		fwrite(&edges_len, sizeof(uint8_t), 1, f);
		fwrite(n->edges, sizeof(uint32_t), n->edges_len, f);
		uint8_t reference = n->reference;
		fwrite(&reference, sizeof(uint8_t), 1, f);
	}

	fclose(f);
}

void Graph::SetEncoding(const char *encoding) {
	memcpy(this->encoding, encoding, sizeof(char) * 4);
	fill_map_by_encoding(this->encoding_map, encoding);
}

