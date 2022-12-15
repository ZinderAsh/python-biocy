#include "Graph.hpp"

#include <stdlib.h>
#include <string.h>
#include <iostream>

#include "GFA.hpp"

#define LINE_BUF_LEN 1024
#define DEFAULT_ENCODING "ACGT"

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
		node->reference = gfa->reference_nodes[index];
	}

	delete gfa;
	
	return graph;
}

uint32_t *Graph::CreateCompressedIDMap(uint32_t *return_compressed_node_count) {
	bool visited[nodes_len];
	memset(visited, false, sizeof(bool) * nodes_len);

	uint32_t *id_map = (uint32_t *) malloc(sizeof(uint32_t) * nodes_len);
	memset(id_map, 0, sizeof(uint32_t) * nodes_len);

	uint32_t compressed_node_count = 0;
	for (uint32_t node_id = 0; node_id < nodes_len; node_id++) {
		struct node *node = (nodes + node_id);
		if (visited[node_id]) continue;
		visited[node_id] = true;
		id_map[node_id] = compressed_node_count;
		if (node->reference) {
			uint32_t edge_id = node_id;
			struct node *edge = (nodes + edge_id);
			while (edge->edges_len == 1) {
				edge_id = edge->edges[0];
				edge = (nodes + edge_id);
				if (visited[edge_id] || !(edge->reference)) break;
				if (edge->edges_in_len > 1) break;
				id_map[edge_id] = compressed_node_count;
				visited[edge_id] = true;
			}
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
		if (node->reference) {
			uint32_t node_length = node->length;
			uint32_t edge_id = node_id;
			struct node *edge = node;
			while (edge->edges_len == 1) {
				edge_id = edge->edges[0];
				edge = (nodes + edge_id);
				if (visited[edge_id] || !(edge->reference)) break;
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
					if (visited[temp_edge_id] || !(temp_edge->reference)) break;
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
		} else {
			compressed_node->length = node->length;
			if (compressed_node->length > 0) {
				compressed_node->sequences = (uint64_t *) malloc(sizeof(uint64_t));
				compressed_node->sequences[0] = node->sequences[0];
				compressed_node->sequences_len = 1;
			} else {
				compressed_node->sequences = 0;
				compressed_node->sequences_len = 0;
			}
			compressed_node->edges_len = node->edges_len;
			if (compressed_node->edges_len == 0) {
				compressed_node->edges = NULL;
			} else {
				compressed_node->edges = (uint32_t *) malloc(sizeof(uint32_t) * node->edges_len);
				for (uint8_t i = 0; i < node->edges_len; i++) {
					compressed_node->edges[i] = id_map[node->edges[i]];
				}
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

	free(id_map);
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
	printf("File version: %d\n", version_number);

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

