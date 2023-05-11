#ifndef KIVS_FILE_READER_GFA
#define KIVS_FILE_READER_GFA

#include <cstdlib>
#include <stdio.h>
#include <stdint.h>
#include <cstring>
#include "hashing.hpp"

class GFA {
public:
	char encoding[4];
	uint32_t node_count;
	uint64_t *sequences;
	uint32_t *sequence_lengths;
	uint32_t **edges_out;
	uint32_t **edges_in;
	uint8_t *edges_out_lengths;
	uint8_t *edges_in_lengths;
	uint32_t *reference_indices;
	bool *reference_nodes;

private:
	uint8_t encoding_map[256];
	char *filepath;
	FILE *source_file;
	uint32_t min_node_id;
	uint32_t max_node_id;
	uint32_t *id_map;
	bool use_id_map;
	uint32_t id_offset;
	char line_buffer[1024];
	char id_buffer[64];
	char edge_buffer[64];
	char sequence_buffer[64];
	uint8_t id_buffer_length;

public:
	~GFA() {
		if (sequences) free(sequences);
		if (sequence_lengths) free(sequence_lengths);
		if (edges_out) {
			for (uint32_t i = 0; i < node_count; i++) {
				if (edges_out[i]) free(edges_out[i]);
			}
			free(edges_out);
		}
		if (edges_in) {
			for (uint32_t i = 0; i < node_count; i++) {
				if (edges_in[i]) free(edges_in[i]);
			}
			free(edges_in);
		}
		if (edges_out_lengths) free(edges_out_lengths);
		if (edges_in_lengths) free(edges_in_lengths);
		if (reference_indices) free(reference_indices);
		if (reference_nodes) free(reference_nodes);
		if (id_map) free(id_map);
		free(filepath);
	}

	static GFA *ReadFile(char *filepath, const char *encoding);
private:

	GFA(char *filepath, const char *encoding) {
		this->filepath = strdup(filepath);
		memcpy(this->encoding, encoding, sizeof(char) * 4);
		fill_map_by_encoding(this->encoding_map, encoding);
		sequences = NULL;
		sequence_lengths = NULL;
		edges_out = NULL;
		edges_in = NULL;
		edges_out_lengths = NULL;
		edges_in_lengths = NULL;
		reference_indices = NULL;
		reference_nodes = NULL;
		id_map = NULL;
	}

	void ReadNodeCountAndIDRange();
	void InitializeArrays();
	void ReadSequences();
	void ReadReferenceNodes();
	void CountEdges();
	void InitializeEdges();
	void ReadEdges();
	int64_t ReadNextID();
	uint32_t GetMappedID(uint32_t);
};

#endif
