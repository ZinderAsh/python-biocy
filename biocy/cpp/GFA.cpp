#include "GFA.hpp"
#include <iostream>

GFA *GFA::ReadFile(char *filepath, const char *encoding) {
	GFA *gfa = new GFA(filepath, encoding);
	gfa->source_file = fopen(filepath, "rb");

	gfa->ReadNodeCountAndIDRange();

	gfa->InitializeArrays();

	gfa->ReadSequences();

	gfa->ReadReferenceNodes();

	gfa->CountEdges();

	gfa->InitializeEdges();

	gfa->ReadEdges();
	
	fclose(gfa->source_file);

	return gfa;
}

void GFA::InitializeArrays() {
	sequences          = (uint64_t *) malloc(sizeof(uint64_t) * node_count);
	sequence_lengths   = (uint32_t *) malloc(sizeof(uint32_t) * node_count);
	edges_out          = (uint32_t **) malloc(sizeof(uint32_t *) * node_count);
	edges_in           = (uint32_t **) malloc(sizeof(uint32_t *) * node_count);
	edges_out_lengths  = (uint8_t *) malloc(sizeof(uint8_t) * node_count);
	edges_in_lengths   = (uint8_t *) malloc(sizeof(uint8_t) * node_count);
	reference_indices  = (uint32_t *) malloc(sizeof(uint32_t) * node_count);
	reference_nodes    = (bool *) malloc(sizeof(bool) * node_count);
	id_map             = (uint32_t *) malloc(sizeof(uint32_t) * node_count);

	memset(reference_indices, 0, sizeof(uint32_t) * node_count);
	memset(reference_nodes, false, sizeof(uint8_t) * node_count);
	memset(edges_in_lengths, 0, sizeof(uint8_t) * node_count);
	memset(edges_out_lengths, 0, sizeof(uint8_t) * node_count);
}

void GFA::ReadNodeCountAndIDRange() {
	bool newline = false;
	int c;
	node_count = 0;
	min_node_id = -1L;
	max_node_id = 0;
	
	fseek(source_file, 0, SEEK_SET);
	while ((c = fgetc(source_file)) != EOF) {
		if (c == 'S' && newline) {
			fgetc(source_file);
			fgets(line_buffer, 1024, source_file);
			sscanf(line_buffer, "%s\t", id_buffer);
			uint32_t id = strtoull(id_buffer, NULL, 10);
			node_count++;
			if (id > max_node_id) max_node_id = id;
			if (id < min_node_id) min_node_id = id;
		} else {
			newline = (c == '\n');
		}
	}

	if (min_node_id + node_count - 1 == max_node_id) {
		printf("Node ID space is continuous from ID %u to %u.\n", min_node_id, max_node_id);
		if (min_node_id != 0) printf("Remapping IDs to start from 0.\n");
		id_offset = min_node_id;
		use_id_map = false;
	} else {
		printf("Node ID space is not continuous. Remapping IDs completely.\n");
		use_id_map = true;
	}
}

void GFA::ReadSequences() {
	bool newline = false;
	int c;
	int32_t id_map_index = 0;
	fseek(source_file, 0, SEEK_SET);
	while ((c = fgetc(source_file)) != EOF) {
		if (c == 'S' && newline) {
			fgetc(source_file);
			fgets(line_buffer, 1024, source_file);
			sscanf(line_buffer, "%s\t%s\n", id_buffer, sequence_buffer);
			uint32_t id = strtoull(id_buffer, NULL, 10);
			if (use_id_map) {
				id_map[id_map_index] = id;
				id = id_map_index;
				id_map_index++;
			} else {
				id -= id_offset;
			}
			sequence_lengths[id] = strlen(sequence_buffer);
			sequences[id] = hash_max_kmer_by_map(sequence_buffer, sequence_lengths[id], encoding_map);
		} else {
			newline = (c == '\n');
		}
	}
}

void GFA::ReadReferenceNodes() {
	bool newline = false;
	int c;
	fseek(source_file, 0, SEEK_SET);
	while ((c = fgetc(source_file)) != EOF) {
		if (c == 'P' && newline) {
			fgetc(source_file);
			while ((c = fgetc(source_file)) != '\t') {}
			int64_t id;
			uint32_t ref_index = 0;
			while ((id = ReadNextID()) != -1) {
				reference_nodes[id] = true;
				reference_indices[id] = ref_index++;
			}
		} else {
			newline = (c == '\n');
		}
	}
}

void GFA::CountEdges() {
	bool newline = false;
	int c;
	fseek(source_file, 0, SEEK_SET);
	while ((c = fgetc(source_file)) != EOF) {
		if (c == 'L' && newline) {
			fgetc(source_file);
			fgets(line_buffer, 1024, source_file);
			sscanf(line_buffer, "%s\t+\t%s\t+\n", id_buffer, edge_buffer);
			uint32_t id = strtoull(id_buffer, NULL, 10);
			uint32_t edge = strtoull(edge_buffer, NULL, 10);
			id = (use_id_map ? GetMappedID(id) : (id - id_offset));
			edge = (use_id_map ? GetMappedID(edge) : (edge - id_offset));
			edges_out_lengths[id]++;
			edges_in_lengths[edge]++;
		} else {
			newline = (c == '\n');
		}
	}
}

void GFA::InitializeEdges() {
	for (uint32_t index = 0; index < node_count; index++) {
		if (edges_out_lengths[index] > 0) {
			edges_out[index] = (uint32_t *) malloc(sizeof(uint32_t) * edges_out_lengths[index]);
		} else {
			edges_out[index] = NULL;
		}

		if (edges_in_lengths[index] > 0) {
			edges_in[index] = (uint32_t *) malloc(sizeof(uint32_t) * edges_in_lengths[index]);
		} else {
			edges_in[index] = NULL;
		}
	}
}

void GFA::ReadEdges() {
	memset(edges_out_lengths, 0, sizeof(uint8_t) * node_count);
	memset(edges_in_lengths, 0, sizeof(uint8_t) * node_count);

	bool newline;
	int c;
	fseek(source_file, 0, SEEK_SET);
	while ((c = fgetc(source_file)) != EOF) {
		if (c == 'L' && newline) {
			fgetc(source_file);
			fgets(line_buffer, 1024, source_file);
			sscanf(line_buffer, "%s\t+\t%s\t+\n", id_buffer, edge_buffer);
			uint32_t id = strtoull(id_buffer, NULL, 10);
			uint32_t edge = strtoull(edge_buffer, NULL, 10);
			id = (use_id_map ? GetMappedID(id) : (id - id_offset));
			edge = (use_id_map ? GetMappedID(edge) : (edge - id_offset));
			edges_out[id][edges_out_lengths[id]++] = edge;
			edges_in[edge][edges_in_lengths[edge]++] = id;
		} else {
			newline = (c == '\n');
		}
	}
}

int64_t GFA::ReadNextID() {
	int c;
	id_buffer_length = 0;
	while ((c = fgetc(source_file)) != '\n') {
		if (c >= '0' && c <= '9') {
			id_buffer[id_buffer_length++] = c;
		} else if (id_buffer_length > 0) {
			id_buffer[id_buffer_length] = '\0';
			int64_t id = strtoll(id_buffer, NULL, 10);
			id = (use_id_map ? GetMappedID(id) : (id - id_offset));
			return id;
		}
	}
	return -1;
}

uint32_t GFA::GetMappedID(uint32_t id) {
	for (uint32_t mapped_id = 0; mapped_id < node_count; mapped_id++) {
		if (id_map[mapped_id] == id) {
			return mapped_id;
		}
	}
	std::cout << "Fatal Error: Did not find a mapped ID." << std::endl;
	return -1;
}


