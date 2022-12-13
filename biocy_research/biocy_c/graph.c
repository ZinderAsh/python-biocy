#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "graph.h"
#include "utils.h"

#define LINE_BUF_LEN 1024
#define DEFAULT_ENCODING "ACGT"

void get_file_str_between(FILE *f, char *buf, long start_pos, char *from, char *to);

int from_file_gfa(char *filepath, struct graph *graph, char flags) {
	return from_file_gfa_encoded(filepath, graph, DEFAULT_ENCODING, flags);
}

int from_file_gfa_encoded(char *filepath, struct graph *graph, const char *encoding, char flags) {
	set_encoding(graph, encoding);
	return from_file_gfa_by_map(filepath, graph, graph->encoding_map, flags);
}

int from_file_gfa_by_map(char *filepath, struct graph *graph, unsigned char *map, char flags) {
	FILE *f = fopen(filepath, "rb");

	unsigned long node_count = 0;
	unsigned long min_node_id = -1L;
	unsigned long max_node_id = 0;

	char newline = 0;
	int c, i;
	char line_buf[LINE_BUF_LEN];
	char id_buf[64], edge_buf[64], sequence_buf[64];
	unsigned char id_buf_len;
	unsigned long id, edge;

	memset(id_buf, 0, sizeof(char) * 64);
	memset(edge_buf, 0, sizeof(char) * 64);
	memset(sequence_buf, 0, sizeof(char) * 64);

	// Find node ID space
	while ((c = fgetc(f)) != EOF) {
		if (c == 'S' && newline) {
			fgetc(f);
			fgets(line_buf, LINE_BUF_LEN, f);
			sscanf(line_buf, "%s\t", id_buf);
			id = strtoull(id_buf, NULL, 10);
			node_count++;
			if (id > max_node_id) max_node_id = id;
			if (id < min_node_id) min_node_id = id;
		} else {
			newline = (c == '\n');
		}
	}

	unsigned long long *sequences = malloc(sizeof(unsigned long long) * node_count);
	unsigned char *sequence_lens = malloc(sizeof(unsigned char) * node_count);
	unsigned char *ref_nodes = malloc(sizeof(unsigned char) * node_count);
	memset(ref_nodes, 0, sizeof(unsigned char) * node_count);
	unsigned long **edges_in = malloc(sizeof(unsigned long *) * node_count);
	unsigned char *edges_in_lens = malloc(sizeof(unsigned char) * node_count);
	unsigned long **edges_out = malloc(sizeof(unsigned long *) * node_count);
	unsigned char *edges_out_lens = malloc(sizeof(unsigned char) * node_count);
	memset(edges_in_lens, 0, sizeof(unsigned char) * node_count);
	memset(edges_out_lens, 0, sizeof(unsigned char) * node_count);
	unsigned long *id_map = NULL;
	
	if (min_node_id + node_count - 1 == max_node_id) {
		printf("Node ID space is continuous from ID %lu to %lu.\n", min_node_id, max_node_id);
		if (min_node_id != 0) printf("Remapping IDs to start from 0.\n");
	} else {
		printf("Node ID space is not continuous. Remapping IDs completely.\n");
		id_map = malloc(sizeof(unsigned long) * node_count);
	}

	// Get all sequences
	i = 0;
	fseek(f, 0, SEEK_SET);
	while ((c = fgetc(f)) != EOF) {
		if (c == 'S' && newline) {
			fgetc(f);
			fgets(line_buf, LINE_BUF_LEN, f);
			sscanf(line_buf, "%s\t%s\n", id_buf, sequence_buf);
			id = strtoul(id_buf, NULL, 10);
			if (id_map) {
				id_map[i] = id;
				sequence_lens[i] = strlen(sequence_buf);
				sequences[i] = hash_max_kmer_by_map(sequence_buf, sequence_lens[i], map);
				i++;
			} else {
				sequence_lens[id - min_node_id] = strlen(sequence_buf);
				sequences[id - min_node_id] = hash_max_kmer_by_map(sequence_buf, sequence_lens[id - min_node_id], map);
			}
		} else {
			newline = (c == '\n');
		}
	}

	// Find and label ref nodes
	fseek(f, 0, SEEK_SET);
	while ((c = fgetc(f)) != EOF) {
		if (c == 'P' && newline) {
			fgetc(f);
			while ((c = fgetc(f)) != '\t') {}
			id_buf_len = 0;
			while ((c = fgetc(f)) != '\n') {
				if (c >= '0' && c <= '9') {
					id_buf[id_buf_len++] = c;
				} else if (id_buf_len > 0) {
					id_buf[id_buf_len] = '\0';
					id = strtoul(id_buf, NULL, 10);
					if (id_map) {
						for (unsigned long index = 0; index < node_count; index++) {
							if (id_map[index] == id) {
								ref_nodes[index] = 1;
								break;
							}
						}
					} else {
						ref_nodes[id - min_node_id] = 1;
					}
					id_buf_len = 0;
				}
			}
		} else {
			newline = (c == '\n');
		}
	}	

	// Count edges to allocate arrays
	fseek(f, 0, SEEK_SET);
	while ((c = fgetc(f)) != EOF) {
		if (c == 'L' && newline) {
			fgetc(f);
			fgets(line_buf, LINE_BUF_LEN, f);
			sscanf(line_buf, "%s\t+\t%s\t+\n", id_buf, edge_buf);
			id = strtoul(id_buf, NULL, 10);
			edge = strtoul(edge_buf, NULL, 10);
			if (id_map) {
				char found = 0;
				for (unsigned long index = 0; index < node_count; index++) {
					if (id_map[index] == id && ((found & 1) == 0)) {
						edges_out_lens[index]++;
						found |= 1;
						if (found == 3) break;
					}
					if (id_map[index] == edge && ((found & 2) == 0)) {
						edges_in_lens[index]++;
						found |= 2;
						if (found == 3) break;
					}
				}
			} else {
				edges_in_lens[edge - min_node_id]++;
				edges_out_lens[id - min_node_id]++;
			}
		} else {
			newline = (c == '\n');
		}
	}

	// Allocating edge arrays
	for (unsigned long index = 0; index < node_count; index++) {
		if (edges_in_lens[index] > 0) {
			edges_in[index] = malloc(sizeof(unsigned long) * edges_in_lens[index]);
		} else {
			edges_in[index] = NULL;
		}

		if (edges_out_lens[index] > 0) {
			edges_out[index] = malloc(sizeof(unsigned long) * edges_out_lens[index]);
		} else {
			edges_out[index] = NULL;
		}
	}

	memset(edges_in_lens, 0, sizeof(unsigned char) * node_count);
	memset(edges_out_lens, 0, sizeof(unsigned char) * node_count);

	// Fill edge arrays
	fseek(f, 0, SEEK_SET);
	while ((c = fgetc(f)) != EOF) {
		if (c == 'L' && newline) {
			fgetc(f);
			fgets(line_buf, LINE_BUF_LEN, f);
			sscanf(line_buf, "%s\t+\t%s\t+\n", id_buf, edge_buf);
			id = strtoul(id_buf, NULL, 10);
			edge = strtoul(edge_buf, NULL, 10);
			if (id_map) {
				char found = 0;
				for (unsigned long index = 0; index < node_count; index++) {
					if (id_map[index] == id && ((found & 1) == 0)) {
						id = index;
						found |= 1;
						if (found == 3) break;
					}
					if (id_map[index] == id && ((found & 2) == 0)) {
						edge = index;
						found |= 2;
						if (found == 3) break;
					}
				}
				if (found == 3) {
					edges_in[edge][edges_in_lens[edge - min_node_id]++] = id - min_node_id;
					edges_out[id][edges_out_lens[id - min_node_id]++] = edge - min_node_id;
				}
			} else {
				edges_in[edge - min_node_id][edges_in_lens[edge - min_node_id]++] = id - min_node_id;
				edges_out[id - min_node_id][edges_out_lens[id - min_node_id]++] = edge - min_node_id;
			}
		} else {
			newline = (c == '\n');
		}
	}

	// Done reading file
	fclose(f);

	if (flags & GRAPH_OPTIMIZE) {

		char *visited = malloc(sizeof(char) * node_count);
		memset(visited, 0, sizeof(char) * node_count);

		if (!id_map) id_map = malloc(sizeof(unsigned long) * node_count);
		memset(id_map, 0, sizeof(unsigned long) * node_count);

		unsigned long optimized_node_count = 0;
		for (unsigned long index = 0; index < node_count; index++) {
			if (visited[index]) continue;
			visited[index] = 1;
			id_map[index] = optimized_node_count;
			if (ref_nodes[index]) {
				unsigned long edge = index;
				while (edges_out_lens[edge] == 1) {
					edge = edges_out[edge][0];
					if (visited[edge] || !ref_nodes[edge]) break;
					if (edges_in_lens[edge] > 1) break;
					id_map[edge] = optimized_node_count;
					visited[edge] = 1;
				}
			}
			optimized_node_count++;
		}

		printf("Optimizing from %lu nodes to %lu nodes\n", node_count, optimized_node_count);

		graph->nodes_len = optimized_node_count;
		graph->nodes = malloc(sizeof(struct node) * optimized_node_count);
		memset(graph->nodes, 0, sizeof(struct node) * optimized_node_count);
		memset(visited, 0, sizeof(char) * node_count);

		optimized_node_count = 0;
		for (unsigned long index = 0; index < node_count; index++) {
			if (visited[index]) continue;
			visited[index] = 1;
			struct node *node = (graph->nodes + optimized_node_count);
			if (ref_nodes[index]) {
				unsigned long node_length = sequence_lens[index];
				unsigned long edge = index;
				while (edges_out_lens[edge] == 1) {
					edge = edges_out[edge][0];
					if (visited[edge] || !ref_nodes[edge]) break;
					node_length += sequence_lens[edge];
				}
				if (node_length == 0) {
					node->length = 0;
					node->sequences = NULL;
					node->sequences_len = 0;
				} else {
					node->length = node_length;
					node->sequences_len = (31 + node_length) / 32;
					node->sequences = malloc(sizeof(unsigned long long) * node->sequences_len);
					memset(node->sequences, 0, sizeof(unsigned long long) * node->sequences_len);
					node->sequences[0] = sequences[index];
					node_length = sequence_lens[index];
					edge = index;
					while (edges_out_lens[edge] == 1) {
						if (visited[edges_out[edge][0]] || !ref_nodes[edges_out[edge][0]]) break;
						if (edges_in_lens[edges_out[edge][0]] > 1) break;
						edge = edges_out[edge][0];
						if (node_length % 32 == 0) {
							node->sequences[node_length / 32] = sequences[edge];
							node_length += sequence_lens[edge];
						} else {
							node->sequences[node_length / 32] |= (sequences[edge] >> ((node_length % 32) * 2));
							if (node_length / 32 < (node_length + sequence_lens[edge]) / 32 &&
									(node_length + sequence_lens[edge]) % 32 != 0)
								node->sequences[node_length / 32 + 1] = (sequences[edge] << ((31 - (node_length % 32)) * 2));
							node_length += sequence_lens[edge];
						}
						visited[edge] = 1;
					}
				}
				node->edges = malloc(sizeof(unsigned long) * edges_out_lens[edge]);
				node->edges_len = edges_out_lens[edge];
				for (int i = 0; i < node->edges_len; i++) {
					node->edges[i] = id_map[edges_out[edge][i]];
				}
				node->edges_in = malloc(sizeof(unsigned long) * edges_in_lens[index]);
				node->edges_in_len = edges_in_lens[index];
				for (int i = 0; i < node->edges_in_len; i++) {
					node->edges_in[i] = id_map[edges_in[index][i]];
				}
				node->reference = ref_nodes[index];
			} else {
				node->length = sequence_lens[index];
				if (node->length > 0) {
					node->sequences = malloc(sizeof(unsigned long long));
					node->sequences[0] = sequences[index];
					node->sequences_len = 1;
				} else {
					node->sequences = 0;
					node->sequences_len = 0;
				}
				node->edges = malloc(sizeof(unsigned long) * edges_out_lens[index]);
				node->edges_len = edges_out_lens[index];
				for (int i = 0; i < node->edges_len; i++) {
					node->edges[i] = id_map[edges_out[index][i]];
				}
				node->edges_in = malloc(sizeof(unsigned long) * edges_in_lens[index]);
				node->edges_in_len = edges_in_lens[index];
				for (int i = 0; i < node->edges_in_len; i++) {
					node->edges_in[i] = id_map[edges_in[index][i]];
				}
				node->reference = ref_nodes[index];
			}
			optimized_node_count++;
		}
		free(visited);
	} else { // Do not optimize
		graph->nodes_len = node_count;
		graph->nodes = malloc(sizeof(struct node) * node_count);
		memset(graph->nodes, 0, sizeof(struct node) * node_count);

		for (unsigned long index = 0; index < node_count; index++) {
			struct node *node = (graph->nodes + index);
			node->length = sequence_lens[index];
			if (node->length > 0) {
				node->sequences = malloc(sizeof(unsigned long long));
				node->sequences[0] = sequences[index];
				node->sequences_len = 1;
			} else {
				node->sequences = 0;
				node->sequences_len = 0;
			}
			node->edges = malloc(sizeof(unsigned long) * edges_out_lens[index]);
			node->edges_len = edges_out_lens[index];
			for (int i = 0; i < node->edges_len; i++) {
				node->edges[i] = edges_out[index][i];
			}
			node->edges_in = malloc(sizeof(unsigned long) * edges_in_lens[index]);
			node->edges_in_len = edges_in_lens[index];
			for (int i = 0; i < node->edges_in_len; i++) {
				node->edges_in[i] = edges_in[index][i];
			}
			node->reference = ref_nodes[index];
		}
	}

	/*
	for (unsigned long index = 0; index < node_count; index++) {
		struct node *node = (graph->nodes + index);
		unsigned long node_id = index;
		if (node->edges_in_len == 0) {
			for (int a = 0; a < 1000; a++) {
				if (node->edges_len > 1) printf("%lu\n", node_id);
				if (node->sequences[0] != 0) {
					printf("%lu: ", node_id);
					for (int i = 0; i < 32; i++) {
						char c = (node->sequences[0] >> ((31 - i) * 2)) & 3;
						if (c == 0) printf("A");
						if (c == 1) printf("C");
						if (c == 2) printf("G");
						if (c == 3) printf("T");
					}
					printf("\n");
				}
				//printf("%d\n", node->edges_len);
				node = (graph->nodes + (node->edges[0]));
				node_id = node->edges[0];
			}
		}
	}
	*/
	
	free(sequences);
	free(sequence_lens);
	free(ref_nodes);
	for (unsigned long index = 0; index < node_count; index++) {
		if (edges_in[index]) free(edges_in[index]);
		if (edges_out[index]) free(edges_out[index]);
	}
	free(edges_in);
	free(edges_in_lens);
	free(edges_out);
	free(edges_out_lens);
	if (id_map) free(id_map);
	
	return 0;
}

/*

void init_node_hash(struct node *n, char *sequence, unsigned int sequence_len, unsigned int *edges, long long edges_len) {
	n->reference = 0;
	n->length = sequence_len;
	if (sequence_len != 0) {
		n->sequences_len = 1 + sequence_len / 31;
	} else {
		n->sequences_len = 0;
	}
	n->sequences = malloc(n->sequences_len * sizeof(unsigned long long));
	for (unsigned long i = 0; i < n->sequences_len; i++) {
		unsigned long remaining = sequence_len - i * 31;
		//n->sequences[i] = hash_max_kmer_by_map(sequence + i * 31, (remaining < 31) ? remaining : 31);
	}
	n->edges = malloc(edges_len * sizeof(unsigned int));
	for (long long i = 0; i < edges_len; i++) n->edges[i] = edges[i];
	n->edges_len = edges_len;
}

void init_node_pack(struct node *n, unsigned char *sequence, unsigned int sequence_len, unsigned int *edges, long long edges_len) {
	n->reference = 0;
	n->length = sequence_len;
	if (sequence_len != 0) {
		n->sequences_len = 1 + sequence_len / 31;
	} else {
		n->sequences_len = 0;
	}
	n->sequences = malloc(n->sequences_len * sizeof(unsigned long long));
	for (unsigned long i = 0; i < n->sequences_len; i++) {
		unsigned long remaining = sequence_len - i * 31;
		n->sequences[i] = pack_max_kmer(sequence + i * 31, (remaining < 31) ? remaining : 31);
	}
	n->edges = malloc(edges_len * sizeof(unsigned int));
	for (long long i = 0; i < edges_len; i++) n->edges[i] = edges[i];
	n->edges_len = edges_len;
}

void from_obgraph(struct graph *graph, unsigned int node_count,
		  unsigned char *sequences, unsigned int *sequence_lens,
		  unsigned int *edges, long long *edges_lens) {
	graph->nodes = malloc(node_count * sizeof(struct node));
	graph->nodes_len = node_count;
	unsigned char *subsequence = sequences;
	unsigned int *subedges = edges;
	unsigned int node_id;
	for (node_id = 0; node_id < node_count; node_id++) {
		init_node_pack(graph->nodes + node_id,
				subsequence,
				sequence_lens[node_id],
				subedges,
				edges_lens[node_id]);
		subsequence += sequence_lens[node_id];
		subedges += edges_lens[node_id];
	}
}

void from_obgraph_npz(char *filepath) {
	char *filename_start = filepath;
	char *a;
	while (1) {
		if ((a = strstr(filename_start, "/"))) {
			filename_start = a + 1;
		} else if ((a = strstr(filename_start, "\\"))) {
			filename_start = a + 1;
		} else {
			break;
		}
	}
	int filename_len = strlen(filename_start) - 4;
	char *filename = malloc(sizeof(char) * (filename_len + 1));
	memcpy(filename, filename_start, sizeof(char) * filename_len);
	filename[filename_len] = '\0';
	FILE *f = fopen(filepath, "rb");
	char buffer_arr[1024];
	char name_buffer[1024];
	short buffer;
	char c;
	const char magic_string[] = "\x93NUMPY";
	int magic_index = 0;
	long last_newline = 0;

	int counter = 0;
	while ((buffer = getc(f)) != EOF) {
		c = buffer;
		if (c == magic_string[magic_index]) {
			magic_index++;
			// if (magic_index > 2) printf("%d: %c %c\n", magic_index, c, magic_string[magic_index - 1]);
			if (magic_index == strlen(magic_string)) {
				char major = getc(f);
				char minor = getc(f);
				if (major != 1) {
					printf("Expected numpy arrays to be of version 1.0. Trying anyway, but may segfault.\n");
				}
				short header_len;
				fread(&header_len, sizeof(short), 1, f);
				printf("Found array! %d.%d, %d\n", major, minor, header_len);
				get_file_str_between(f, name_buffer, last_newline, filename, ".");
				fread(buffer_arr, sizeof(char), 1024, f);
				printf("%s: %s", name_buffer, buffer_arr);
				magic_index = 0;
			}
		} else {
			magic_index = 0;
		}
		if (c == '\n') last_newline = ftell(f);
		counter++;
	}
	printf("%d\n", counter);

	free(filename);

	fclose(f);
}

void get_file_str_between(FILE *f, char *buf, long start_pos, char *from, char *to) {
	long cur_pos = ftell(f);
	char buffer[1024];

	fseek(f, start_pos, SEEK_SET);
	int bytes = fread(buffer, sizeof(char), 1024, f);
	fseek(f, cur_pos, SEEK_SET);

	char *from_pos = strstr(buffer, from);

	printf("test\n");

	if (from_pos == NULL) return;
	from_pos += strlen(from);

	printf("test2\n");
	char *to_pos = strstr(from_pos, to);

	if (to_pos == NULL) return;

	printf("%s\n", from_pos);
	
	//printf("%ld - %ld = %ld\n", to_pos, from_pos, (sizeof(char)) * (to_pos - from_pos));
	memcpy(buf, from_pos, sizeof(char) * (to_pos - from_pos));
	buf[to_pos - from_pos] = '\0';
}
*/

void free_graph(struct graph *graph) {
	for (unsigned long i = 0; i < graph->nodes_len; i++) {
		struct node *n = (graph->nodes + i);
		free(n->sequences);
		free(n->edges);
		free(n->edges_in);
	}
	free(graph->nodes);
	graph->nodes = NULL;
}

int from_file(char *filepath, struct graph *graph) {
	FILE *f = fopen(filepath, "rb");
	int format_code_len = strlen(BCG_FORMAT_CODE);

	for (int i = 0; i < format_code_len; i++) {
		if (fgetc(f) != BCG_FORMAT_CODE[i]) {
			return 1;
		}
	}

	unsigned char version_number = fgetc(f);
	printf("File version: %d\n", version_number);

	for (int i = 0; i < 4; i++) {
		graph->encoding[i] = fgetc(f);
	}
	set_encoding(graph, graph->encoding);

	fread(&(graph->nodes_len), sizeof(unsigned long), 1, f);
	graph->nodes = malloc(sizeof(struct node) * graph->nodes_len);
	memset(graph->nodes, 0, sizeof(struct node) * graph->nodes_len);

	for (unsigned long i = 0; i < graph->nodes_len; i++) {
		struct node *n = (graph->nodes + i);
		fread(&(n->length), sizeof(unsigned long), 1, f);
		fread(&(n->sequences_len), sizeof(unsigned long), 1, f);
		n->sequences = malloc(sizeof(unsigned long long) * n->sequences_len);
		memset(n->sequences, 0, sizeof(unsigned long long) * n->sequences_len);
		fread(n->sequences, sizeof(unsigned long long), n->sequences_len, f);
		unsigned char edges_len = 0;
		fread(&edges_len, sizeof(unsigned char), 1, f);
		n->edges_len = edges_len;
		n->edges = malloc(sizeof(unsigned long) * edges_len);
		memset(n->edges, 0, sizeof(unsigned long) * edges_len);
		fread(n->edges, sizeof(unsigned char), edges_len, f);
		unsigned char reference = 0;
		fread(&reference, sizeof(unsigned char), 1, f);
		n->reference = reference;
	}

	fclose(f);

	return 0;
}

void to_file(char *filepath, struct graph *graph) {
	FILE *f = fopen(filepath, "wb");

	// Indicator for file format
	fwrite(BCG_FORMAT_CODE, sizeof(char), strlen(BCG_FORMAT_CODE), f);

	// Format version number
	fputc(1, f);

	// Encoding
	fwrite(graph->encoding, sizeof(char), 4, f);

	// Write graph info (currently only number of nodes)
	fwrite(&(graph->nodes_len), sizeof(unsigned long), 1, f);

	// Write all nodes
	for (unsigned long i = 0; i < graph->nodes_len; i++) {
		struct node *n = (graph->nodes + i);
		fwrite(&(n->length), sizeof(unsigned long), 1, f);
		fwrite(&(n->sequences_len), sizeof(unsigned long), 1, f);
		fwrite(n->sequences, sizeof(unsigned long long), n->sequences_len, f);
		unsigned char edges_len = n->edges_len;
		fwrite(&edges_len, sizeof(unsigned char), 1, f);
		fwrite(n->edges, sizeof(unsigned char), n->edges_len, f);
		unsigned char reference = n->reference;
		fwrite(&reference, sizeof(unsigned char), 1, f);
	}

	fclose(f);
}

void set_encoding(struct graph *graph, const char *encoding) {
	memcpy(graph->encoding, encoding, sizeof(char) * 4);
	fill_map_by_encoding(graph->encoding_map, encoding);
}

/*
int main(int argc, char** argv) {
	if (argc != 2) {
		printf("This program requires one argument: A filename for an npz file.\n");
	}
	struct graph graph;
	from_file(argv[1], &graph);
	printf("%ld nodes\n", graph.nodes_len);
	return 0;
}
*/
