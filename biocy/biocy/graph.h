#ifndef GRAPH_H
#define GRAPH_H

#define BCG_FORMAT_CODE "BIOCYGRAPH"
#define GRAPH_OPTIMIZE 0x01
#define GRAPH_ADD_EMTPY 0x02

struct node {
	unsigned long length;
	unsigned long long *sequences;
	unsigned long sequences_len;
	unsigned long *edges;
	unsigned long *edges_in;
	unsigned char edges_len;
	unsigned char edges_in_len;
	unsigned char reference:1;
};

struct graph {
	struct node *nodes;
	unsigned long nodes_len;
	char encoding[4];
	unsigned char encoding_map[256];
};

void from_obgraph(struct graph *graph, unsigned int node_count,
		  unsigned char *sequences, unsigned int *sequence_lens,
		  unsigned int *edges, long long *edges_lens);

void free_graph(struct graph *graph);

void set_encoding(struct graph *graph, const char *encoding);

int from_file(char *filepath, struct graph *graph);
void to_file(char *filepath, struct graph *graph);

int from_file_gfa(char *filepath, struct graph *graph, char flags);
int from_file_gfa_encoded(char *filepath, struct graph *graph, const char *encoding, char flags);
int from_file_gfa_by_map(char *filepath, struct graph *graph, unsigned char *map, char flags);

#endif
