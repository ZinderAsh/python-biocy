#ifndef GRAPH_H
#define GRAPH_H

struct node {
	unsigned int len;
	unsigned long long *sequences;
	unsigned short sequences_len;
	unsigned long *edges;
	unsigned char edges_len:7;
	unsigned char reference:1;
};

struct graph {
	struct node *nodes;
	unsigned long nodes_len;
};

#endif
