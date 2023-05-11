#ifndef KIVS_NODE_H
#define KIVS_NODE_H

#include <stdint.h>

struct node {
	uint32_t length;
	uint64_t *sequences;
	uint32_t sequences_len;
	uint32_t *edges;
	uint32_t *edges_in;
	uint8_t edges_len;
	uint8_t edges_in_len;
	uint32_t reference_index;
	bool reference;
};

#endif
