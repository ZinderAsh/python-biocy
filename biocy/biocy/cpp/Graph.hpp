#ifndef BIOCY_GRAPH_H
#define BIOCY_GRAPH_H

#include "node.hpp"
#include <cstring>
#include <cstdlib>
#include <stdint.h>
#include <tuple>
#include <stdio.h>

#define BCG_FORMAT_CODE "BIOCYGRAPH"

class Graph {
public:
	struct node *nodes;
	uint32_t nodes_len;
	char encoding[4];
	uint8_t encoding_map[256];	

	~Graph() {
		if (nodes != nullptr) {
			for (uint32_t i = 0; i < nodes_len; i++) {
				free((nodes + i)->sequences);
				free((nodes + i)->edges);
				free((nodes + i)->edges_in);
			}
		}
		free(nodes);
	}
	
	static Graph *FromFile(char *filepath);
	static Graph *FromGFAFile(char *filepath);
	static Graph *FromGFAFileEncoded(char *filepath, const char *encoding);

	void Compress();
	void AddEmptyNodes();

	void ToFile(char *filepath);

private:
	Graph(const char *encoding) {
		nodes = nullptr;
		nodes_len = 0;
		this->SetEncoding(encoding);
	}

	void SetEncoding(const char *encoding);
	static std::tuple<uint32_t, uint32_t> GFAGetNodeIDRange(FILE *f);

	uint32_t *CreateCompressedIDMap(uint32_t *return_compressed_node_count);
};

#endif
