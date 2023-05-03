#ifndef BIOCY_GRAPH_H
#define BIOCY_GRAPH_H

#include "node.hpp"

#include <cstring>
#include <cstdlib>
#include <stdint.h>
#include <tuple>
#include <stdio.h>

#include "hashing.hpp"

#define BCG_FORMAT_CODE "BIOCYGRAPH"

class Graph {
public:
	struct node *nodes;
	uint32_t nodes_len;
	char encoding[4];
	uint8_t encoding_map[256];	
	
	Graph(const char *encoding) {
		nodes = NULL;
		nodes_len = 0;
		this->SetEncoding(encoding);
	}

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
	static Graph *FromFastaVCF(char *fasta_filepath, char *vcf_filepath, int16_t chromosome);
	static Graph *FromFastaVCFEncoded(char *fasta_filepath, char *vcf_filepath, int16_t chromosome, const char *encoding);

	struct node *Get(uint32_t node_id) {
		return nodes + node_id;
	}

	void Compress();
	uint32_t AddEmptyNodes();

	void AddInEdges() {
		struct node *nodes_end = nodes + nodes_len;
		struct node *n;
		uint8_t edge_index;
		for (n = nodes; n < nodes_end; n++) {
			if (n->edges_in) free(n->edges_in);
			n->edges_in_len = 0;
		}
		for (n = nodes; n < nodes_end; n++) {
			for (edge_index = 0; edge_index < n->edges_len; edge_index++) {
				(nodes + n->edges[edge_index])->edges_in_len++;
			}
		}
		for (n = nodes; n < nodes_end; n++) {
			n->edges_in = (uint32_t *) malloc(sizeof(uint32_t) * n->edges_in_len);
			n->edges_in_len = 0;
		}
		for (uint32_t i = 0; i < nodes_len; i++) {
			n = nodes + i;
			for (edge_index = 0; edge_index < n->edges_len; edge_index++) {
				struct node *in_node = nodes + n->edges[edge_index];
				in_node->edges_in[in_node->edges_in_len++] = i;
			}
		}
	}

	void ToFile(char *filepath);

	uint64_t HashMinKmer(const char *str, uint8_t k) {
		return hash_min_kmer_by_map(str, k, encoding_map);
	}
	uint64_t HashMaxKmer(const char *str, uint8_t k) {
		return hash_max_kmer_by_map(str, k, encoding_map);
	}
	uint64_t HashKmer(const char *str, uint8_t k) { return HashMinKmer(str, k); }

	char *DecodeKmer(uint64_t hash, uint8_t k) {
		return decode_kmer_by_map(hash, k, encoding_map);
	}

	uint32_t GetRootNodeID();
	uint32_t GetLastNodeID();
	uint32_t GetReferenceNodeID(uint32_t reference_index);
	uint32_t GetNextReferenceNodeID(uint32_t previous_id);

	uint32_t AddNode(const char *sequence);
	void AddEdge(uint32_t from_node_id, uint32_t to_node_id);
	uint32_t AppendEmptyNode();

private:
	void SetEncoding(const char *encoding);
	static std::tuple<uint32_t, uint32_t> GFAGetNodeIDRange(FILE *f);

	uint32_t *CreateCompressedIDMap(uint32_t *return_compressed_node_count);
	uint32_t GetRequiredEmptyNodeCount();
	uint32_t GetRequiredEmptyNodesFromNode(uint32_t from_node_id);
	bool NodeHasEdge(uint32_t node_id, uint32_t edge_id);
	void InitializeEmptyNode(uint32_t node_id);
	uint32_t CreateEmptyNodes();
	void MoveEdgesToIntermediateNode(uint32_t from_node_id, uint32_t to_node_id, uint32_t mid_node_id);
};

#endif
