#include <stdlib.h>
#include <stdio.h>
#include <stack>
#include <vector>
#include <algorithm>
#include <cmath>

#include "KmerFinder.hpp"
#include "node.hpp"

int main(int argc, char** argv) {
	if (argc != 2) {
		printf("This program requires one argument: A filename for an npz file.\n");
		return 1;
	}
	Graph *graph = Graph::FromGFAFile(argv[1]);
	printf("%d nodes\n", graph->nodes_len);	

	graph->Compress();
	printf("%d nodes\n", graph->nodes_len);

	KmerFinder *kf = new KmerFinder(graph, 31, 31);
	kf->Find();
	printf("Done. %ld kmers found\n", kf->found_count);
	
	auto kmer_index = kf->CreateKmerIndex();
	printf("Kmer Count Size: %ld\n", kmer_index.size());
	printf("Kmer Count Bucket Count: %ld\n", kmer_index.bucket_count());
	
	delete kf;

	//uint32_t empty_node_count = graph->AddEmptyNodes();
	//printf("Added %u empty nodes.\n", empty_node_count);

	//empty_node_count = graph->AddEmptyNodes();
	//printf("Added another %u empty nodes.\n", empty_node_count);
	
	delete graph;
	return 0;
	graph->ToFile(argv[2]);
	
	kf = new KmerFinder(graph, 31, 31);
	kf->Find();
	printf("Done. %ld kmers found\n", kf->found_count);
	delete kf;
	
	delete graph;
	
	return 0;
}
