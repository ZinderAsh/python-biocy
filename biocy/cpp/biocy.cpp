#include <stdlib.h>
#include <stdio.h>
#include <stack>
#include <vector>
#include <algorithm>
#include <cmath>

#include "KmerFinder.hpp"
#include "node.hpp"
#include "VCF.hpp"
#include "FASTA.hpp"

int main(int argc, char** argv) {
	/*
	if (argc != 4) return 1;
	Graph *graph = Graph::FromFastaVCF(argv[1], argv[2], 21);
	//graph->ToFile(argv[3]);

	uint8_t k = 31;
	KmerFinder *kf = new KmerFinder(graph, k, 200);
	kf->Find();
	printf("Found %lu kmers\n", kf->found_count);
	kf->FindVariantSignatures();
	delete kf;

	delete graph;

	graph = Graph::FromFile(argv[3]);

	kf = new KmerFinder(graph, k, 200);
	kf->Find();
	printf("Found %lu kmers\n", kf->found_count);
	delete kf;

	delete graph;
	*/
	
	Graph *graph = new Graph("ACGT");

	uint32_t node_0 = graph->AddNode("ACTGACTGACTG");
	uint32_t node_1 = graph->AddNode("G");
	uint32_t node_2 = graph->AddNode("");
	uint32_t node_3 = graph->AddNode("AT");
	uint32_t node_4 = graph->AddNode("AACTG");
	uint32_t node_5 = graph->AddNode("CTA");
	uint32_t node_6 = graph->AddNode("CTGCTTTTTTGTATA");

	graph->Get(node_0)->reference = true;
	graph->Get(node_1)->reference = true;
	graph->Get(node_3)->reference = true;
	graph->Get(node_4)->reference = true;
	graph->Get(node_6)->reference = true;

	graph->AddEdge(node_0, node_1);
	graph->AddEdge(node_0, node_2);
	graph->AddEdge(node_1, node_3);
	graph->AddEdge(node_2, node_3);
	graph->AddEdge(node_3, node_4);
	graph->AddEdge(node_3, node_5);
	graph->AddEdge(node_4, node_6);
	graph->AddEdge(node_5, node_6);

	graph->AddInEdges();
	
	uint8_t k = 5;
	KmerFinder *kf = new KmerFinder(graph, k, 31);
	kf->SetFlag(FLAG_TO_STDOUT, false);
	kf->FindKmersSpanningNode(node_4);
	//auto windows = kf->FindWindowsForVariant(node_4, node_5);
	//kf->SetFlag(FLAG_ALIGN_SIGNATURE_WINDOWS, false);
	//VariantWindow *min_window = kf->FindVariantSignatures(node_4, node_5);
	//kf->FindKmersForVariant(node_4, node_5);

	//delete min_window;
	delete kf;
	delete graph;
	
	/*
	printf("window kmers\n");
	for (uint64_t i = 0; i < kf->found_count; i++) {
		uint64_t kmer_hash = kf->found_kmers[i];
		char *kmer = graph->DecodeKmer(kmer_hash, k);
		printf("%s: %2u %2u %2u %2u\n",
				kmer,
				kf->found_nodes[i],
				kf->found_node_sequence_start_positions[i],
				kf->found_node_sequence_kmer_positions[i],
				index[kmer_hash]);
		free(kmer);
	}
	*/
	/*
	printf("variant windows\n");
	for (uint32_t i = 0; i < windows.size(); i++) {
		windows[i]->Print(kf);
	}

	printf("rarest window\n");
	min_window->Print(kf);

	for (uint32_t i = 0; i < windows.size(); i++) {
		delete windows[i];
	}
	delete min_window;

	delete kf;
	delete graph;

	return 0;
	*/
	/*
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
	*/
}
