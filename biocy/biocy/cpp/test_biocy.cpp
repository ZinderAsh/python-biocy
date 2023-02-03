#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <stdint.h> 
#include <string.h>
#include <stdio.h>
#include <unordered_map>

#include "Graph.hpp"
#include "hashing.hpp"
#include "KmerFinder.hpp"
#include "node.hpp"

void fill_index(Graph *graph, std::unordered_map<uint64_t, uint32_t> *index, const char **kmers, const uint32_t *counts, uint32_t len) {
	for (uint32_t i = 0; i < len; i++) {
		(*index)[graph->HashMinKmer(kmers[i], strlen(kmers[i]))] = counts[i];
	}
}

uint32_t count_kmer(KmerFinder *kf, const char *kmer) {
	uint64_t kmer_hash = kf->graph->HashMinKmer(kmer, strlen(kmer));
	uint32_t count = 0;
	for (uint64_t i = 0; i < kf->found_count; i++) {
		if (kf->found_kmers[i] == kmer_hash) {
			count++;
		}
	}
	return count;
}

TEST_CASE("Fill map and decode kmer hashes with the map.") {
	uint8_t map[256];

	uint64_t kmers[] = {
		0b00011011,
		0b1100,
		0b1111111100000000,
		0b01010101101010100101010110101010,
		0b0101010110101010010101011010101011111111000000001111111100000000
	};
	
	SUBCASE("ACGT") {
		fill_map_by_encoding(map, "ACGT");
		CHECK(map['A'] == 0);
		CHECK(map['a'] == 0);
		CHECK(map['C'] == 1);
		CHECK(map['c'] == 1);
		CHECK(map['G'] == 2);
		CHECK(map['g'] == 2);
		CHECK(map['T'] == 3);
		CHECK(map['t'] == 3);

		CHECK(map['N'] == 0);
		CHECK(map['n'] == 0);

		CHECK(map[0] == 'A');
		CHECK(map[1] == 'C');
		CHECK(map[2] == 'G');
		CHECK(map[3] == 'T');

		CHECK(strcmp(decode_kmer_by_map(kmers[0], 4, map), "ACGT") == 0);
		CHECK(strcmp(decode_kmer_by_map(kmers[1], 2, map), "TA") == 0);
		CHECK(strcmp(decode_kmer_by_map(kmers[2], 8, map), "TTTTAAAA") == 0);
		CHECK(strcmp(decode_kmer_by_map(kmers[3], 16, map), "CCCCGGGGCCCCGGGG") == 0);
		CHECK(strcmp(decode_kmer_by_map(kmers[4], 32, map), "CCCCGGGGCCCCGGGGTTTTAAAATTTTAAAA") == 0);
	}

	SUBCASE("cGaT") {
		fill_map_by_encoding(map, "cGaT");
		CHECK(map['A'] == 2);
		CHECK(map['a'] == 2);
		CHECK(map['C'] == 0);
		CHECK(map['c'] == 0);
		CHECK(map['G'] == 1);
		CHECK(map['g'] == 1);
		CHECK(map['T'] == 3);
		CHECK(map['t'] == 3);

		CHECK(map['N'] == 0);
		CHECK(map['n'] == 0);

		CHECK(map[0] == 'C');
		CHECK(map[1] == 'G');
		CHECK(map[2] == 'A');
		CHECK(map[3] == 'T');

		CHECK(strcmp(decode_kmer_by_map(kmers[0], 4, map), "CGAT") == 0);
		CHECK(strcmp(decode_kmer_by_map(kmers[1], 2, map), "TC") == 0);
		CHECK(strcmp(decode_kmer_by_map(kmers[2], 8, map), "TTTTCCCC") == 0);
		CHECK(strcmp(decode_kmer_by_map(kmers[3], 16, map), "GGGGAAAAGGGGAAAA") == 0);
		CHECK(strcmp(decode_kmer_by_map(kmers[4], 32, map), "GGGGAAAAGGGGAAAATTTTCCCCTTTTCCCC") == 0);
	}

	SUBCASE("gATc") {
		fill_map_by_encoding(map, "gATc");
		CHECK(map['A'] == 1);
		CHECK(map['a'] == 1);
		CHECK(map['C'] == 3);
		CHECK(map['c'] == 3);
		CHECK(map['G'] == 0);
		CHECK(map['g'] == 0);
		CHECK(map['T'] == 2);
		CHECK(map['t'] == 2);

		CHECK(map['N'] == 0);
		CHECK(map['n'] == 0);

		CHECK(map[0] == 'G');
		CHECK(map[1] == 'A');
		CHECK(map[2] == 'T');
		CHECK(map[3] == 'C');

		CHECK(strcmp(decode_kmer_by_map(kmers[0], 4, map), "GATC") == 0);
		CHECK(strcmp(decode_kmer_by_map(kmers[1], 2, map), "CG") == 0);
		CHECK(strcmp(decode_kmer_by_map(kmers[2], 8, map), "CCCCGGGG") == 0);
		CHECK(strcmp(decode_kmer_by_map(kmers[3], 16, map), "AAAATTTTAAAATTTT") == 0);
		CHECK(strcmp(decode_kmer_by_map(kmers[4], 32, map), "AAAATTTTAAAATTTTCCCCGGGGCCCCGGGG") == 0);
	}
}

TEST_CASE("Test kmer hashing.") {
	int encodings_len = 4;
	const char *encodings[] = {
		"ACTG",
		"ctAG",
		"gAtC",
		"cGat"
	};
	int kmers_len = 6;
	const char *kmers[] = {
		"ACgt",
		"TA",
		"CCCCGGGG",
		"AaaATTttAAAATTTT",
		"AAAATTttAAAATTTTCCCCGGGGCCCCGGGG",
		"GCTgcATgCGATCGTAACGTtcACTGCAATCT"
	};
	char upper_kmer[33];

	for (int i = 0; i < encodings_len; i++) {
		uint8_t map[256];
		fill_map_by_encoding(map, encodings[i]);
		for (int j = 0; j < kmers_len; j++) {
			CAPTURE(i);
			CAPTURE(j);
			const char *kmer = kmers[j];
			uint8_t kmer_len = strlen(kmer);
			for (int k = 0; k < kmer_len; k++) {
				upper_kmer[k] = toupper(kmer[k]);
			}
			uint64_t hash = hash_min_kmer_by_map(kmer, kmer_len, map);
			char *decoded = decode_kmer_by_map(hash, kmer_len, map);
			CHECK(strncmp(upper_kmer, decoded, kmer_len) == 0);
			free(decoded);
			hash = hash_max_kmer_by_map(kmer, kmer_len, map);
			decoded = decode_kmer_by_map(hash >> (64 - kmer_len * 2), kmer_len, map);
			CHECK(strncmp(upper_kmer, decoded, kmer_len) == 0);
			free(decoded);
		}
	}
}

TEST_CASE("Kmer packing") {
	char kmers[6][38] {
		{0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{1, 1, 1, 1, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{2, 2, 1, 1, 3, 3, 0, 0, 2, 2, 1, 1, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0},
		{2, 2, 1, 1, 3, 3, 0, 0, 2, 2, 1, 1, 3, 3, 0, 0, 1, 1, 1, 1, 3, 3, 3, 3, 1, 1, 1, 1, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0},
		{0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0}
	};
	
	CHECK(pack_min_kmer(kmers[0], 4) == 0b00011011);
	CHECK(pack_max_kmer(kmers[0], 4) == 0b0001101100000000000000000000000000000000000000000000000000000000);
	CHECK(pack_min_kmer(kmers[1], 2) == 0b1001);
	CHECK(pack_max_kmer(kmers[1], 2) == 0b1001000000000000000000000000000000000000000000000000000000000000);
	CHECK(pack_min_kmer(kmers[2], 8) == 0b0101010111111111);
	CHECK(pack_max_kmer(kmers[2], 8) == 0b0101010111111111000000000000000000000000000000000000000000000000);
	CHECK(pack_min_kmer(kmers[3], 16) == 0b10100101111100001010010111110000);
	CHECK(pack_max_kmer(kmers[3], 16) == 0b1010010111110000101001011111000000000000000000000000000000000000);
	CHECK(pack_min_kmer(kmers[4], 32) == 0b1010010111110000101001011111000001010101111111110101010111111111);
	CHECK(pack_max_kmer(kmers[4], 32) == 0b1010010111110000101001011111000001010101111111110101010111111111);
	CHECK(pack_max_kmer_with_offset(kmers[0], 1, 3) == 0b0110110000000000000000000000000000000000000000000000000000000000);
	CHECK(pack_max_kmer_with_offset(kmers[5], 16, 22) == 0b0000000000000000000000000000000001010101000000000000000000000000);
}

TEST_CASE("Test small graph.") {
	
	Graph *graph = new Graph("ACGT");

	uint32_t node_0 = graph->AddNode("ACTGACTGACTG");
	uint32_t node_1 = graph->AddNode("G");
	uint32_t node_2 = graph->AddNode("T");
	uint32_t node_3 = graph->AddNode("AT");
	uint32_t node_4 = graph->AddNode("ACT");
	uint32_t node_5 = graph->AddNode("CTA");
	uint32_t node_6 = graph->AddNode("CTGCTTTTTTTT");

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

	REQUIRE(graph->nodes_len == 7);

	SUBCASE("Test kmer finder for k = 4") {
		KmerFinder *kf = new KmerFinder(graph, 4, 31);
		kf->Find();
		
		CHECK(count_kmer(kf, "ACTG") == 5);
		CHECK(count_kmer(kf, "CTGG") == 2);
		CHECK(count_kmer(kf, "ATAC") == 2);
		CHECK(count_kmer(kf, "AAAA") == 0);
		CHECK(count_kmer(kf, "CTGC") == 1);

		delete kf;
	}

	SUBCASE("Testing variant windows for k = 5") {
		KmerFinder *kf = new KmerFinder(graph, 5, 31);

		SUBCASE("Testing for variant on node 1 and 2") {
			kf->FindKmersForVariant(node_1, node_2);

			CHECK(count_kmer(kf, "TGTAT") == 1);
			CHECK(count_kmer(kf, "GATAC") == 1);
			CHECK(count_kmer(kf, "GATCT") == 1);
			CHECK(count_kmer(kf, "ACTGT") == 1);
			CHECK(count_kmer(kf, "TGGAT") == 1);
		}

		SUBCASE("Testing for variant on node 4 and 5") {
			kf->FindKmersForVariant(node_4, node_5);

			CHECK(count_kmer(kf, "ACTCT") == 1);
			CHECK(count_kmer(kf, "CTACT") == 1);
			CHECK(count_kmer(kf, "GATAC") == 1);
			CHECK(count_kmer(kf, "TATCT") == 1);
			CHECK(count_kmer(kf, "TCTAC") == 1);
			CHECK(count_kmer(kf, "GTATA") == 1);
		}

		delete kf;
	}

	SUBCASE("Test kmer index for k = 4") {
		KmerFinder *kf = new KmerFinder(graph, 4, 31);
		kf->Find();
		auto kmer_index = kf->CreateKmerIndex();

		CHECK(kmer_index[graph->HashKmer("ACTG", 4)] == 5);
		CHECK(kmer_index[graph->HashKmer("CTGG", 4)] == 2);
		CHECK(kmer_index[graph->HashKmer("ATAC", 4)] == 2);
		CHECK(kmer_index[graph->HashKmer("AAAA", 4)] == 0);
		CHECK(kmer_index[graph->HashKmer("CTGC", 4)] == 1);
	}

	delete graph;
}

TEST_CASE("Test finding minimal variant windows.") {

	SUBCASE("Variant and reference of equal length.") {
		Graph *graph = new Graph("ACGT");

		uint32_t node_0 = graph->AddNode("ACTGACTGACTG");
		uint32_t node_1 = graph->AddNode("G");
		uint32_t node_2 = graph->AddNode("T");
		uint32_t node_3 = graph->AddNode("AT");
		uint32_t node_4 = graph->AddNode("ACT");
		uint32_t node_5 = graph->AddNode("CTA");
		uint32_t node_6 = graph->AddNode("CTGCTTTTTTTT");

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

		REQUIRE(graph->nodes_len == 7);

		std::unordered_map<uint64_t, uint32_t> index;

		const char *kmers[] = {
			"GGATA", "GTATA", "GGATC", "GTATC",
			"GATAC", "TATAC", "GATCT", "TATCT",
			"ATACT", "ATCTA",
			"TACTC", "TCTAC",
			"ACTCT", "CTACT",
			"CTCTG", "TACTG",
			"TCTGC", "ACTGC"
		};
		
		SUBCASE("Test 1") {
			const uint32_t counts[] = {
				2, 4, 2, 3,
				7, 7, 7, 7,
				7, 7,
				7, 7,
				7, 7,
				7, 7,
				7, 7
			};
			fill_index(graph, &index, kmers, counts, 18);
			KmerFinder *kf = new KmerFinder(graph, 5, 31);
			kf->SetKmerFrequencyIndex(index);

			VariantWindow *min_window = kf->FindRarestWindowForVariant(node_4, node_5);

			CHECK(min_window->max_frequency == 7);
			CHECK(min_window->reference_kmers_len == 2);
			CHECK(min_window->variant_kmers_len == 2);
			CHECK(min_window->reference_kmers[0] == graph->HashKmer("GGATA", 5));
			CHECK(min_window->reference_kmers[1] == graph->HashKmer("GTATA", 5));
			CHECK(min_window->variant_kmers[0] == graph->HashKmer("GGATC", 5));
			CHECK(min_window->variant_kmers[1] == graph->HashKmer("GTATC", 5));

			delete min_window;

			delete kf;
		}
		
		SUBCASE("Test 2") {
			const uint32_t counts[] = {
				2, 5, 2, 3,
				7, 7, 7, 7,
				7, 7,
				7, 7,
				4, 3,
				7, 7,
				7, 7
			};
			fill_index(graph, &index, kmers, counts, 18);
			KmerFinder *kf = new KmerFinder(graph, 5, 31);
			kf->SetKmerFrequencyIndex(index);

			VariantWindow *min_window = kf->FindRarestWindowForVariant(node_4, node_5);

			CHECK(min_window->max_frequency == 7);
			CHECK(min_window->reference_kmers_len == 1);
			CHECK(min_window->variant_kmers_len == 1);
			CHECK(min_window->reference_kmers[0] == graph->HashKmer("ACTCT", 5));
			CHECK(min_window->variant_kmers[0] == graph->HashKmer("CTACT", 5));

			delete min_window;

			delete kf;
		}

		delete graph;
	}

	SUBCASE("Variant and reference of different length.") {
		Graph *graph = new Graph("ACGT");

		uint32_t node_0 = graph->AddNode("ACTGACTGACTG");
		uint32_t node_1 = graph->AddNode("G");
		uint32_t node_2 = graph->AddNode("");
		uint32_t node_3 = graph->AddNode("AT");
		uint32_t node_4 = graph->AddNode("AACTG");
		uint32_t node_5 = graph->AddNode("CTA");
		uint32_t node_6 = graph->AddNode("CTGCTTTTTTTT");

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

		REQUIRE(graph->nodes_len == 7);

		std::unordered_map<uint64_t, uint32_t> index;

		const char *kmers[] = {
			"GGATA", "TGATA", "GGATC", "TGATC",
			"GATAA", "GATAA", "GATCT",
			"ATAAC", "ATCTA",
			"TAACT", "TCTAC",
			"AACTG", "CTACT",
			"ACTGC", "TACTG",
			"CTGCT", "ACTGC",
			"GCTGC", /*"ACTGC",*/
			"TGCTG", /*"TACTG",*/
			"CTGCT", /*"CTACT",*/
			"ACTGC", /*"TCTAC",*/
			"AACTG", /*"ATCTA",*/
			"TAACT", /*"GATCT",*/
			"ATAAC", /*"GGATC", "TGATC"*/
		};
		
		SUBCASE("Test 1") {
			const uint32_t counts[] = {
				2, 4, 2, 3,
				7, 7, 7,
				7, 7,
				7, 7,
				7, 7,
				7, 7,
				7, 7,
				7,
				7,
				7,
				7,
				7,
				7,
				7
			};
			fill_index(graph, &index, kmers, counts, 24);
			KmerFinder *kf = new KmerFinder(graph, 5, 31);
			kf->SetKmerFrequencyIndex(index);

			VariantWindow *min_window = kf->FindRarestWindowForVariant(node_4, node_5);

			CHECK(min_window->max_frequency == 7);
			CHECK(min_window->reference_kmers_len == 2);
			CHECK(min_window->variant_kmers_len == 2);
			CHECK(min_window->reference_kmers[0] == graph->HashKmer("GGATA", 5));
			CHECK(min_window->reference_kmers[1] == graph->HashKmer("TGATA", 5));
			CHECK(min_window->variant_kmers[0] == graph->HashKmer("GGATC", 5));
			CHECK(min_window->variant_kmers[1] == graph->HashKmer("TGATC", 5));

			delete min_window;

			delete kf;
		}
		
	}
}



