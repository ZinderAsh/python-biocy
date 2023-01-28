#include "hashing.hpp"
#include <stdio.h>
#include <stdlib.h>

void fill_map_by_encoding(uint8_t *map, const char *encoding) {
	for (uint8_t i = 0; i < 4; i++) {
		map[encoding[i] | 0x20] = i;
		map[encoding[i] & 0xDF] = i;
		map[i] = encoding[i] & 0xDF;
	}
	map['N'] = 0;
	map['n'] = 0;
}

// Hashes a string of bases (max 31) to a 2-bit encoded long long
uint64_t hash_min_kmer_by_map(const char *str, uint8_t k, uint8_t *map) {
	uint64_t hashed = 0;
	for (uint8_t i = 0; i < k; i++) {
		uint64_t val = map[str[i]];
		hashed |= val << ((k - i - 1) * 2);
	}
	return hashed;
}

// Same as hash_kmer, except the bits are left-aligned in the long long.
// This is faster, but also practical for the kmer_finder algorithm
uint64_t hash_max_kmer_by_map(const char *str, uint8_t k, uint8_t *map) {
	uint64_t hashed = 0;
	for (uint8_t i = 0; i < k; i++) {
		uint64_t val = map[str[i]];
		hashed |= val << (62 - i * 2);
	}
	return hashed;
}

// Alias for intuitive usage
uint64_t hash_kmer_by_map(char *str, uint8_t k, uint8_t *map) {
	return hash_min_kmer_by_map(str, k, map);
}

uint64_t hash_min_kmer_by_encoding(const char *str, uint8_t k, const char *encoding) {
	uint8_t map[256];
	fill_map_by_encoding(map, encoding);
	return hash_min_kmer_by_map(str, k, map);
}

char *decode_kmer_by_map(uint64_t hash, uint8_t k, uint8_t *map) {
	char *kmer = (char *) malloc(sizeof(char) * (k + 1));
	for (uint8_t i = 0; i < k; i++) {
		uint8_t base = (hash >> ((k - i - 1) * 2)) & 3;
		kmer[i] = map[base];
	}
	kmer[k] = 0;
	return kmer;
}

// Packs an array of 2-bit encoded values into a single long long, right-aligned
uint64_t pack_min_kmer(char *arr, uint8_t k) {
	uint64_t packed = 0;
	for (uint8_t i = 0; i < k; i++) {
		uint64_t val = arr[i];
		packed |= (val << ((k - i - 1) * 2));
	}
	return packed;
}

// Packs an array of 2-bit encoded values into a single long long, left-aligned
// This is more practical in the kmer_finder algorithm
uint64_t pack_max_kmer(char *arr, uint8_t k) {
	uint64_t packed = 0;
	for (uint8_t i = 0; i < k; i++) {
		uint64_t val = arr[i];
		packed |= (val << (62 - i * 2));
	}
	return packed;
}

// Alias for intuitive usage
uint64_t pack_kmer(char *arr, uint8_t k) {
	return pack_min_kmer(arr, k);
}

uint64_t pack_max_kmer_with_offset(char *arr, uint32_t offset, uint8_t k) {
	return pack_max_kmer(arr + offset, k);
}

