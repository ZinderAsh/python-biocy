#include "utils.h"
#include <stdio.h>

void fill_map_by_encoding(unsigned char *map, const char *encoding) {
	for (int i = 0; i < 4; i++) {
		map[encoding[i] | 0x20] = i;
		map[encoding[i] & 0xDF] = i;
	}
	map['N'] = 0;
	map['n'] = 0;
}

// Hashes a string of bases (max 31) to a 2-bit encoded long long
unsigned long long hash_min_kmer_by_map(char *str, unsigned char k, unsigned char *map) {
	unsigned long long hashed = 0;
	for (unsigned char i = 0; i < k; i++) {
		hashed |= map[str[i]] << ((k - i) * 2);
	}
	return hashed >> 1;
}

// Same as hash_kmer, except the bits are left-aligned in the long long.
// This is faster, but also practical for the kmer_finder algorithm
unsigned long long hash_max_kmer_by_map(char *str, unsigned char k, unsigned char *map) {
	unsigned long long hashed = 0;
	for (unsigned char i = 0; i < k; i++) {
		char c = str[i];
		if (c != 'A' && c != 'a' && c != 'C' && c != 'c' &&
				c != 'G' && c != 'g' && c != 'T' && c != 't' &&
				c != 'N' && c != 'n') printf("%c\n", c);
		hashed |= (map[str[i]] * 1L) << (62 - i * 2);
	}
	return hashed;
}

// Alias for intuitive usage
unsigned long long hash_kmer_by_map(char *str, unsigned char k, unsigned char *map) {
	return hash_min_kmer_by_map(str, k, map);
}

// Packs an array of 2-bit encoded values into a single long long, right-aligned
unsigned long long pack_min_kmer(unsigned char *arr, unsigned char k) {
	unsigned long long packed = 0;
	for (unsigned char i = 0; i < k; i++) {
		packed |= (arr[i] << (i * 2));
	}
	return packed;
}

// Packs an array of 2-bit encoded values into a single long long, left-aligned
// This is more practical in the kmer_finder algorithm
unsigned long long pack_max_kmer(unsigned char *arr, unsigned char k) {
	unsigned long long packed = 0;
	for (unsigned char i = 0; i < k; i++) {
		packed |= (arr[i] << (62 - i * 2));
	}
	return packed;
}

// Alias for intuitive usage
unsigned long long pack_kmer(char *arr, unsigned char k) {
	return pack_min_kmer(arr, k);
}

