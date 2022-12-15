#ifndef BIOCY_HASHING_H
#define BIOCY_HASHING_H

#include <stdint.h>

void fill_map_by_encoding(uint8_t *map, const char *encoding);
uint64_t hash_min_kmer_by_map(char *str, uint8_t k, uint8_t *map);
uint64_t hash_max_kmer_by_map(char *str, uint8_t k, uint8_t *map);
uint64_t hash_kmer_by_map(char *str, uint8_t k, uint8_t *map);
uint64_t pack_min_kmer(char *arr, uint8_t k);
uint64_t pack_max_kmer(char *arr, uint8_t k);
uint64_t pack_kmer(char *arr, uint8_t k);
uint64_t pack_max_kmer_with_offset(char *arr, uint32_t offset, uint8_t k);

#endif
