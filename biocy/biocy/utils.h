#ifndef BIOCY_UTILS_H
#define BIOCY_UTILS_H

void fill_map_by_encoding(unsigned char *map, const char *encoding);
unsigned long long hash_min_kmer_by_map(char *str, unsigned char k, unsigned char *map);
unsigned long long hash_max_kmer_by_map(char *str, unsigned char k, unsigned char *map);
unsigned long long hash_kmer_by_map(char *str, unsigned char k, unsigned char *map);
unsigned long long pack_min_kmer(unsigned char *arr, unsigned char k);
unsigned long long pack_max_kmer(unsigned char *arr, unsigned char k);
unsigned long long pack_kmer(char *arr, unsigned char k);

#endif
