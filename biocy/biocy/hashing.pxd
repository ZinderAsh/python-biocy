from libc.stdint cimport uint8_t, uint32_t, uint64_t

cdef extern from "cpp/hashing.hpp":
    void fill_map_by_encoding(uint8_t *map, char *encoding)
    char *decode_kmer_by_map(uint64_t hash, uint8_t k, char *map)
    uint64_t hash_min_kmer_by_map(char *str, uint8_t k, uint8_t *map)
    uint64_t hash_max_kmer_by_map(char *str, uint8_t k, uint8_t *map)
    uint64_t hash_kmer_by_map(char *str, uint8_t k, uint8_t *map)
    uint64_t hash_min_kmer_by_encoding(char *str, uint8_t k, char *encoding)
    uint64_t pack_min_kmer(char *arr, uint8_t k)
    uint64_t pack_max_kmer(char *arr, uint8_t k)
    uint64_t pack_kmer(char *arr, uint8_t k)
    uint64_t pack_max_kmer_with_offset(char *arr, uint32_t offset, uint8_t k)
