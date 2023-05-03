#ifndef BIOCY_FILE_READER_VCF
#define BIOCY_FILE_READER_VCF

#include <cstdlib>
#include <stdio.h>
#include <stdint.h>
#include <cstring>

class VCF {
public:
	uint64_t length;
	int16_t *chromosomes;
	uint64_t *positions;
	char **references;
	char **variants;

private:
	char *filepath;
	FILE *source_file;
	char line_buffer[2048];
	char chromosome_buffer[64];
	char position_buffer[64];
	char reference_buffer[512];
	char variant_buffer[512];
	char dummy_buffer[512];
	uint16_t id_buffer_length;

public:
	~VCF() {
		if (chromosomes) free(chromosomes);
		if (positions) free(positions);
		if (references) {
			for (uint64_t i = 0; i < length; i++) {
				free(references[i]);
			}
			free(references);
		}
		if (variants) {
			for (uint64_t i = 0; i < length; i++) {
				free(variants[i]);
			}
			free(variants);
		}
		free(filepath);
	}

	static VCF *ReadFile(char *filepath, int16_t chromosome);
private:

	VCF(char *filepath) {
		this->filepath = strdup(filepath);
		length = 0;
		chromosomes = NULL;
		positions = NULL;
		references = NULL;
		variants = NULL;
	}

	void InitializeArrays();
	uint64_t ReadChromosomeRowCount(int16_t chromosome);
	void ReadChromosome(int16_t chromosome);
	void AddVariant(uint64_t index, int16_t chromosome, uint64_t position, char *reference, uint16_t ref_len, char *variant, uint16_t var_len);
	void AddSingleVariant(uint64_t index, int16_t chromosome, uint64_t position, char *reference, uint16_t ref_len, char *variant, uint16_t var_len);
	void AddMultiVariant(uint64_t index, int16_t chromosome, uint64_t position, char *reference, uint16_t ref_len, char *variant_buffer, uint16_t var_len);
};

#endif
