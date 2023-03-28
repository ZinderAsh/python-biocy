#ifndef BIOCY_FILE_READER_VCF
#define BIOCY_FILE_READER_VCF

#include <cstdlib>
#include <stdio.h>
#include <stdint.h>
#include <cstring>

class VCF {
public:
	uint64_t length;
	uint32_t *chromosomes;
	uint64_t *positions;
	char **references;
	char **variants;

private:
	char *filepath;
	FILE *source_file;
	char line_buffer[2048];
	char chromosome_buffer[64];
	char position_buffer[64];
	char reference_buffer[64];
	char variant_buffer[64];
	char dummy_buffer[64];
	uint8_t id_buffer_length;

public:
	~VCF() {
		if (chromosomes) free(chromosomes);
		if (positions) free(positions);
		if (references) free(references);
		if (variants) free(variants);
		free(filepath);
	}

	static VCF *ReadFile(char *filepath, uint32_t chromosome);
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
	uint64_t ReadChromosomeRowCount(uint32_t chromosome);
	void ReadChromosome(uint32_t chromosome);
	
	/*
	void ReadNodeCountAndIDRange();
	void InitializeArrays();
	void ReadSequences();
	void ReadReferenceNodes();
	void CountEdges();
	void InitializeEdges();
	void ReadEdges();
	int64_t ReadNextID();
	uint32_t GetMappedID(uint32_t);
	*/
};

#endif
