#ifndef BIOCY_FILE_READER_FASTA
#define BIOCY_FILE_READER_FASTA

#include <cstdlib>
#include <stdio.h>
#include <stdint.h>
#include <cstring>

class FASTA {
public:
	char *buffer;
	uint32_t buffer_len;
private:
	char *filepath;
	FILE *source_file;
	char line_buffer[2048];
public:
	~FASTA() {
		if (buffer) free(buffer);
		fclose(source_file);
		free(filepath);
	}

	static FASTA *ReadFile(char *filepath);
	void GoToStart();
	bool GoToChromosome(int16_t chromosome);
	char *ReadNext(uint32_t count);
private:

	FASTA(char *filepath) {
		this->filepath = strdup(filepath);
		buffer_len = 32;
		buffer = (char *) malloc(sizeof(char) * (buffer_len + 1));
	}
};

#endif
