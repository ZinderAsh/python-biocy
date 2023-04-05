#include "FASTA.hpp"
#include <iostream>

FASTA *FASTA::ReadFile(char *filepath) {
	FASTA *fasta = new FASTA(filepath);
	fasta->source_file = fopen(filepath, "rb");

	return fasta;
}

bool FASTA::GoToChromosome(int16_t chromosome) {
	fseek(source_file, 0, SEEK_SET);
	int16_t current_chromosome = -1;
	char c;
	while (current_chromosome != chromosome) {
		while ((c = fgetc(source_file)) != '>') {
			if (c == EOF) {
				printf("Failed to find chromosome #%d in file.\n", chromosome);
				return false;
			}
		}
		fgets(line_buffer, 2047, source_file);
		if (strtol(line_buffer, NULL, 10) == chromosome) {
			current_chromosome = chromosome;
		}
	}
	return (current_chromosome != -1);
}

char *FASTA::ReadNext(uint32_t count) {
	if (buffer_len < count) {
		buffer_len = count;
		buffer = (char *) realloc(buffer, sizeof(char) * (buffer_len + 1));
	}
	
	uint32_t buffer_pos = 0;
	char c;
	while (buffer_pos < count) {
		c = fgetc(source_file);
		if (c == EOF || c == '>') break;
		if (c == 'A' || c == 'a' || c == 'C' || c == 'c' || c == 'G' || c == 'g' ||
				c == 'T' || c == 't' || c == 'N' || c == 'n') {
			buffer[buffer_pos++] = c;
		}
	}
	if (buffer_pos == 0) return NULL;

	buffer[buffer_pos] = '\0';

	return buffer;
}
