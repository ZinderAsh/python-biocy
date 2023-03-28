#include "VCF.hpp"
#include <iostream>

VCF *VCF::ReadFile(char *filepath, uint32_t chromosome) {
	VCF *vcf = new VCF(filepath);
	vcf->source_file = fopen(filepath, "rb");

	vcf->length = vcf->ReadChromosomeRowCount(chromosome);

	vcf->InitializeArrays();

	vcf->ReadChromosome(chromosome);
	
	fclose(vcf->source_file);

	return vcf;
}

void VCF::InitializeArrays() {
	uint32_t *chromosomes;
	uint64_t *positions;
	char **references;
	char **variants;
	chromosomes = (uint32_t *) malloc(sizeof(uint32_t) * length);
	positions   = (uint64_t *) malloc(sizeof(uint64_t) * length);
	references  = (char **) malloc(sizeof(char *) * length);
	variants    = (char **) malloc(sizeof(char *) * length);

	memset(chromosomes, 0, sizeof(uint32_t) * length);
	memset(positions, 0, sizeof(uint64_t) * length);
	memset(references, 0, sizeof(char *) * length);
	memset(variants, 0, sizeof(char *) * length);
}

uint64_t VCF::ReadChromosomeRowCount(uint32_t chromosome) {
	uint64_t count = 0;
	uint64_t ignored = 0;

	fseek(source_file, 0, SEEK_SET);
	while (fgets(line_buffer, 2048, source_file) != NULL) {
		if (line_buffer[0] != '#') {
			sscanf(line_buffer, "%s\t%s\t%s\t%s\t%s",
				chromosome_buffer, position_buffer, dummy_buffer, reference_buffer, variant_buffer);
			uint32_t row_chromosome = strtoull(chromosome_buffer, NULL, 10);
			if (chromosome == -1 || chromosome == row_chromosome) {
				uint8_t variant_length = strlen(variant_buffer);
				bool structural_variant = false;
				for (uint8_t i = 0; i < variant_length; i++) {
					char base = variant_buffer[i];
					if (base != 'A' && base != 'a' && base != 'C' && base != 'c' && base != 'G' && base != 'g' &&
						base != 'T' && base != 't' && base != 'N' && base != 'n') {
						ignored++;
						structural_variant = true;
						break;
					}
				}
				if (!structural_variant) count++;
			};
		}
	}

	printf("Structural Variants ignored: %lu\n", ignored);
	printf("Count: %lu\n", count);

	return count;
}

void VCF::ReadChromosome(uint32_t chromosome) {
	uint64_t count = 0;

	fseek(source_file, 0, SEEK_SET);
	while (fgets(line_buffer, 2048, source_file) != NULL) {
		if (line_buffer[0] != '#') {
			sscanf(line_buffer, "%s\t%s\t%s\t%s\t%s",
				chromosome_buffer, position_buffer, dummy_buffer, reference_buffer, variant_buffer);
			uint32_t row_chromosome = strtoull(chromosome_buffer, NULL, 10);
			if (chromosome == -1 || chromosome == row_chromosome) {
				uint8_t variant_length = strlen(variant_buffer);
				bool structural_variant = false;
				for (uint8_t i = 0; i < variant_length; i++) {
					char base = variant_buffer[i];
					if (base != 'A' && base != 'a' && base != 'C' && base != 'c' && base != 'G' && base != 'g' &&
						base != 'T' && base != 't' && base != 'N' && base != 'n') {
						structural_variant = true;
						break;
					}
				}
				if (structural_variant) continue;
				uint64_t position = strtoull(position_buffer, NULL, 10);
				count++;
			}
		}
	}
}
