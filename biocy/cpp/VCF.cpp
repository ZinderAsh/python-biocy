#include "VCF.hpp"
#include <iostream>

VCF *VCF::ReadFile(char *filepath, int16_t chromosome) {
	VCF *vcf = new VCF(filepath);
	vcf->source_file = fopen(filepath, "rb");

	vcf->length = vcf->ReadChromosomeRowCount(chromosome);

	vcf->InitializeArrays();

	vcf->ReadChromosome(chromosome);
	
	fclose(vcf->source_file);

	return vcf;
}

void VCF::InitializeArrays() {
	chromosomes = (int16_t *) malloc(sizeof(int16_t) * length);
	positions   = (uint64_t *) malloc(sizeof(uint64_t) * length);
	references  = (char **) malloc(sizeof(char *) * length);
	variants    = (char **) malloc(sizeof(char *) * length);

	memset(chromosomes, 0, sizeof(int16_t) * length);
	memset(positions, 0, sizeof(uint64_t) * length);
	memset(references, 0, sizeof(char *) * length);
	memset(variants, 0, sizeof(char *) * length);
}

bool is_structural_variant(char *ref, uint8_t ref_len, char *var, uint8_t var_len) {
	for (uint8_t i = 0; i < ref_len; i++) {
		char base = ref[i];
		if (i == 0 && base == '.') break;
		if (base != 'A' && base != 'a' && base != 'C' && base != 'c' &&
		    base != 'G' && base != 'g' && base != 'T' && base != 't' &&
		    base != 'N' && base != 'n' && base != ',') {
			return true;
		}
	}
	for (uint8_t i = 0; i < var_len; i++) {
		char base = var[i];
		if (i == 0 && base == '.') break;
		if (base != 'A' && base != 'a' && base != 'C' && base != 'c' &&
		    base != 'G' && base != 'g' && base != 'T' && base != 't' &&
		    base != 'N' && base != 'n' && base != ',') {
			return true;
		}
	}
	return false;
}

uint64_t VCF::ReadChromosomeRowCount(int16_t chromosome) {
	uint64_t count = 0;
	uint64_t ignored = 0;

	fseek(source_file, 0, SEEK_SET);
	while (fgets(line_buffer, 2048, source_file) != NULL) {
		if (line_buffer[0] != '#') {
			sscanf(line_buffer, "%s\t%s\t%s\t%s\t%s",
				chromosome_buffer, position_buffer, dummy_buffer, reference_buffer, variant_buffer);
			int16_t row_chromosome = strtoull(chromosome_buffer, NULL, 10);
			if (chromosome == -1 || chromosome == row_chromosome) {
				uint8_t reference_length = strlen(reference_buffer);
				uint8_t variant_length = strlen(variant_buffer);
				if (is_structural_variant(reference_buffer, reference_length,
							  variant_buffer, variant_length)) {
					ignored++;
					continue;
				}
				count++;
			};
		}
	}

	printf("Structural Variants ignored: %lu\n", ignored);
	printf("Count: %lu\n", count);

	return count;
}

void VCF::AddMultiVariant(uint64_t index, int16_t chromosome, uint64_t position, char *reference, uint8_t ref_len, char *variant_buffer, uint8_t var_len) {
	uint8_t variant_count = 1;
	for (uint8_t i = 0; i < var_len; i++) {
		if (variant_buffer[i] == ',') variant_count++;
	}

	char **all_variants = (char **) malloc(sizeof(char *) * variant_count);

	char *variant = strtok(variant_buffer, ",");

	variant_count = 0;

	while (variant != NULL) {
		all_variants[variant_count] = variant;
		variant = strtok(NULL, ",");
		variant_count++;
	}

	uint8_t variant_start_pos = 0;
	bool keep_going = true;
	while (keep_going) {
		if (variant_start_pos >= ref_len) {
			keep_going = false;
			break;
		}
		for (uint8_t i = 0; i < variant_count; i++) {
			if (variant_start_pos >= strlen(all_variants[i])) {
				keep_going = false;
				break;
			}
		}
		if (!keep_going) break;
		char ref_base = reference[variant_start_pos];
		for (uint8_t i = 0; i < variant_count; i++) {
			if (all_variants[i][variant_start_pos] != ref_base) {
				keep_going = false;
				break;
			}
		}
		if (!keep_going) break;
		variant_start_pos++;
	}

	if (variant_start_pos > 0) {
		for (uint8_t i = variant_start_pos; i < ref_len; i++) {
			reference[i - variant_start_pos] = reference[i];
		}
		ref_len -= variant_start_pos;
		reference[ref_len] = '\0';
		for (uint8_t i = 0; i < variant_count; i++) {
			uint8_t this_var_len = strlen(all_variants[i]);
			for (uint8_t j = variant_start_pos; j < this_var_len; j++) {
				all_variants[i][j - variant_start_pos] = all_variants[i][j];
			}
			all_variants[i][this_var_len - variant_start_pos] = '\0';
		}
		position += variant_start_pos;
	}

	uint8_t variant_buffer_len = 0;
	for (uint8_t i = 0; i < variant_count; i++) {
		uint8_t this_var_len = strlen(all_variants[i]);
		for (uint8_t j = 0; j < this_var_len; j++) {
			variant_buffer[variant_buffer_len++] = all_variants[i][j];
		}
		variant_buffer[variant_buffer_len++] = ((i == variant_count - 1) ? '\0' : ',');
	}

	AddVariant(index, chromosome, position, reference, ref_len, variant_buffer, variant_buffer_len);

	free(all_variants);
}

void VCF::AddSingleVariant(uint64_t index, int16_t chromosome, uint64_t position, char *reference, uint8_t ref_len, char *variant, uint8_t var_len) {
	uint8_t variant_start_pos = 0;
	while (variant_start_pos < ref_len && variant_start_pos < var_len &&
	       reference[variant_start_pos] == variant[variant_start_pos]) {
		variant_start_pos++;
	}
	//variant_start_pos = 0;
	if (variant_start_pos > 0) {
		for (uint8_t i = variant_start_pos; i < ref_len; i++) {
			reference[i - variant_start_pos] = reference[i];
		}
		ref_len -= variant_start_pos;
		reference[ref_len] = '\0';
		for (uint8_t i = variant_start_pos; i < var_len; i++) {
			variant[i - variant_start_pos] = variant[i];
		}
		var_len -= variant_start_pos;
		variant[var_len] = '\0';
		position += variant_start_pos;
	}

	AddVariant(index, chromosome, position, reference, ref_len, variant, var_len);
}

void VCF::AddVariant(uint64_t index, int16_t chromosome, uint64_t position, char *reference, uint8_t ref_len, char *variant, uint8_t var_len) {
	if (reference[0] == '.') {
		references[index] = (char *) malloc(sizeof(char));
		references[index][0] = '\0';
	} else {
		references[index] = (char *) malloc(sizeof(char) * (ref_len + 1));
		strncpy(references[index], reference, ref_len);
		references[index][ref_len] = '\0';
	}
	if (variant[0] == '.') {
		variants[index] = (char *) malloc(sizeof(char));
		variants[index][0] = '\0';
	} else {
		variants[index] = (char *) malloc(sizeof(char) * (var_len + 1));
		strncpy(variants[index], variant, var_len);
		variants[index][var_len] = '\0';
	}

	chromosomes[index] = chromosome;
	positions[index] = position;
}

void VCF::ReadChromosome(int16_t chromosome) {
	uint64_t count = 0;

	fseek(source_file, 0, SEEK_SET);
	while (fgets(line_buffer, 2048, source_file) != NULL) {
		if (line_buffer[0] != '#') {
			sscanf(line_buffer, "%s\t%s\t%s\t%s\t%s",
				chromosome_buffer, position_buffer, dummy_buffer, reference_buffer, variant_buffer);
			int16_t row_chromosome = strtoull(chromosome_buffer, NULL, 10);
			if (chromosome == -1 || chromosome == row_chromosome) {
				uint8_t reference_length = strlen(reference_buffer);
				uint8_t variant_length = strlen(variant_buffer);
				if (is_structural_variant(reference_buffer, reference_length,
							  variant_buffer, variant_length)) {
					continue;
				}
				uint64_t position = strtoull(position_buffer, NULL, 10);

				bool multiple_variant = (strstr(variant_buffer, ",") != NULL);

				if (multiple_variant) {
					AddMultiVariant(count, row_chromosome, position,
							reference_buffer, reference_length,
							variant_buffer, variant_length);
				} else {
					AddSingleVariant(count, row_chromosome, position,
							 reference_buffer, reference_length,
							 variant_buffer, variant_length);
				}

				count++;
			}
		}
	}
}
