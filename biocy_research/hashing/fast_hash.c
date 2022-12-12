#include <stdio.h>
#include <stdlib.h>

int main() {
	const char in[] = {
		'A', 'a',
		'C', 'c',
		'G', 'g',
		'T', 't'
	};
	const char out[] = {
		0, 0,
		1, 1,
		2, 2,
		3, 3
	};
	int num_tests = 8;

	int min_add = 999;
	int min_xor = 999;
	int min_shift1 = 999;
	int min_shift2 = 999;
	int min_operations = 999;

	char tmp;

	for (int add = 0; add < 256; add++) {
		for (int xor = 0; xor < 256; xor++) {
			for (int shift1 = 0; shift1 < 8; shift1++) {
				for (int shift2 = 0; shift2 < 8; shift2++) {
					int passed = 1;
					for (int i = 0; i < num_tests; i++) {
						tmp = (in[i] ^ xor) + add;
						tmp = (((tmp >> shift1) & 1) << 1) | ((tmp >> shift2) & 1);
						if (tmp != out[i]) {
							passed = 0;
							break;
						}
					}
					if (passed) {
						int operations = 0;
						if (add != 0) operations++;
						if (xor != 0) operations++;
						if (shift1 != 0) operations++;
						if (shift2 != 0) operations++;
						if (operations < min_operations) {
							min_operations = operations;
							min_add = add;
							min_xor = xor;
							min_shift1 = shift1;
							min_shift2 = shift2;
						}
					}
				}
			}
		}
	}

	printf("Formula Found! Operations: %d\n", min_operations);
	printf("x = (in ^ %d) + %d\n", min_xor, min_add);
	printf("hash = ((((x >> %d) & 1) << 1)) | ((x >> %d) & 1)\n\n", min_shift1, min_shift2);
	return 0;
}
