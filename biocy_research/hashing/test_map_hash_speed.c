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

	char res;
	
	char map[255];
	for (int i = 0; i < num_tests; i++) {
		map[in[i]] = out[i];
	}

	for (long long a = 0; a < 500000000; a++) {
		for (int i = 0; i < num_tests; i++) {
			res = map[in[i]];
		}
	}

	printf("Success!\n");

	return 0;
}
