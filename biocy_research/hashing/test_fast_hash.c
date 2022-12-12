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

	for (int i = 0; i < num_tests; i++) {
		res = in[i] + 9;
		res = ((res >> 3) & 2) | ((res >> 2) & 1);
		if (res != out[i]) {
			printf("Failed for %c! Got %d, expected %d\n", in[i], res, out[i]);
			return 1;
		}
	}

	printf("Success!\n");

	return 0;
}
