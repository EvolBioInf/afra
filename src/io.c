#include <err.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

#include "matrix.h"

matrix read_matrix(FILE *in, char ***out_matrix_names) {
	size_t matrix_size;
	char **matrix_names;

	int check = fscanf(in, "%zu\n", &matrix_size);
	if (check < 1) err(errno, "derp1");

	matrix distance;
	int l = matrix_init(&distance, matrix_size);
	if (l != 0) err(errno, "huh?");

	matrix_names = malloc(matrix_size * sizeof(char *));
	if (!matrix_names) err(errno, "derp2");

	size_t i, j;
	for (i = 0; i < matrix_size; i++) {
		check = fscanf(in, "%ms ", &matrix_names[i]);
		if (check < 1) err(errno, "derp3");
		for (j = 0; j < matrix_size; j++) {
			check = fscanf(in, "%lf ", &MATRIX_CELL(distance, i, j));
			if (check < 1) err(errno, "derp4");
		}
	}

	*out_matrix_names = matrix_names;
	return distance;
}
