#include <err.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

#include "matrix.h"

matrix read_matrix(FILE *in, char ***out_matrix_names) {
	size_t matrix_size;
	char **matrix_names;

	int check = fscanf(in, "%zu\n", &matrix_size);
	if (check < 1) goto format_error;

	matrix distance;
	int l = matrix_init(&distance, matrix_size);
	if (l != 0) goto format_error;

	matrix_names = malloc(matrix_size * sizeof(char *));
	if (!matrix_names) goto format_error;

	size_t i, j;
	for (i = 0; i < matrix_size; i++) {
		check = fscanf(in, "%ms ", &matrix_names[i]);
		if (check < 1) goto format_error;
		for (j = 0; j < matrix_size; j++) {
			check = fscanf(in, "%lf ", &MATRIX_CELL(distance, i, j));
			if (check < 1) goto format_error;
		}
	}

	*out_matrix_names = matrix_names;
	return distance;

format_error:
	errx(1, "format error: expected phylip-style matrix");
}
