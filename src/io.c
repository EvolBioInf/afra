#include <err.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

#define BUFFERSIZE 4096

#define M(I, J) (matrix_ptr[(I)*matrix_size + (J)])

static char *buffer = NULL;
int read_matrix(FILE *in, size_t *out_matrix_size, double **out_matrix_ptr,
                char ***out_matrix_names) {
	size_t matrix_size;
	double *matrix_ptr;
	char **matrix_names;

	int check = fscanf(in, "%zu\n", &matrix_size);
	if (check < 1) err(errno, "derp1");

	matrix_ptr = malloc(matrix_size * matrix_size * sizeof(double));
	matrix_names = malloc(matrix_size * sizeof(char *));
	if (!matrix_ptr || !matrix_names) err(errno, "derp2");

	size_t i, j;
	for (i = 0; i < matrix_size; i++) {
		check = fscanf(in, "%ms ", &matrix_names[i]);
		if (check < 1) err(errno, "derp3");
		for (j = 0; j < matrix_size; j++) {
			check = fscanf(in, "%lf ", &M(i, j));
			if (check < 1) err(errno, "derp4");
		}
	}

	*out_matrix_size = matrix_size;
	*out_matrix_ptr = matrix_ptr;
	*out_matrix_names = matrix_names;
	return 0;
}
