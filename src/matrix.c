#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "matrix.h"

void matrix_free(matrix *mx) {
	if (!mx) return;
	free(mx->data);
	*mx = (struct matrix){0, NULL};
}

int matrix_init(matrix *mx, size_t size) {
	if (!mx || !size) return -1;
	mx->data = malloc(size * size * sizeof(double));
	if (!mx->data) return errno;
	mx->size = size;

	return 0;
}

int matrix_copy(matrix *dest, matrix *src) {
	int check = matrix_init(dest, src->size);
	if (check == 0) {
		memcpy(dest->data, src->data, src->size * src->size * sizeof(double));
	}
	return check;
}
