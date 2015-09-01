#ifndef _MATRIX_H_
#define _MATRIX_H_ 1

typedef struct matrix {
	size_t size;
	double *data;
} matrix;

int matrix_init(matrix *, size_t);
void matrix_free(matrix *);
int matrix_copy(matrix *dest, matrix *src);

#define MATRIX_CELL(MATRIX, I, J) ((MATRIX).data[(I) * (MATRIX).size + (J)])

#endif
