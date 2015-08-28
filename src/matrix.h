#ifndef _MATRIX_H_
#define _MATRIX_H_ 1

typedef struct matrix {
	size_t size;
	double *data;
} matrix;

int matrix_init( matrix*, size_t);
void matrix_free( matrix*);

#define MATRIX_CELL(MATRIX,I,J) ((MATRIX).data[(I)*(MATRIX).size+(J)])

#endif
