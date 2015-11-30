/*
 * Copyright (C) 2015  Fabian Klötzl
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

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
