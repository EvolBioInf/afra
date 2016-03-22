/** @file This module implements a distance matrix.
 *
 * Copyright (C) 2015 - 2016  Fabian Klötzl
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

#include <err.h>
#include <errno.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "matrix.h"

/** @brief Free all memory aquired by a distance matrix.
 *
 * @param mx - The matrix to destroy. All fields will be zeroed.
 */
void matrix_free(matrix *mx) {
	if (!mx) return;
	if (mx->names) {
		for (size_t i = 0; i < mx->size; i++) {
			free(mx->names[i]);
		}
		free(mx->names);
	}
	free(mx->data);
	*mx = (struct matrix){0, NULL, NULL};
}

/** @brief Create a new distance matrix. Will allocate enough space for the data
 * and matrix names array. The names themselves have to be stored somewhere
 * else (i.e. heap).
 *
 * @param mx - Space for the new matrix.
 * @param size - The matrix' size.
 * @returns 0 on success.
 */
int matrix_init(matrix *mx, size_t size) {
	if (!mx || !size) return -1;
	mx->size = size;
	// potential integer overflow.
	mx->data = malloc(size * size * sizeof(double));
	mx->names = malloc(size * sizeof(char *));
	CHECK_MALLOC(mx->data);
	CHECK_MALLOC(mx->names);

	return 0;
}

/** @brief Creates a copy of a matrix. Does *not* copy the names, only data.
 *
 * @param dest - The destination matrix.
 * @param src - The source matrix.
 * @returns 0 on success.
 */
int matrix_copy(matrix *dest, const matrix *src) {
	if (!dest || !src || !src->size) return -1;
	size_t size = src->size;

	matrix_init(dest, size);
	memcpy(dest->data, src->data, size * size * sizeof(double));
	memset(dest->names, 0, size * sizeof(char *));

	return 0;
}
