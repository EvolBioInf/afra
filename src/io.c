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

#include <err.h>
#include <errno.h>
#include <stdlib.h>
#include <stdio.h>

#include "global.h"
#include "matrix.h"

matrix read_matrix(FILE *in, char ***out_matrix_names) {
	size_t matrix_size;

	int check = fscanf(in, "%zu\n", &matrix_size);
	if (check < 1) goto format_error;

	matrix distance;
	int l = matrix_init(&distance, matrix_size);
	if (l != 0) goto format_error;

	size_t i, j;
	for (i = 0; i < matrix_size; i++) {
		check = fscanf(in, "%ms ", &distance.names[i]);
		if (check < 1) goto format_error;
		for (j = 0; j < matrix_size; j++) {
			check = fscanf(in, "%lf ", &MATRIX_CELL(distance, i, j));
			if (check < 1) goto format_error;
		}
	}

	*out_matrix_names = distance.names;
	return distance;

format_error:
	errx(1, "format error: expected phylip-style matrix");
}
