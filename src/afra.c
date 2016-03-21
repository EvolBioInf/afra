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

#include <assert.h>
#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "io.h"
#include "matrix.h"
#include "graph.h"
#include "quartet.h"

void consense(char **matrix_names, matrix distance, tree_root root);

int main(int argc, const char *argv[]) {

	if (argc < 2 || (argc == 2 && strcmp(argv[1],"-h") == 0)) {
		fprintf(stderr, "Usage: %s quartet|consense [MATRIX...]\n", argv[0]);
		return 1;
	}

	argv += 1;

	enum {QUARTET,CONSENSE} mode;

	if(strcmp(*argv,"quartet") ==  0){
		mode = QUARTET;
	} else if (strcmp(*argv,"consense") == 0){
		mode = CONSENSE;
	} else {
		errx(1, "invalid mode. Should be one of 'quartet' or 'consense'.");
	}

	argv += 1;

	int firsttime = 1;
	int exit_code = EXIT_SUCCESS;

	for (;; firsttime = 0) {
		FILE *file_ptr;
		const char *file_name;
		if (!*argv) {
			if (!firsttime) exit(exit_code);

			file_ptr = stdin;
			file_name = "stdin";
		} else {
			file_name = *argv++;
			file_ptr = fopen(file_name, "r");
			if (!file_ptr) err(1, "%s", file_name);
		}

		char **matrix_names;
		matrix distance = read_matrix(file_ptr, &matrix_names);

		if (distance.size < 4) {
			errx(1, "this program requires at least four taxa.");
		}

		tree_s tree;
		neighbor_joining(&distance, &tree);

		if( mode == CONSENSE){
			consense(matrix_names, distance, tree.root);
		} else {
			quad_root(&distance, &tree.root);
			newick_sv(&tree.root, matrix_names);
		}

		fclose(file_ptr);
		tree_free(&tree);
		for (size_t i = 0; i < distance.size; i++) {
			free(matrix_names[i]);
		}
		free(matrix_names);
	}

	return exit_code;
}
