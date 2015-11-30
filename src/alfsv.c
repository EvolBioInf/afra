#include <assert.h>
#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io.h"
#include "matrix.h"
#include "graph.h"

int quad_root(matrix *distance, tree_node *root);
void newick_sv(tree_node *, char**);

void print_species(char **names, size_t n);
void consense( char **matrix_names, matrix distance, tree_node root );

int main(int argc, const char *argv[]) {

	if (argc == 2 && argv[1][0] == '-' && argv[1][1] == 'h') {
		fprintf(stderr, "Usage: %s [FASTA...]\n", argv[0]);
		return 1;
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

		tree_node root;
		neighbor_joining(&distance, &root);

		consense( matrix_names, distance, root );

		fclose(file_ptr);

		free(matrix_names);
	}

	return exit_code;
}
