#include <assert.h>
#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io.h"
#include "matrix.h"
#include "graph.h"

int quadruple_stats(matrix *distance, tree_node *root);

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
		newick(&root);

		quadruple_stats(&distance, &root);

		fclose(file_ptr);
	}

	return exit_code;
}

enum { SET_D, SET_A, SET_B, SET_C };
char *types, type;

void consume(tree_node *current);

int quadruple_stats(matrix *distance, tree_node *root) {
	// left branch
	size_t size = distance->size;

	types = malloc(size);
	memset(types, SET_D, size);

	type = SET_A;
	consume(root->left_branch->left_branch);
	type = SET_B;
	consume(root->left_branch->right_branch);
	type = SET_C;
	consume(root->right_branch);
	// D = not A, B, C;

	size_t set_counts[4] = {0};
	for (size_t i = 0; i < size; i++) {
		set_counts[(int)types[i]]++;
	}

	printf("%zu %zu %zu %zu\n", set_counts[SET_A], set_counts[SET_B],
	       set_counts[SET_C], set_counts[SET_D]);

	size_t non_supporting_counter = 0;

	size_t A = 0, B, C, D;
	for (; A < size; A++) {
		for (; types[A] != SET_A && A < size; A++)
			;
		if (A >= size) {
			break;
		}

		for (B = 0; B < size; B++) {
			for (; types[B] != SET_B && B < size; B++)
				;
			if (B >= size) {
				break;
			}

			for (C = 0; C < size; C++) {
				for (; types[C] != SET_C && C < size; C++)
					;
				if (C >= size) {
					break;
				}

				for (D = 0; D < size; D++) {
					for (; types[D] != SET_D && D < size; D++)
						;
					if (D >= size) {
						break;
					}

#define M(I, J) (MATRIX_CELL(*distance, I, J))
					printf("%zu %zu %zu %zu\n", A, B, C, D);

					double D_abcd = M(A, B) + M(C, D);
					if (M(A, C) + M(B, D) < D_abcd ||
					    M(A, D) + M(B, C) < D_abcd) {
						non_supporting_counter++;
					}
				}
			}
		}
	}

	printf("%zu\n", non_supporting_counter);
}

void consume_process(tree_node *current) {
	if (!current->left_branch) {
		types[current->index] = type;
	}
}

void consume(tree_node *current) {
	if (!current) return;
	visitor v = {.pre = NULL, .process = consume_process, .post = NULL};

	traverse(current, &v);
}
