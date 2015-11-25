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

void newick_sv(tree_root *, char **);

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

		if (distance.size < 4) {
			errx(1, "this program requires at least four taxa.");
		}

		tree_s tree;
		neighbor_joining(&distance, &tree);
		// newick(&root);

		quad_root(&distance, &tree.root);
		newick_sv(&tree.root, matrix_names);

		fclose(file_ptr);
		tree_free(&tree);
		for (size_t i = 0; i < distance.size; i++) {
			free(matrix_names[i]);
		}
		free(matrix_names);
	}

	return exit_code;
}

void newick_sv_pre(tree_node *current, void *ctx) {
	if (current->left_branch) {
		printf("(");
	}
}

void newick_sv_process(tree_node *current, void *ctx) {
	if (current->left_branch) {
		if (current->left_branch->left_branch) {
			printf("%d:%lf,", (int)(current->left_support * 100),
			       current->left_dist);
		} else {
			printf(":%lf,", current->left_dist);
		}
	} else {
		printf("%s", ((char **)ctx)[current->index]);
	}
}

void newick_sv_post(tree_node *current, void *ctx) {
	if (!current->right_branch) return;
	if (current->right_branch->right_branch) {
		printf("%d:%lf)", (int)(current->right_support * 100),
		       current->right_dist);
	} else {
		printf(":%lf)", current->right_dist);
	}
}

void newick_sv(tree_root *root, char **names) {
	visitor_ctx v = {.pre = newick_sv_pre,
	                 .process = newick_sv_process,
	                 .post = newick_sv_post};

	printf("(");
	traverse_ctx(root->left_branch, &v, names);
	newick_sv_process(&(root->as_tree_node), names);

	traverse_ctx(root->right_branch, &v, names);
	if (root->right_branch && root->right_branch->right_branch) {
		printf("%d:%lf,", (int)(root->right_support * 100), root->right_dist);
	} else {
		printf(":%lf,", root->right_dist);
	}

	traverse_ctx(root->extra_branch, &v, names);
	if (root->extra_branch && root->extra_branch->left_branch) {
		printf("%d:%lf)", (int)(root->extra_support * 100), root->extra_dist);
	} else {
		printf(":%lf)", root->extra_dist);
	}
	printf(";\n");
}
