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

		printf("\nConsensus tree program, version 3.695\n\n");

		print_species(matrix_names, distance.size);

		printf("\n\n\nSets included in the consensus tree\n\n"
			"Set (species in order)     How many times out of  100.00\n\n");

		quad_root(&distance, &root);

		printf("\n\nSets NOT included in consensus tree: NONE.\n\n");

		newick_sv(&root, matrix_names);

		fclose(file_ptr);

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
		printf("%s", ((char**)ctx)[current->index]);
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

void newick_sv(tree_node *root, char **names) {
	ctx_visitor v = {.pre = newick_sv_pre,
	                 .process = newick_sv_process,
	                 .post = newick_sv_post};

	printf("(");
	ctx_traverse(root->left_branch, &v, names);
	newick_sv_process(root, names);

	ctx_traverse(root->right_branch, &v, names);
	if (root->right_branch && root->right_branch->right_branch) {
		printf("%d:%lf,", (int)(root->right_support * 100), root->right_dist);
	} else {
		printf(":%lf,", root->right_dist);
	}

	ctx_traverse(root->extra_branch, &v, names);
	if (root->extra_branch && root->extra_branch->left_branch) {
		printf("%d:%lf)", (int)(root->extra_support * 100), root->extra_dist);
	} else {
		printf(":%lf)", root->extra_dist);
	}
	printf(";\n");
}

enum { SET_D, SET_A, SET_B, SET_C };

double support(const matrix *distance, const char *types) {
	const size_t size = distance->size;

	/*
	size_t set_counts[4] = {0};

	for (size_t i = 0; i < size; i++) {
	    set_counts[(int)types[i]]++;
	}

	printf("%zu %zu %zu %zu\n", set_counts[SET_A], set_counts[SET_B],
	       set_counts[SET_C], set_counts[SET_D]); */

	size_t non_supporting_counter = 0;
	size_t quadruple_counter = 0;

	size_t A = 0, B, C, D;
	for (; A < size; A++) {
		if (types[A] != SET_A) continue;

		for (B = 0; B < size; B++) {
			if (types[B] != SET_B) continue;

			for (C = 0; C < size; C++) {
				if (types[C] != SET_C) continue;

				for (D = 0; D < size; D++) {
					if (types[D] != SET_D) continue;

 #define M(I, J) (MATRIX_CELL(*distance, I, J))

					quadruple_counter++;

					double D_abcd = M(A, B) + M(C, D);
					if (((M(A, C) + M(B, D)) < D_abcd) ||
					    ((M(A, D) + M(B, C)) < D_abcd)) {
						// printf("%zu %zu %zu %zu\n", A, B, C, D);
						non_supporting_counter++;
					}
				}
			}
		}
	}

	// printf("%zu of %zu\n", non_supporting_counter, quadruple_counter);
	return 1 - ((double)non_supporting_counter / quadruple_counter);
}

typedef struct color_context { char *types, color; } color_context;

void colorize(tree_node *current, color_context *);

void quad_left(tree_node *current, matrix *distance) {
	if (!current->left_branch || !current->left_branch->left_branch) return;
	color_context cctx;

	cctx.types = malloc(distance->size);
	memset(cctx.types, SET_D, distance->size);

	cctx.color = SET_A;
	colorize(current->left_branch->left_branch, &cctx);

	cctx.color = SET_B;
	colorize(current->left_branch->right_branch, &cctx);

	cctx.color = SET_C;
	colorize(current->right_branch, &cctx);
	// D = not A, B, C;

	double d = support(distance, cctx.types);
	current->left_support = d;
	// printf("%lf\n", d);

	for(size_t i=0; i<distance->size; i++){
		printf("%c", (cctx.types[i] == SET_A || cctx.types[i] == SET_B) ? '*' : '.');
	}
	printf("                     %2.1lf\n", d*100);

	free(cctx.types);
}

void quad_right(tree_node *current, matrix *distance) {
	if (!current->left_branch || !current->right_branch->left_branch) return;
	color_context cctx;

	cctx.types = malloc(distance->size);
	memset(cctx.types, SET_D, distance->size);

	cctx.color = SET_A;
	colorize(current->right_branch->left_branch, &cctx);

	cctx.color = SET_B;
	colorize(current->right_branch->right_branch, &cctx);

	cctx.color = SET_C;
	colorize(current->left_branch, &cctx);
	// D = not A, B, C;

	double d = support(distance, cctx.types);
	current->right_support = d;

	for(size_t i=0; i<distance->size; i++){
		printf("%c", (cctx.types[i] == SET_A || cctx.types[i] == SET_B) ? '*' : '.');
	}
	printf("                     %2.1lf\n", d*100);

	free(cctx.types);
}

void quad_node(tree_node *current, void *ctx) {
	if (!current->left_branch) return;

	// left branch
	quad_left(current, (matrix *)ctx);

	// right branch
	quad_right(current, (matrix *)ctx);
}

int quad_root(matrix *distance, tree_node *root) {

	ctx_visitor v = {.pre = NULL, .process = quad_node, .post = NULL};

	ctx_traverse(root, &v, distance);
	ctx_traverse(root->extra_branch, &v, distance);

	if (root->extra_branch->left_branch) {
		// Support Value for Root→Extra
		color_context cctx;

		cctx.types = malloc(distance->size);
		memset(cctx.types, SET_D, distance->size);

		cctx.color = SET_A;
		colorize(root->extra_branch->left_branch, &cctx);

		cctx.color = SET_B;
		colorize(root->extra_branch->right_branch, &cctx);

		cctx.color = SET_C;
		colorize(root->left_branch, &cctx);
		// D = not A, B, C;

		double d = support(distance, cctx.types);
		root->extra_support = d;
		
		for(size_t i=0; i<distance->size; i++){
			printf("%c", (cctx.types[i] == SET_A || cctx.types[i] == SET_B) ? '*' : '.');
		}
		printf("                     %2.1lf\n", d*100);

		free(cctx.types);
	}

	return 0;
}

void colorize_process(tree_node *current, void *vctx) {
	color_context *ctx = (color_context *)vctx;

	if (!current->left_branch) {
		ctx->types[current->index] = ctx->color;
	}
}

void colorize(tree_node *current, color_context *cctx) {
	if (!current) return;
	ctx_visitor v = {.pre = NULL, .process = colorize_process, .post = NULL};

	ctx_traverse(current, &v, cctx);
}

void print_species(char **names, size_t n){
	printf("Species in order:\n\n");
	for(size_t i=0; i< n; i++){
		printf("  %zu. %s\n",i+1, names[i] );
	}
}

void print_sets(tree_root *root, color_context *cctx){
}

void print_sets_node(tree_node *current, void *ctx){
	size_t size;
	color_context cctx;
	cctx.types = malloc(size);
	memset(cctx.types, SET_D, size);

	cctx.color = SET_A;
	colorize(current->left_branch, &cctx);

	for(size_t i=0; i<size; i++){
		printf("%c", cctx.types[i] == SET_A ? '*' : '.');
	}
	printf("                     %lf\n", current->left_dist);

	free(cctx.types);
}
