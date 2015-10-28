#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "quadruple.h"

enum { SET_D, SET_A, SET_B, SET_C };

double support(const matrix *distance, const char *types) {
	const size_t size = distance->size;

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
	// printf("%lf\n", d);

	free(cctx.types);
}

void quad_node(tree_node *current, void *ctx) {
	if (!current->left_branch) return;

	// left branch
	quad_left(current, (matrix *)ctx);

	// right branch
	quad_right(current, (matrix *)ctx);
}

int quad_root(matrix *distance, tree_root *root) {

	visitor_ctx v = {.pre = NULL, .process = quad_node, .post = NULL};

	traverse_ctx(root, &v, distance);
	traverse_ctx(root->extra_branch, &v, distance);

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
		// printf("%lf\n", d);

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
	visitor_ctx v = {.pre = NULL, .process = colorize_process, .post = NULL};

	traverse_ctx(current, &v, cctx);
}
