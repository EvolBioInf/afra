/*
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
#include <sys/types.h>

#include "global.h"
#include "quartet.h"

double support(const matrix *distance, const char *types) {
	const size_t size = distance->size;

	size_t non_supporting_counter = 0;
	size_t quartet_counter = 0;

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

					quartet_counter++;

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

	// printf("%zu of %zu\n", non_supporting_counter, quartet_counter);
	return 1 - ((double)non_supporting_counter / quartet_counter);
}

void quartet_left(tree_node *current, matrix *distance) {
	if (!current->left_branch || !current->left_branch->left_branch) return;
	color_context cctx = {.size = distance->size,
	                      .types = malloc(distance->size)};
	CHECK_MALLOC(cctx.types);

	colorize_dry(current->left_branch, current->right_branch, &cctx);

	double d = support(distance, cctx.types);
	current->left_support = d;

	free(cctx.types);
}

void quartet_right(tree_node *current, matrix *distance) {
	if (!current->left_branch || !current->right_branch->left_branch) return;
	color_context cctx = {.size = distance->size,
	                      .types = malloc(distance->size)};
	CHECK_MALLOC(cctx.types);

	colorize_dry(current->right_branch, current->left_branch, &cctx);

	double d = support(distance, cctx.types);
	current->right_support = d;

	free(cctx.types);
}

void quartet_node(tree_node *current, void *ctx) {
	if (!current->left_branch) return;

	// left branch
	quartet_left(current, (matrix *)ctx);

	// right branch
	quartet_right(current, (matrix *)ctx);
}

void quartet_all(matrix *distance, tree_s *baum) {
	// iterate over all nodes
	size_t size = distance->size;
	tree_node *inner_nodes = baum->pool + size;

#pragma omp parallel for schedule(dynamic) num_threads(THREADS)
	for (size_t i = 0; i < size - 2; i++) {
		quartet_node(&inner_nodes[i], distance);
	}

	tree_root *root = &baum->root;
	quartet_node(&root->as_tree_node, distance);

	if (root->extra_branch->left_branch) {
		// Support Value for Root→Extra
		color_context cctx = {.size = distance->size,
		                      .types = malloc(distance->size)};
		CHECK_MALLOC(cctx.types);

		colorize_dry(root->extra_branch, root->left_branch, &cctx);

		double d = support(distance, cctx.types);
		root->extra_support = d;

		free(cctx.types);
	}
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

	traverse_all(current, &v, cctx);
}

/** Colorize according to the following scheme.
 *
 *  A -left--             --y---- C, bar
 *           \           /
 *          foo --x-- current
 *           /           \
 *  B -right-             -extra- D
 *
 * current is the node of the caller. It has two branches pointing down (left
 * and right). The extra branch points to the parent. foo is supposed to be one
 * of the down pointers and bar the other.
 *
 * @param foo - The node connecting the subsets A and B.
 * @param bar - The node of subset C.
 */
void colorize_dry(tree_node *foo, tree_node *bar, color_context *cctx) {
	memset(cctx->types, SET_D, cctx->size);

	cctx->color = SET_A;
	colorize(foo->left_branch, cctx);

	cctx->color = SET_B;
	colorize(foo->right_branch, cctx);

	cctx->color = SET_C;
	colorize(bar, cctx);

	// D is not A, B or C
}
