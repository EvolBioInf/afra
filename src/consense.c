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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "matrix.h"
#include "global.h"
#include "graph.h"
#include "quartet.h"

void print_species(char **names, size_t n);
int set_root(matrix *distance, tree_root *root);

void consense(char **matrix_names, matrix distance, tree_root root) {
	printf("\nConsensus tree program, version 3.695\n\n");

	print_species(matrix_names, distance.size);

	printf("\n\n\nSets included in the consensus tree\n\n"
	       "Set (species in order)     How many times out of  100.00\n\n");

	set_root(&distance, &root);

	printf("\n\nSets NOT included in consensus tree: NONE.\n\n");

	newick_sv(&root, matrix_names);
}

void set_left(tree_node *current, matrix *distance) {
	if (!current->left_branch || !current->left_branch->left_branch) return;
	color_context cctx = {.size = distance->size,
	                      .types = malloc(distance->size),
	                      .color = SET_A};
	CHECK_MALLOC(cctx.types);

	memset(cctx.types, SET_D, distance->size);
	colorize(current->left_branch, &cctx);

	double d = current->left_support;

	for (size_t i = 0; i < distance->size; i++) {
		printf("%c", cctx.types[i] == SET_A ? '*' : '.');
	}
	printf("                     %2.1lf\n", d * 100);

	free(cctx.types);
}

void set_right(tree_node *current, matrix *distance) {
	if (!current->left_branch || !current->right_branch->left_branch) return;
	color_context cctx = {.size = distance->size,
	                      .types = malloc(distance->size),
	                      .color = SET_A};
	CHECK_MALLOC(cctx.types);

	memset(cctx.types, SET_D, distance->size);
	colorize(current->right_branch, &cctx);

	double d = current->right_support;

	for (size_t i = 0; i < distance->size; i++) {
		printf("%c", cctx.types[i] == SET_A ? '*' : '.');
	}
	printf("                     %2.1lf\n", d * 100);

	free(cctx.types);
}

void set_node(tree_node *current, void *ctx) {
	if (!current->left_branch) return;

	// left branch
	set_left(current, (matrix *)ctx);

	// right branch
	set_right(current, (matrix *)ctx);
}

int set_root(matrix *distance, tree_root *root) {

	visitor_ctx v = {.pre = NULL, .process = set_node, .post = NULL};

	traverse_all(&root->as_tree_node, &v, distance);
	traverse_all(root->extra_branch, &v, distance);

	if (root->extra_branch->left_branch) {
		// Support Value for Root→Extra
		color_context cctx = {.size = distance->size,
		                      .types = malloc(distance->size),
		                      .color = SET_A};
		CHECK_MALLOC(cctx.types);

		memset(cctx.types, SET_D, distance->size);
		colorize(root->extra_branch, &cctx);

		double d = root->extra_support;

		for (size_t i = 0; i < distance->size; i++) {
			printf("%c", cctx.types[i] == SET_A ? '*' : '.');
		}
		printf("                     %2.1lf\n", d * 100);

		free(cctx.types);
	}

	return 0;
}

void print_species(char **names, size_t n) {
	printf("Species in order:\n\n");
	for (size_t i = 0; i < n; i++) {
		printf("  %zu. %s\n", i + 1, names[i]);
	}
}
