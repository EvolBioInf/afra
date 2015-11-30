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

#ifndef GRAPH_H
#define GRAPH_H

#include <stdio.h>

#include "matrix.h"

typedef struct tree_node {
	struct tree_node *left_branch, *right_branch;
	double left_dist, right_dist;
	double left_support, right_support;
	ssize_t index;
} tree_node;

typedef struct tree_root {
	union {
		struct tree_node;
		struct tree_node as_tree_node;
	};
	tree_node *extra_branch;
	double extra_dist;
	double extra_support;
} tree_root;

#define LEAF(I) ((struct tree_node){.index = (I)})
#define BRANCH(...) ((struct tree_node){__VA_ARGS__})

typedef struct tree_s {
	size_t size;
	tree_node *pool;
	tree_root root;
} tree_s;

int tree_init(tree_s *baum, size_t size);
void tree_free(tree_s *baum);

int neighbor_joining(matrix *distance, tree_s *out_tree);

typedef void (*tree_node_processor_context)(tree_node *, void *);
typedef struct visitor_ctx {
	tree_node_processor_context pre, process, post;
} visitor_ctx;
void traverse_ctx(tree_node *current, visitor_ctx *v, void *);

void newick_sv(tree_root *, char **);

#endif
