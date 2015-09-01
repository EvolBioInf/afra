#include <stdio.h>

#include "matrix.h"

typedef struct tree_node {
	struct tree_node *left_branch, *right_branch, *extra_branch;
	double left_dist, right_dist, extra_dist;
	double left_support, right_support, extra_support;
	ssize_t index;
} tree_node;

typedef struct tree_root {
	tree_node *left_branch, *right_branch, *extra_branch;
	double left_dist, right_dist, extra_dist;
} tree_root;

#define LEAF() ((struct tree_node){0})
#define BRANCH(...) ((struct tree_node){__VA_ARGS__})

int neighbor_joining(matrix *distance);
int matrix_from_tree(size_t matrix_size, tree_root root);
