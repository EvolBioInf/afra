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
	double left_support, right_support, extra_support;
} tree_root;

#define LEAF(I) ((struct tree_node){.index = (I)})
#define BRANCH(...) ((struct tree_node){__VA_ARGS__})

typedef struct tree {
	size_t size;
	tree_node *pool;
	tree_root root;
} tree;

int neighbor_joining(matrix *distance, tree_node *out_root, tree *out_tree) ;

typedef void (*tree_node_processor_context)(tree_node *, void *);
typedef struct visitor_ctx {
	tree_node_processor_context pre, process, post;
} visitor_ctx;
void traverse_ctx(tree_node *current, visitor_ctx *v, void *);
