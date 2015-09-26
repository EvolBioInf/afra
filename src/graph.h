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

#define LEAF(I) ((struct tree_node){.index = (I)})
#define BRANCH(...) ((struct tree_node){__VA_ARGS__})

int neighbor_joining(matrix *, tree_node *);
int matrix_from_tree(size_t matrix_size, tree_root root);

void newick(tree_node *root);

typedef void (*tree_node_processor)(tree_node *);
typedef struct visitor { tree_node_processor pre, process, post; } visitor;
void traverse(tree_node *current, visitor *v);

typedef void (*tree_node_processor_context)(tree_node *, void *);
typedef struct visitor_ctx {
	tree_node_processor_context pre, process, post;
} visitor_ctx;
void traverse_ctx(tree_node *current, visitor_ctx *v, void *);
