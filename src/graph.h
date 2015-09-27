#include <stdio.h>

#include "matrix.h"

typedef struct tree_node {
	struct tree_node *left_branch, *right_branch;
	double left_dist, right_dist;
	double left_support, right_support;
	ssize_t index;
} tree_node;

typedef struct tree_root {
	struct tree_node;
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
