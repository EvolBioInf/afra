#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "graph.h"
#include "matrix.h"

tree_node *node_pool;

void newick(tree_node *root);

int neighbor_joining(matrix *distance) {
	size_t matrix_size = distance->size;
	if (matrix_size < 3) return -2;

	node_pool = malloc(2 * matrix_size * sizeof(tree_node));
	tree_node *empty_node_ptr = &node_pool[matrix_size];
	tree_node **unjoined_nodes = malloc(matrix_size * sizeof(tree_node *));

	size_t n = matrix_size;
	size_t i, j;

	for (i = 0; i < n; i++) {
		node_pool[i] = LEAF();
		unjoined_nodes[i] = &node_pool[i];
	}

	double r[matrix_size];

	matrix local_copy;
	matrix_copy(&local_copy, distance);

#define M(I, J) (MATRIX_CELL(local_copy, I, J))

	while (n > 3) {
		for (i = 0; i < n; i++) {
			double rr = 0.0;
			for (j = 0; j < n; j++) {
				if (i == j) assert(M(i, j) == 0.0);
				assert(M(i, j) == M(j, i));

				rr += M(i, j);
			}
			r[i] = rr / (double)(n - 2);
		}

		size_t min_i = 0, min_j = 1;
		double min_value = M(0, 1) - r[0] - r[1];

		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				if (i == j) continue;

				double value = M(i, j) - r[i] - r[j];
				if (value < min_value) {
					min_i = i;
					min_j = j;
					min_value = value;
				}
			}
		}

		// force i < j
		if (min_j < min_i) {
			size_t temp = min_i;
			min_i = min_j;
			min_j = temp;
		}

		tree_node branch = {
		    .left_branch = unjoined_nodes[min_i],
		    .right_branch = unjoined_nodes[min_j],
		    .left_dist = (M(min_i, min_j) + r[min_i] - r[min_j]) / 2.0,
		    .right_dist = (M(min_i, min_j) - r[min_i] + r[min_j]) / 2.0};

		*empty_node_ptr++ = branch;
		unjoined_nodes[min_i] = empty_node_ptr - 1;
		unjoined_nodes[min_j] = unjoined_nodes[n - 1];

		double row_k[matrix_size];
		double M_ij = M(min_i, min_j);

		for (size_t m = 0; m < n; m++) {
			if (m == min_i || m == min_j) continue;

			row_k[m] = (M(min_i, m) + M(min_j, m) - M_ij) / 2.0;
			// if( row_k[m] < 0) row_k[m] = 0;
		}

		// row_k[min_i] and row_k[min_j] are undefined!
		row_k[min_i] = 0.0;
		row_k[min_j] = row_k[n - 1];

		memcpy(&M(min_i, 0), row_k, matrix_size * sizeof(double));
		memcpy(&M(min_j, 0), &M((n - 1), 0), matrix_size * sizeof(double));

		M(min_i, min_i) = M(min_j, min_j) = 0.0;

		for (i = 0; i < n; i++) {
			M(i, min_i) = M(min_i, i);
		}

		for (i = 0; i < n; i++) {
			M(i, min_j) = M(min_j, i);
		}

		n--;
	}

	// join three remaining nodes
	tree_node root = {0};

	root.left_branch = unjoined_nodes[0];
	root.right_branch = unjoined_nodes[1];
	root.extra_branch = unjoined_nodes[2];

	root.left_dist = (M(0, 1) + M(0, 2) - M(1, 2)) / 2.0;
	root.right_dist = (M(0, 1) + M(1, 2) - M(0, 2)) / 2.0;
	root.extra_dist = (M(0, 2) + M(1, 2) - M(0, 1)) / 2.0;

	*empty_node_ptr++ = root;

	newick(&root);

	/*
	// arbitrary root:
	size_t min_i = 0;
	size_t min_j = 1;
	tree_node branch = {.left_branch = unjoined_nodes[min_i],
	                    .right_branch = unjoined_nodes[min_j],
	                    .left_dist = (M(min_i, min_j)) / 2.0,
	                    .right_dist = (M(min_i, min_j)) / 2.0};

	*empty_node_ptr++ = branch;
	unjoined_nodes[min_i] = empty_node_ptr - 1;*/

	// matrix_from_tree(matrix_size, root);

	free(unjoined_nodes);
	return -1;
}

typedef void (*tree_node_processor)(tree_node *);

typedef struct visitor { tree_node_processor pre, process, post; } visitor;

void traverse(tree_node *current, visitor *v) {
	if (v->pre) {
		v->pre(current);
	}
	if (current->left_branch) {
		traverse(current->left_branch, v);
	}
	if (v->process) {
		v->process(current);
	}
	if (current->right_branch) {
		traverse(current->right_branch, v);
	}
	if (v->post) {
		v->post(current);
	}
}

void newick_pre(tree_node *current) {
	if (current->left_branch) printf("(");
}

void newick_post(tree_node *current) {
	if (current->left_branch) printf(":%lf)", current->right_dist);
}

void newick_process(tree_node *current) {
	if (current->left_branch) {
		printf(":%lf,", current->left_dist);
	} else {
		printf("%zu", current - node_pool);
	}
}

void newick(tree_node *root) {
	visitor v = {
	    .pre = newick_pre, .process = newick_process, .post = newick_post};

	printf("(");
	traverse(root->left_branch, &v);
	printf(":%lf,", root->left_dist);
	traverse(root->right_branch, &v);
	printf(":%lf,", root->right_dist);
	traverse(root->extra_branch, &v);
	printf(":%lf);\n", root->extra_dist);
}

#define INDUCED(I, J) MATRIX_CELL(induced, I, J)

int matrix_from_tree(size_t matrix_size, tree_root root) {
	matrix induced;
	matrix_init(&induced, matrix_size * 2 - 1);

	fprintf(stderr, "%zu\n", induced.size);

	size_t i, j;
	for (i = 0; i < induced.size; i++) {
		for (j = 0; j < induced.size; j++) {
			INDUCED(i, j) = i == j ? 0 : 99999.; // zero? save_inf?
		}
	}

	// infer induced distances from tree
	tree_node *derp = &node_pool[matrix_size * 2 - 3];
	while (derp > node_pool) {
		i = derp - node_pool;
		if (derp->left_branch) {
			j = derp->left_branch - node_pool;
			INDUCED(i, j) = INDUCED(j, i) = derp->left_dist;
		}

		if (derp->right_branch) {
			j = derp->right_branch - node_pool;
			INDUCED(i, j) = INDUCED(j, i) = derp->right_dist;
		}
		derp--;
	}

	size_t A = root.left_branch - node_pool;
	size_t B = root.right_branch - node_pool;
	size_t C = root.extra_branch - node_pool;

	INDUCED(A, B) = INDUCED(B, A) = root.left_dist + root.right_dist;
	INDUCED(A, C) = INDUCED(C, A) = root.left_dist + root.extra_dist;
	INDUCED(C, B) = INDUCED(B, C) = root.right_dist + root.extra_dist;

	// Floyd-Warshall
	size_t k;

	for (k = 0; k < induced.size; k++) {
		for (i = 0; i < induced.size; i++) {
			for (j = 0; j < induced.size; j++) {
				double d = INDUCED(k, i) + INDUCED(k, j);
				if (d < INDUCED(i, j)) {
					INDUCED(i, j) = INDUCED(j, i) = d;
				}
			}
		}
	}

	for (i = 0; i < matrix_size; i++) {
		for (j = 0; j < matrix_size; j++) {
			printf("%lf ", INDUCED(i, j));
		}
		printf("\n");
	}

	// midpoint root
	size_t max_i = 0, max_j = 0;
	double max_value = 0.0;

	for (i = 0; i < matrix_size; i++) {
		for (j = 0; j < matrix_size; j++) {
			double d = INDUCED(i, j);
			if (d > max_value) {
				max_value = d;
				max_i = i;
				max_j = j;
			}
		}
	}

	fprintf(stderr, "most distant pair: %zu %zu %lf\n", max_i, max_j,
	        max_value);

	// mroot:
	// find the halfway point.
	// add node there and reroot.

	matrix_free(&induced);
}
