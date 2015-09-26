#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "graph.h"
#include "matrix.h"

tree_node *node_pool;

int neighbor_joining(matrix *distance, tree_node *out_root) {
	size_t matrix_size = distance->size;
	if (matrix_size < 3) return -2;

	node_pool = malloc(2 * matrix_size * sizeof(tree_node));
	tree_node *empty_node_ptr = &node_pool[matrix_size];
	tree_node **unjoined_nodes = malloc(matrix_size * sizeof(tree_node *));

	size_t n = matrix_size;
	size_t i, j;

	for (i = 0; i < n; i++) {
		node_pool[i] = LEAF(i);
		unjoined_nodes[i] = &node_pool[i];
	}

	double r[matrix_size];

	matrix local_copy;
	int check = matrix_copy(&local_copy, distance);
	assert(check == 0);
	assert(local_copy.size == distance->size);
	assert(memcmp(local_copy.data, distance->data,
	              matrix_size * matrix_size * sizeof(double)) == 0);

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
		    .right_dist = (M(min_i, min_j) - r[min_i] + r[min_j]) / 2.0,
		    .index = -1};

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
	tree_node root = {.left_branch = unjoined_nodes[0],
	                  .right_branch = unjoined_nodes[1],
	                  .extra_branch = unjoined_nodes[2],

	                  .left_dist = (M(0, 1) + M(0, 2) - M(1, 2)) / 2.0,
	                  .right_dist = (M(0, 1) + M(1, 2) - M(0, 2)) / 2.0,
	                  .extra_dist = (M(0, 2) + M(1, 2) - M(0, 1)) / 2.0};

	*empty_node_ptr++ = root;
	*out_root = root;

	free(unjoined_nodes);
	return 0;
}

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

void traverse_ctx(tree_node *current, visitor_ctx *v, void *context) {
	if (v->pre) {
		v->pre(current, context);
	}
	if (current->left_branch) {
		traverse_ctx(current->left_branch, v, context);
	}
	if (v->process) {
		v->process(current, context);
	}
	if (current->right_branch) {
		traverse_ctx(current->right_branch, v, context);
	}
	if (v->post) {
		v->post(current, context);
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
