#include <assert.h>
#include <err.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "io.h"
#include "matrix.h"

typedef struct tree_node {
	struct tree_node *left_branch, *right_branch;
	double left_dist, right_dist;
} tree_node;

typedef struct tree_root {
	tree_node *left_branch, *right_branch, *extra_branch;
	double left_dist, right_dist, extra_dist;
} tree_root;

#define LEAF() ((struct tree_node){0})
#define BRANCH(...) ((struct tree_node){__VA_ARGS__})

#define M(I, J) (matrix_ptr[(I)*matrix_size + (J)])

int neighbor_joining(size_t matrix_size, double *matrix_ptr);
int matrix_from_tree(size_t matrix_size, tree_root root);

int main(int argc, const char *argv[]) {

	if (argc == 2 && argv[1][0] == '-' && argv[1][1] == 'h') {
		fprintf(stderr, "Usage: %s [FASTA...]\n", argv[0]);
		return 1;
	}

	argv += 1;

	int firsttime = 1;
	int exit_code = EXIT_SUCCESS;

	for (;; firsttime = 0) {
		FILE *file_ptr;
		const char *file_name;
		if (!*argv) {
			if (!firsttime) exit(exit_code);

			file_ptr = stdin;
			file_name = "stdin";
		} else {
			file_name = *argv++;
			file_ptr = fopen(file_name, "r");
			if (!file_ptr) err(1, "%s", file_name);
		}

		size_t matrix_size;
		double *matrix_ptr;
		char **matrix_names;

		int l;
		if ((l = read_matrix(file_ptr, &matrix_size, &matrix_ptr,
		                     &matrix_names)) != 0) {
			exit_code = EXIT_FAILURE;
			goto fail;
		}

		// print matrix
		// printf("size: %zu\n", matrix_size);

		neighbor_joining(matrix_size, matrix_ptr);

	fail:
		fclose(file_ptr);
	}

	return exit_code;
}


tree_node *node_pool;

void stringify2(size_t n, tree_node *node) {
	if (!node->left_branch) {
		printf("%zu", node - node_pool);
		return;
	}
	printf("(");
	if (node->left_branch) {
		stringify2(n, node->left_branch);
		printf(":%lf", node->left_dist);
	}
	printf(",");
	if (node->right_branch) {
		stringify2(n, node->right_branch);
		printf(":%lf", node->right_dist);
	}
	printf(")");
}

void stringify(size_t n, tree_root root) {
	printf("(");

	// root: left
	stringify2(n, root.left_branch);
	printf(":%lf,", root.left_dist);

	// root: right
	stringify2(n, root.right_branch);
	printf(":%lf,", root.right_dist);

	// root: extra
	stringify2(n, root.extra_branch);
	printf(":%lf", root.extra_dist);

	printf(");\n");
}

int neighbor_joining(size_t matrix_size, double *matrix_ptr) {
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
	tree_root root;

	root.left_branch = unjoined_nodes[0];
	root.right_branch = unjoined_nodes[1];
	root.extra_branch = unjoined_nodes[2];

	root.left_dist = (M(0, 1) + M(0, 2) - M(1, 2)) / 2.0;
	root.right_dist = (M(0, 1) + M(1, 2) - M(0, 2)) / 2.0;
	root.extra_dist = (M(0, 2) + M(1, 2) - M(0, 1)) / 2.0;

	stringify(999, root);

	matrix_from_tree(matrix_size, root);

	free(unjoined_nodes);
	return -1;
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
			printf("%lf ", INDUCED(i,j));
		}
		printf("\n");
	}

	matrix_free(&induced);
}
