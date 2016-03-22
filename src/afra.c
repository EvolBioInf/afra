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

#include <assert.h>
#include <err.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#include "config.h"
#include "io.h"
#include "matrix.h"
#include "graph.h"
#include "quartet.h"

void usage(int);
void version(void);

void consense(char **matrix_names, matrix distance, tree_root root);

int main(int argc, char *argv[]) {

	const struct option long_options[] = {
	    {"version", no_argument, NULL, 'V'},
	    {"mode", required_argument, NULL, 'm'},
	    {"help", no_argument, NULL, 'h'},
	    {0, 0, 0, 0}};

	enum { QUARTET, CONSENSE } mode = QUARTET;

	while (1) {
		int c = getopt_long(argc, argv, "Vhm:", long_options, NULL);
		if (c == -1) {
			break;
		}

		switch (c) {
		case 'V':
			version();
		case 'h':
			usage(EXIT_SUCCESS);
		case 'm':
			if (strcmp(optarg, "quartet") == 0) {
				mode = QUARTET;
			} else if (strcmp(optarg, "consense") == 0) {
				mode = CONSENSE;
			} else {
				errx(1,
				     "invalid mode. Should be one of 'quartet' or 'consense'.");
			}
			break;
		case '?': /* intentional fall-through */
		default:
			usage(EXIT_FAILURE);
		}
	}

	argv += optind;

	int firsttime = 1;

	for (;; firsttime = 0) {
		FILE *file_ptr;
		const char *file_name;
		if (!*argv) {
			if (!firsttime) break;

			file_ptr = stdin;
			file_name = "stdin";
		} else {
			file_name = *argv++;
			file_ptr = fopen(file_name, "r");
			if (!file_ptr) err(1, "%s", file_name);
		}

		matrix distance = read_matrix(file_ptr);

		if (distance.size < 4) {
			errx(1, "this program requires at least four taxa.");
		}

		tree_s tree;
		neighbor_joining(&distance, &tree);

		if (mode == CONSENSE) {
			consense(distance.names, distance, tree.root);
		} else {
			quad_all(&distance, &tree);
			newick_sv(&tree.root, distance.names);
		}

		fclose(file_ptr);
		tree_free(&tree);
		matrix_free(&distance);
	}

	return EXIT_SUCCESS;
}

void usage(int exit_code) {
	static const char *str = {
	    "Usage: afra [-Vh] [-m quartet|consense] [MATRIX...]\n"
	    "\tMATRIX... can be any sequence of matrices in PHYLIP format. If no "
	    "files are supplied, stdin is used instead.\n"
	    "Options:\n"
	    "  -m, --mode <quartet|consense>\n"
	    "                    Analysis mode; default: quartet\n"
	    "  -h, --help        Display this help and exit\n"
	    "  -V, --version     Output version information\n"};

	printf("%s", str);
	exit(exit_code);
}

void version(void) {
	static const char *str = {
	    "afra " VERSION "\n"
	    "Copyright (C) 2015 - 2016 Fabian Klötzl\n"
	    "Licenses GPLv3+: GNU GPL version 3 or later "
	    "<http://gnu.org/licenses/gpl.html>\n"
	    "This is free software: you are free to change and redistribute it.\n"
	    "There is NO WARRANTY, to the extent permittet by law.\n"};

	printf("%s", str);
	exit(EXIT_SUCCESS);
}
