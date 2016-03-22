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
#include <errno.h>
#include <getopt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "config.h"
#include "io.h"
#include "matrix.h"
#include "graph.h"
#include "quartet.h"

int THREADS = 1;

void usage(int);
void version(void);

void consense(char **matrix_names, matrix distance, tree_root root);

int main(int argc, char *argv[]) {

	const struct option long_options[] = {
	    {"version", no_argument, NULL, 'V'},
	    {"mode", required_argument, NULL, 'm'},
	    {"help", no_argument, NULL, 'h'},
	    {"threads", required_argument, NULL, 't'},
	    {0, 0, 0, 0}};

#ifdef _OPENMP
	// Use all available processors by default.
	THREADS = omp_get_num_procs();
#endif

	enum { QUARTET, CONSENSE } mode = QUARTET;

	while (1) {
		int c = getopt_long(argc, argv, "Vhm:t:", long_options, NULL);
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
		case 't': {
#ifdef _OPENMP
			errno = 0;
			char *end;
			long unsigned int threads = strtoul(optarg, &end, 10);

			if (errno || end == optarg || *end != '\0') {
				warnx("Expected a number for -t argument, but '%s' was "
				      "given. Ignoring -t argument.",
				      optarg);
				break;
			}

			if (threads > (long unsigned int)omp_get_num_procs()) {
				warnx("The number of threads to be used, is greater then the "
				      "number of available processors; Ignoring -t %lu "
				      "argument.",
				      threads);
				break;
			}

			THREADS = threads;
#else
			warnx("This version of afra was built without OpenMP and thus "
			      "does not support multi threading. Ignoring -t argument.");
#endif
			break;
		}

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
		quartet_all(&distance, &tree);

		if (mode == CONSENSE) {
			consense(distance.names, distance, tree.root);
		} else {
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
	    "Usage: afra [-Vh] [-t INT] [-m quartet|consense] [MATRIX...]\n"
	    "\tMATRIX... can be any sequence of matrices in PHYLIP format. If no "
	    "files are supplied, stdin is used instead.\n"
	    "Options:\n"
	    "  -m, --mode <quartet|consense>\n"
	    "                    Analysis mode; default: quartet\n"
	    "  -t, --threads int Number of threads; by default all processors are "
	    "used.\n"
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
