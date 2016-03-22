/*
 * Copyright (C) 2015 - 2016  Fabian Klötzl
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

#ifndef QUARTET_H
#define QUARTET_H

#include "graph.h"
#include "matrix.h"

int quartet_root(matrix *distance, tree_root *root);
void quartet_all(matrix *distance, tree_s *baum);
double support(const matrix *distance, const char *types);

// A set of four colors.
enum { SET_D, SET_A, SET_B, SET_C };

typedef struct color_context {
	char *types;
	size_t size;
	char color;
} color_context;

void colorize(tree_node *current, color_context *);
void colorize_dry(tree_node *foo, tree_node *bar, color_context *cctx);

#endif
