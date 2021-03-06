/** @file This file contains some global definitions and thus should be kept as
 * small as possible.
 *
 * Copyright (C) 2016  Fabian Klötzl
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

#pragma once

#define CHECK_MALLOC(PTR)                                                      \
	do {                                                                       \
		if ((PTR) == NULL) {                                                   \
			err(errno, "Out of memory");                                       \
		}                                                                      \
	} while (0);

extern int THREADS;
