/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Angus Lepper (angus.lepper@ed.ac.uk)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_CELL_CACHE_H
#define SWIFT_CELL_CACHE_H

#include <stdbool.h>

typedef float cell_cache_float;

struct cell_cache_data {
  cell_cache_float *parts_xs[3];
};

void cell_cache_init(int nr_threads, int nr_slots, int nr_ways);
void cell_cache_invalidate(void);

struct runner;
struct cell_cache_data *cell_cache_single(struct runner *r, struct cell *c,
    int sid);
void cell_cache_pair(struct runner* r, struct cell_cache_data **di,
    struct cell *ci, struct cell_cache_data **dj, struct cell *cj, int sid);

#endif
