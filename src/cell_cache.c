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

/* Standard headers */
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>

/* Local headers */
#include "atomic.h"
#include "cell.h"
#include "cell_cache.h"
#include "error.h"
#include "runner.h"
#include "space.h"

struct cell_cache_slot {
  struct cell *cell;
  int sid, age;

  int capacity, count;
  void* data;
};

static struct cell_cache_data *cell_cache_view;
static bool *cell_cache_view_off;

static struct cell_cache_slot *cell_cache;
static int cell_cache_nr_threads, cell_cache_nr_slots, cell_cache_nr_ways;

#ifdef CELL_CACHE_STATISTICS
size_t cell_cache_hit_count = 0, cell_cache_miss_count = 0;
int cell_cache_nr_sets;
size_t *cell_cache_set_count;
#endif

void cell_cache_init(int nr_threads, int nr_slots, int nr_ways) {
  int i, j, k;

  if (nr_ways < 2) {
    // DOPAIR1 and DOPAIR2 expect to cache a pair of cells.
    error("cell_cache_nr_ways must be at least 2.");
  }

  if (nr_slots < 2 || nr_slots % nr_ways != 0) {
    // Expect enough slots to divide into complete n-way sets.
    error("cell_cache_nr_slots must be a multiple of cell_cache_nr_ways.");
  }

  int last_set_bit_only = nr_ways & ~(nr_ways-1);
  if (last_set_bit_only != nr_ways) {
    // cell_cache_single expects a power of two to mask out the 'way' hash bits.
    error("cell_cache_nr_ways must be a power of two.");
  }

  if ((cell_cache_view = malloc(2*nr_threads * sizeof *cell_cache_view))
      == NULL) {
    error("Failed to allocate cell cache view.");
  }

  if ((cell_cache_view_off = calloc(nr_threads, sizeof *cell_cache_view_off))
      == NULL) {
    error("Failed to allocate cell cache view.");
  }

#ifdef CELL_CACHE_STATISTICS
  cell_cache_nr_sets = nr_slots / nr_ways;
  cell_cache_set_count = calloc(cell_cache_nr_sets, sizeof *cell_cache_set_count);
#endif

  cell_cache_nr_threads = nr_threads;
  cell_cache_nr_slots = nr_slots;
  cell_cache_nr_ways = nr_ways;

  if ((cell_cache = malloc(nr_threads*nr_slots * sizeof *cell_cache)) == NULL) {
    error("Failed to allocate cell cache.");
  }

  for (i = 0; i < nr_threads; ++i) {
    struct cell_cache_slot *slots = &cell_cache[i * nr_slots];
    for (j = 0; j < nr_slots; ++j) {
      slots[j].cell = NULL;   // Will not match a valid cell.
      slots[j].age = INT_MAX; // (Joint) best candidate for replacement.
      slots[j].capacity = 0;  // Will be expanded: NULL is safe to free.
      slots[j].data = NULL;
    }
  }
}

void cell_cache_invalidate(void) {
  int i;

  for (i = 0; i < cell_cache_nr_threads*cell_cache_nr_slots; ++i) {
    cell_cache[i].cell = NULL;   // Will not match a valid cell.
    cell_cache[i].age = INT_MAX; // (Joint) best candidate for replacement.
  }
}

struct cell_cache_data *
cell_cache_single(struct runner *r, struct cell *c, int sid) {
  int i, j, best_idx, best_age = -1;
  struct cell_cache_slot *slots = &cell_cache[r->id*cell_cache_nr_slots], *best;

  /* Find a slot. */
  unsigned hash = (unsigned)((uintptr_t)c >> 6) & 0x0FFFFFFFu;
  hash |= (unsigned char)sid << 28;

  // TODO: Different hash and/or a fallback from _mm_crc32_u32
  hash = _mm_crc32_u32(0xFFFFFFFFu, hash);
  hash &= 0x100000000ull - cell_cache_nr_ways;

  /* % vs. & avoids strict constraint, but uniform distribution over
   * associative sets only if cell_cache_nr_slots is a power of two i.e. may
   * observe performance degradation otherwise.
   */
  int begin = hash % cell_cache_nr_slots;
  for (i = begin; i < begin + cell_cache_nr_ways; ++i) {
    if (slots[i].cell == c && slots[i].sid == sid) {
#ifdef CELL_CACHE_STATISTICS
      atomic_inc(&cell_cache_hit_count);
#endif
      best_idx = i;
      best_age = 0; // not (necessarily) best->age
      break;
    } else if (slots[i].age > best_age) {
      best_idx = i;
      best_age = slots[i].age; // at least 1
    }
  }

#ifdef CELL_CACHE_STATISTICS
  atomic_inc(&cell_cache_set_count[begin / cell_cache_nr_ways]);
#endif

  /* Age the less-recently-used slots. */
  best = &slots[best_idx];
  if (best->age > 0) {
    for (i = begin; i < begin + cell_cache_nr_ways; ++i) {
      slots[i].age += slots[i].age < INT_MAX ? 1 : 0;
    }

    best->age = 0;
  }

  /* Age the less-recently-used slots. */
  best = &slots[best_idx];
  if (best->age > 0) {
    for (i = begin; i < begin + cell_cache_nr_ways; ++i) {
      slots[i].age += slots[i].age < INT_MAX ? 1 : 0;
    }

    best->age = 0;
  }

  int idx = 2*r->id + cell_cache_view_off[r->id];
  struct cell_cache_data *view = &cell_cache_view[idx];
  cell_cache_view_off[r->id] ^= 1;

  /* Miss? */
  if (best_age > 0) {
#ifdef CELL_CACHE_STATISTICS
    atomic_inc(&cell_cache_miss_count);
#endif
    /* Initialise, with increased capacity if necessary. */
    best->cell = c;
    best->sid = sid;

    best->count = c->count;
    if (best->capacity < c->count) {
      size_t inter_field_space = c->count % VEC_SIZE;
      best->capacity = c->count + inter_field_space;

      free(best->data);

      static const size_t field_count = 3;
      size_t size = field_count * best->capacity * sizeof(cell_cache_float);
      if (posix_memalign(&best->data, 64, size) != 0) {
        error("Failed to expand cell cache slot capacity.");
      }
    }

    view->parts_xs[0] = (cell_cache_float*)best->data;
    view->parts_xs[1] = (cell_cache_float*)best->data + best->capacity;
    view->parts_xs[2] = (cell_cache_float*)best->data + 2*best->capacity;

    if (sid >= 0) {
      struct entry *sort = &c->sort[sid * (c->count + 1)];

      /* Experiment (on COSMA5, with Intel 2015) supports up-front vs. as-we-go
       * prefetch, maybe because the latter inhibits autovectorisation. Also,
       * unexpectedly, no benefit was found in particle prefetch.
       */
      for (i = 0; i < c->count; i += 8) {
        __builtin_prefetch(&sort[i]);
      }

      for (i = 0; i < c->count; ++i) {                                          // c->count x...
        for (j = 0; j < 3; ++j) {                                               //   3x...
          view->parts_xs[j][i] = c->parts[sort[i].i].x[j];                      //     4B in sort[i].i, 8B in x[j], 8B out parts_xs[j][i]
        }
      }
    } else {
      // Just the storage order.
      for (i = 0; i < c->count; ++i) {
        for (j = 0; j < 3; ++j) {
          view->parts_xs[j][i] = c->parts[i].x[j];
        }
      }
    }
  } else {
    view->parts_xs[0] = (cell_cache_float*)best->data;
    view->parts_xs[1] = (cell_cache_float*)best->data + best->capacity;
    view->parts_xs[2] = (cell_cache_float*)best->data + 2*best->capacity;
  }

  return view;
}

void
cell_cache_pair(struct runner *r, struct cell_cache_data **di, struct cell *ci,
                struct cell_cache_data **dj, struct cell *cj, int sid) {
  /* With a low number of slots, or large associative sets, there's a risk the
   * first cell evicts a cache of the second, forcing an otherwise-unnecessary
   * reload.
   */
  *di = cell_cache_single(r, ci, sid);
  *dj = cell_cache_single(r, cj, sid);
}
