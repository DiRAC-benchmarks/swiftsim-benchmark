/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Peter W. Draper (p.w.draper@durham.ac.uk)
 *                    Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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

/**
 *  @file partition.c
 *  @brief file of various techniques for partitioning and repartitioning
 *  a grid of cells into geometrically connected regions and distributing
 *  these around a number of MPI nodes.
 *
 *  Currently supported partitioning types: grid, vectorise and METIS.
 */

/* Config parameters. */
#include "../config.h"

/* Standard headers. */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>
#include <values.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
/* METIS headers only used when MPI is also available. */
#ifdef HAVE_METIS
#include <metis.h>
#endif
#endif

/* Local headers. */
#include "const.h"
#include "debug.h"
#include "error.h"
#include "partition.h"
#include "space.h"

/* Maximum weight used for METIS. */
#define metis_maxweight 10000.0f

/* Simple descriptions of initial partition types for reports. */
const char *initial_partition_name[] = {
    "gridded cells",                 "vectorized point associated cells",
    "METIS particle weighted cells", "METIS unweighted cells"};

/* Simple descriptions of repartition types for reports. */
const char *repartition_name[] = {
    "no",
    "METIS edge and vertex time weighted cells",
    "METIS particle count vertex weighted cells",
    "METIS time edge weighted cells",
    "METIS particle count vertex and time edge cells"};

/* Local functions, if needed. */
static int check_complete(struct space *s, int verbose, int nregions);

/*  Vectorisation support */
/*  ===================== */

#if defined(WITH_MPI)
/**
 *  @brief Pick a number of cell positions from a vectorised list.
 *
 *  Vectorise the cell space and pick positions in it for the number of
 *  expected regions using a single step. Vectorisation is guaranteed
 *  to work, providing there are more cells than regions.
 *
 *  @param s the space.
 *  @param nregions the number of regions
 *  @param samplecells the list of sample cell positions, size of 3*nregions
 */
static void pick_vector(struct space *s, int nregions, int *samplecells) {

  /* Get length of space and divide up. */
  int length = s->cdim[0] * s->cdim[1] * s->cdim[2];
  if (nregions > length) {
    error("Too few cells (%d) for this number of regions (%d)", length,
          nregions);
  }

  int step = length / nregions;
  int n = 0;
  int m = 0;
  int l = 0;

  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        if (n == 0 && l < nregions) {
          samplecells[m++] = i;
          samplecells[m++] = j;
          samplecells[m++] = k;
          l++;
        }
        n++;
        if (n == step) n = 0;
      }
    }
  }
}
#endif

#if defined(WITH_MPI)
/**
 * @brief Partition the space.
 *
 * Using the sample positions as seeds pick cells that are geometrically
 * closest and apply the partition to the space.
 */
static void split_vector(struct space *s, int nregions, int *samplecells) {
  int n = 0;
  for (int i = 0; i < s->cdim[0]; i++) {
    for (int j = 0; j < s->cdim[1]; j++) {
      for (int k = 0; k < s->cdim[2]; k++) {
        int select = -1;
        float rsqmax = FLT_MAX;
        int m = 0;
        for (int l = 0; l < nregions; l++) {
          float dx = samplecells[m++] - i;
          float dy = samplecells[m++] - j;
          float dz = samplecells[m++] - k;
          float rsq = (dx * dx + dy * dy + dz * dz);
          if (rsq < rsqmax) {
            rsqmax = rsq;
            select = l;
          }
        }
        s->cells[n++].nodeID = select;
      }
    }
  }
}
#endif

/* METIS support
 * =============
 *
 * METIS partitions using a multi-level k-way scheme. We support using this in
 * a unweighted scheme, which works well and seems to be guaranteed, and a
 * weighted by the number of particles scheme. Note METIS is optional.
 *
 * Repartitioning is based on METIS and uses weights determined from the times
 * that cell tasks have taken. These weight the graph edges and vertices, or
 * just the edges, with vertex weights from the particle counts or none.
 */

#if defined(WITH_MPI) && defined(HAVE_METIS)
/**
 * @brief Fill the METIS xadj and adjncy arrays defining the graph of cells
 *        in a space.
 *
 * See the METIS manual if you want to understand this format. The cell graph
 * consists of all nodes as vertices with edges as the connections to all
 * neighbours, so we have 26 per vertex.
 *
 * @param s the space of cells.
 * @param adjncy the METIS adjncy array to fill, must be of size
 *               26 * the number of cells in the space.
 * @param xadj the METIS xadj array to fill, must be of size
 *             number of cells in space + 1. NULL for not used.
 */
static void graph_init_metis(struct space *s, idx_t *adjncy, idx_t *xadj) {

  /* Loop over all cells in the space. */
  int cid = 0;
  for (int l = 0; l < s->cdim[0]; l++) {
    for (int m = 0; m < s->cdim[1]; m++) {
      for (int n = 0; n < s->cdim[2]; n++) {

        /* Visit all neighbours of this cell, wrapping space at edges. */
        int p = 0;
        for (int i = -1; i <= 1; i++) {
          int ii = l + i;
          if (ii < 0)
            ii += s->cdim[0];
          else if (ii >= s->cdim[0])
            ii -= s->cdim[0];
          for (int j = -1; j <= 1; j++) {
            int jj = m + j;
            if (jj < 0)
              jj += s->cdim[1];
            else if (jj >= s->cdim[1])
              jj -= s->cdim[1];
            for (int k = -1; k <= 1; k++) {
              int kk = n + k;
              if (kk < 0)
                kk += s->cdim[2];
              else if (kk >= s->cdim[2])
                kk -= s->cdim[2];

              /* If not self, record id of neighbour. */
              if (i || j || k) {
                adjncy[cid * 26 + p] = cell_getid(s->cdim, ii, jj, kk);
                p++;
              }
            }
          }
        }

        /* Next cell. */
        cid++;
      }
    }
  }

  /* If given set xadj. */
  if (xadj != NULL) {
    xadj[0] = 0;
    for (int k = 0; k < s->nr_cells; k++) xadj[k + 1] = xadj[k] + 26;
  }
}
#endif

#if defined(WITH_MPI) && defined(HAVE_METIS)
/**
 * @brief Accumulate the counts of particles per cell.
 *
 * @param s the space containing the cells.
 * @param counts the number of particles per cell. Should be
 *               allocated as size s->nr_parts.
 */
static void accumulate_counts(struct space *s, int *counts) {

  struct part *parts = s->parts;
  int *cdim = s->cdim;
  double ih[3], dim[3];
  ih[0] = s->ih[0];
  ih[1] = s->ih[1];
  ih[2] = s->ih[2];
  dim[0] = s->dim[0];
  dim[1] = s->dim[1];
  dim[2] = s->dim[2];

  bzero(counts, sizeof(int) * s->nr_cells);

  for (int k = 0; k < s->nr_parts; k++) {
    for (int j = 0; j < 3; j++) {
      if (parts[k].x[j] < 0.0)
        parts[k].x[j] += dim[j];
      else if (parts[k].x[j] >= dim[j])
        parts[k].x[j] -= dim[j];
    }
    const int cid = cell_getid(cdim, parts[k].x[0] * ih[0],
                               parts[k].x[1] * ih[1], parts[k].x[2] * ih[2]);
    counts[cid]++;
  }
}
#endif

#if defined(WITH_MPI) && defined(HAVE_METIS)
/**
 * @brief Apply METIS cell list partitioning to a cell structure.
 *
 * Uses the results of part_metis_pick to assign each cell's nodeID to the
 * picked region index, thus partitioning the space into regions.
 *
 * @param s the space containing the cells to split into regions.
 * @param nregions number of regions.
 * @param celllist list of regions for each cell.
 */
static void split_metis(struct space *s, int nregions, int *celllist) {

  for (int i = 0; i < s->nr_cells; i++) s->cells[i].nodeID = celllist[i];
}
#endif

#if defined(WITH_MPI) && defined(HAVE_METIS)
/**
 * @brief Partition the given space into a number of connected regions.
 *
 * Split the space using METIS to derive a partitions using the
 * given edge and vertex weights. If no weights are given then an
 * unweighted partition is performed.
 *
 * @param s the space of cells to partition.
 * @param nregions the number of regions required in the partition.
 * @param vertexw weights for the cells, sizeof number of cells if used,
 *        NULL for unit weights.
 * @param edgew weights for the graph edges between all cells, sizeof number
 *        of cells * 26 if used, NULL for unit weights. Need to be packed
 *        in CSR format, so same as adjncy array.
 * @param celllist on exit this contains the ids of the selected regions,
 *        sizeof number of cells.
 */
static void pick_metis(struct space *s, int nregions, int *vertexw, int *edgew,
                       int *celllist) {

  /* Total number of cells. */
  int ncells = s->cdim[0] * s->cdim[1] * s->cdim[2];

  /* Nothing much to do if only using a single partition. Also avoids METIS
   * bug that doesn't handle this case well. */
  if (nregions == 1) {
    for (int i = 0; i < ncells; i++) celllist[i] = 0;
    return;
  }

  /* Allocate weights and adjacency arrays . */
  idx_t *xadj;
  if ((xadj = (idx_t *)malloc(sizeof(idx_t) * (ncells + 1))) == NULL)
    error("Failed to allocate xadj buffer.");
  idx_t *adjncy;
  if ((adjncy = (idx_t *)malloc(sizeof(idx_t) * 26 * ncells)) == NULL)
    error("Failed to allocate adjncy array.");
  idx_t *weights_v = NULL;
  if (vertexw != NULL)
    if ((weights_v = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
      error("Failed to allocate vertex weights array");
  idx_t *weights_e = NULL;
  if (edgew != NULL)
    if ((weights_e = (idx_t *)malloc(26 * sizeof(idx_t) * ncells)) == NULL)
      error("Failed to allocate edge weights array");
  idx_t *regionid;
  if ((regionid = (idx_t *)malloc(sizeof(idx_t) * ncells)) == NULL)
    error("Failed to allocate regionid array");

  /* Define the cell graph. */
  graph_init_metis(s, adjncy, xadj);

  /* Init the vertex weights array. */
  if (vertexw != NULL) {
    for (int k = 0; k < ncells; k++) {
      if (vertexw[k] > 0) {
        weights_v[k] = vertexw[k];
      } else {
        weights_v[k] = 1;
      }
    }
  }

  /* Init the edges weights array. */
  if (edgew != NULL) {
    for (int k = 0; k < 26 * ncells; k++) {
      if (edgew[k] > 0) {
        weights_e[k] = edgew[k];
      } else {
        weights_e[k] = 1;
      }
    }
  }

  /* Set the METIS options. */
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);
  options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT;
  options[METIS_OPTION_NUMBERING] = 0;
  options[METIS_OPTION_CONTIG] = 1;
  options[METIS_OPTION_NCUTS] = 10;
  options[METIS_OPTION_NITER] = 20;

  /* Call METIS. */
  idx_t one = 1;
  idx_t idx_ncells = ncells;
  idx_t idx_nregions = nregions;
  idx_t objval;

  /* Dump graph in METIS format */
  /* dumpMETISGraph("metis_graph", idx_ncells, one, xadj, adjncy,
   *                weights_v, weights_e, NULL);
   */

  if (METIS_PartGraphKway(&idx_ncells, &one, xadj, adjncy, weights_v, weights_e,
                          NULL, &idx_nregions, NULL, NULL, options, &objval,
                          regionid) != METIS_OK)
    error("Call to METIS_PartGraphKway failed.");

  /* Check that the regionids are ok. */
  for (int k = 0; k < ncells; k++)
    if (regionid[k] < 0 || regionid[k] >= nregions)
      error("Got bad nodeID %" PRIDX " for cell %i.", regionid[k], k);

  /* Set the cell list to the region index. */
  for (int k = 0; k < ncells; k++) {
    celllist[k] = regionid[k];
  }

  /* Clean up. */
  if (weights_v != NULL) free(weights_v);
  if (weights_e != NULL) free(weights_e);
  free(xadj);
  free(adjncy);
  free(regionid);
}
#endif

#if defined(WITH_MPI) && defined(HAVE_METIS)
/**
 * @brief Repartition the cells amongst the nodes using task timings
 *        as edge weights and vertex weights also from task timings
 *        or particle cells counts.
 *
 * @param partweights whether particle counts will be used as vertex weights.
 * @param bothweights whether vertex and edge weights will be used, otherwise
 *                    only edge weights will be used.
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells holding our local particles.
 * @param tasks the completed tasks from the last engine step for our node.
 * @param nr_tasks the number of tasks.
 */
static void repart_edge_metis(int partweights, int bothweights, int nodeID,
                              int nr_nodes, struct space *s, struct task *tasks,
                              int nr_tasks) {

  /* Create weight arrays using task ticks for vertices and edges (edges
   * assume the same graph structure as used in the part_ calls). */
  int nr_cells = s->nr_cells;
  struct cell *cells = s->cells;
  float wscale = 1e-3, vscale = 1e-3, wscale_buff = 0.0;
  int wtot = 0;
  int wmax = 1e9 / nr_nodes;
  int wmin;

  /* Allocate and fill the adjncy indexing array defining the graph of
   * cells. */
  idx_t *inds;
  if ((inds = (idx_t *)malloc(sizeof(idx_t) * 26 * nr_cells)) == NULL)
    error("Failed to allocate the inds array");
  graph_init_metis(s, inds, NULL);

  /* Allocate and init weights. */
  int *weights_v = NULL;
  int *weights_e = NULL;
  if (bothweights) {
    if ((weights_v = (int *)malloc(sizeof(int) * nr_cells)) == NULL)
      error("Failed to allocate vertex weights arrays.");
    bzero(weights_v, sizeof(int) * nr_cells);
  }
  if ((weights_e = (int *)malloc(sizeof(int) * 26 * nr_cells)) == NULL)
    error("Failed to allocate edge weights arrays.");
  bzero(weights_e, sizeof(int) * 26 * nr_cells);

  /* Generate task weights for vertices. */
  int taskvweights = (bothweights && !partweights);

  /* Loop over the tasks... */
  for (int j = 0; j < nr_tasks; j++) {
    /* Get a pointer to the kth task. */
    struct task *t = &tasks[j];

    /* Skip un-interesting tasks. */
    if (t->type != task_type_self && t->type != task_type_pair &&
        t->type != task_type_sub && t->type != task_type_ghost &&
        t->type != task_type_drift && t->type != task_type_kick &&
        t->type != task_type_init)
      continue;

    /* Get the task weight. */
    int w = (t->toc - t->tic) * wscale;
    if (w < 0) error("Bad task weight (%d).", w);

    /* Do we need to re-scale? */
    wtot += w;
    while (wtot > wmax) {
      wscale /= 2;
      wtot /= 2;
      w /= 2;
      for (int k = 0; k < 26 * nr_cells; k++) weights_e[k] *= 0.5;
      if (taskvweights)
        for (int k = 0; k < nr_cells; k++) weights_v[k] *= 0.5;
    }

    /* Get the top-level cells involved. */
    struct cell *ci, *cj;
    for (ci = t->ci; ci->parent != NULL; ci = ci->parent)
      ;
    if (t->cj != NULL)
      for (cj = t->cj; cj->parent != NULL; cj = cj->parent)
        ;
    else
      cj = NULL;

    /* Get the cell IDs. */
    int cid = ci - cells;

    /* Different weights for different tasks. */
    if (t->type == task_type_ghost || t->type == task_type_drift ||
        t->type == task_type_kick) {
      /* Particle updates add only to vertex weight. */
      if (taskvweights) weights_v[cid] += w;

    }

    /* Self interaction? */
    else if ((t->type == task_type_self && ci->nodeID == nodeID) ||
             (t->type == task_type_sub && cj == NULL && ci->nodeID == nodeID)) {
      /* Self interactions add only to vertex weight. */
      if (taskvweights) weights_v[cid] += w;

    }

    /* Pair? */
    else if (t->type == task_type_pair ||
             (t->type == task_type_sub && cj != NULL)) {
      /* In-cell pair? */
      if (ci == cj) {
        /* Add weight to vertex for ci. */
        if (taskvweights) weights_v[cid] += w;

      }

      /* Distinct cells with local ci? */
      else if (ci->nodeID == nodeID) {
        /* Index of the jth cell. */
        int cjd = cj - cells;

        /* Add half of weight to each cell. */
        if (taskvweights) {
          if (ci->nodeID == nodeID) weights_v[cid] += 0.5 * w;
          if (cj->nodeID == nodeID) weights_v[cjd] += 0.5 * w;
        }

        /* Add weights to edge. */
        int kk;
        for (kk = 26 * cid; inds[kk] != cjd; kk++)
          ;
        weights_e[kk] += w;
        for (kk = 26 * cjd; inds[kk] != cid; kk++)
          ;
        weights_e[kk] += w;
      }
    }
  }

  /* Re-calculate the vertices if using particle counts. */
  if (partweights && bothweights) {
    accumulate_counts(s, weights_v);

    /*  Rescale to balance times. */
    float vwscale = (float)wtot / (float)nr_tasks;
    for (int k = 0; k < nr_cells; k++) {
      weights_v[k] *= vwscale;
    }
  }

  /* Get the minimum scaling and re-scale if necessary. */
  int res;
  if ((res = MPI_Allreduce(&wscale, &wscale_buff, 1, MPI_FLOAT, MPI_MIN,
                           MPI_COMM_WORLD)) != MPI_SUCCESS)
    mpi_error(res, "Failed to allreduce the weight scales.");

  if (wscale_buff != wscale) {
    float scale = wscale_buff / wscale;
    for (int k = 0; k < 26 * nr_cells; k++) weights_e[k] *= scale;
    if (bothweights)
      for (int k = 0; k < nr_cells; k++) weights_v[k] *= scale;
  }

  /* Merge the weights arrays across all nodes. */
  if (bothweights) {
    if ((res = MPI_Reduce((nodeID == 0) ? MPI_IN_PLACE : weights_v, weights_v,
                          nr_cells, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)) !=
        MPI_SUCCESS)
      mpi_error(res, "Failed to allreduce vertex weights.");
  }

  if ((res = MPI_Reduce((nodeID == 0) ? MPI_IN_PLACE : weights_e, weights_e,
                        26 * nr_cells, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD)) !=
      MPI_SUCCESS)
    mpi_error(res, "Failed to allreduce edge weights.");

  /* Allocate cell list for the partition. */
  int *celllist = (int *)malloc(sizeof(int) * s->nr_cells);
  if (celllist == NULL) error("Failed to allocate celllist");

  /* As of here, only one node needs to compute the partition. */
  if (nodeID == 0) {
    /* Final rescale of all weights to avoid a large range. Large ranges
     * have been seen to cause an incomplete graph. */
    wmin = wmax;
    wmax = 0;
    for (int k = 0; k < 26 * nr_cells; k++) {
      wmax = weights_e[k] > wmax ? weights_e[k] : wmax;
      wmin = weights_e[k] < wmin ? weights_e[k] : wmin;
    }
    if (bothweights) {
      for (int k = 0; k < nr_cells; k++) {
        wmax = weights_v[k] > wmax ? weights_v[k] : wmax;
        wmin = weights_v[k] < wmin ? weights_v[k] : wmin;
      }
    }

    if ((wmax - wmin) > metis_maxweight) {
      wscale = metis_maxweight / (wmax - wmin);
      for (int k = 0; k < 26 * nr_cells; k++) {
        weights_e[k] = (weights_e[k] - wmin) * wscale + 1;
      }
      if (bothweights) {
        for (int k = 0; k < nr_cells; k++) {
          weights_v[k] = (weights_v[k] - wmin) * wscale + 1;
        }
      }
    }

    /* Make sure there are no zero weights. */
    for (int k = 0; k < 26 * nr_cells; k++)
      if (weights_e[k] == 0) weights_e[k] = 1;
    if (bothweights)
      for (int k = 0; k < nr_cells; k++)
        if ((weights_v[k] *= vscale) == 0) weights_v[k] = 1;

    /* And partition, use both weights or not as requested. */
    if (bothweights)
      pick_metis(s, nr_nodes, weights_v, weights_e, celllist);
    else
      pick_metis(s, nr_nodes, NULL, weights_e, celllist);

    /* Check that all cells have good values. */
    for (int k = 0; k < nr_cells; k++)
      if (celllist[k] < 0 || celllist[k] >= nr_nodes)
        error("Got bad nodeID %d for cell %i.", celllist[k], k);

    /* Check that the partition is complete and all nodes have some work. */
    int present[nr_nodes];
    int failed = 0;
    for (int i = 0; i < nr_nodes; i++) present[i] = 0;
    for (int i = 0; i < nr_cells; i++) present[celllist[i]]++;
    for (int i = 0; i < nr_nodes; i++) {
      if (!present[i]) {
        failed = 1;
        message("Node %d is not present after repartition", i);
      }
    }

    /* If partition failed continue with the current one, but make this
     * clear. */
    if (failed) {
      message(
          "WARNING: METIS repartition has failed, continuing with "
          "the current partition, load balance will not be optimal");
      for (int k = 0; k < nr_cells; k++) celllist[k] = cells[k].nodeID;
    }
  }

  /* Distribute the celllist partition and apply. */
  if ((res = MPI_Bcast(celllist, s->nr_cells, MPI_INT, 0, MPI_COMM_WORLD)) !=
      MPI_SUCCESS)
    mpi_error(res, "Failed to bcast the cell list");

  /* And apply to our cells */
  split_metis(s, nr_nodes, celllist);

  /* Clean up. */
  if (bothweights) free(weights_v);
  free(weights_e);
  free(celllist);
}
#endif

/**
 * @brief Repartition the cells amongst the nodes using vertex weights
 *
 * @param s The space containing the local cells.
 * @param nodeID our MPI node id.
 * @param nr_nodes number of MPI nodes.
 */
#if defined(WITH_MPI) && defined(HAVE_METIS)
static void repart_vertex_metis(struct space *s, int nodeID, int nr_nodes) {

  /* Use particle counts as vertex weights. */
  /* Space for particles per cell counts, which will be used as weights. */
  int *weights = NULL;
  if ((weights = (int *)malloc(sizeof(int) * s->nr_cells)) == NULL)
    error("Failed to allocate weights buffer.");

  /* Check each particle and accumulate the counts per cell. */
  accumulate_counts(s, weights);

  /* Get all the counts from all the nodes. */
  int res;
  if ((res = MPI_Allreduce(MPI_IN_PLACE, weights, s->nr_cells, MPI_INT, MPI_SUM,
                           MPI_COMM_WORLD)) != MPI_SUCCESS)
    mpi_error(res, "Failed to allreduce particle cell weights.");

  /* Main node does the partition calculation. */
  int *celllist = (int *)malloc(sizeof(int) * s->nr_cells);
  if (celllist == NULL) error("Failed to allocate celllist");

  if (nodeID == 0) pick_metis(s, nr_nodes, weights, NULL, celllist);

  /* Distribute the celllist partition and apply. */
  if ((res = MPI_Bcast(celllist, s->nr_cells, MPI_INT, 0, MPI_COMM_WORLD)) !=
      MPI_SUCCESS)
    mpi_error(res, "Failed to bcast the cell list");

  /* And apply to our cells */
  split_metis(s, nr_nodes, celllist);

  free(weights);
  free(celllist);
}
#endif

/**
 * @brief Repartition the space using the given repartition type.
 *
 * Note that at the end of this process all the cells will be re-distributed
 * across the nodes, but the particles themselves will not be.
 *
 * @param reparttype the type of repartition to attempt, see the repart_type
 *enum.
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells holding our local particles.
 * @param tasks the completed tasks from the last engine step for our node.
 * @param nr_tasks the number of tasks.
 */
void partition_repartition(enum repartition_type reparttype, int nodeID,
                           int nr_nodes, struct space *s, struct task *tasks,
                           int nr_tasks) {

#if defined(WITH_MPI) && defined(HAVE_METIS)

  if (reparttype == REPART_METIS_BOTH || reparttype == REPART_METIS_EDGE ||
      reparttype == REPART_METIS_VERTEX_EDGE) {

    int partweights;
    int bothweights;
    if (reparttype == REPART_METIS_VERTEX_EDGE)
      partweights = 1;
    else
      partweights = 0;

    if (reparttype == REPART_METIS_BOTH)
      bothweights = 1;
    else
      bothweights = 0;

    repart_edge_metis(partweights, bothweights, nodeID, nr_nodes, s, tasks,
                      nr_tasks);

  } else if (reparttype == REPART_METIS_VERTEX) {

    repart_vertex_metis(s, nodeID, nr_nodes);

  } else {
    error("Unknown repartition type");
  }
#else
  error("SWIFT was not compiled with METIS support.");
#endif
}

/**
 * @brief Initial partition of space cells.
 *
 * Cells are assigned to a node on the basis of various schemes, all of which
 * should attempt to distribute them in geometrically close regions to
 * minimise the movement of particles.
 *
 * Note that the partition type is a suggestion and will be ignored if that
 * scheme fails. In that case we fallback to a vectorised scheme, that is
 * guaranteed to work provided we have more cells than nodes.
 *
 * @param initial_partition the type of partitioning to try.
 * @param nodeID our nodeID.
 * @param nr_nodes the number of nodes.
 * @param s the space of cells.
 */
void partition_initial_partition(struct partition *initial_partition,
                                 int nodeID, int nr_nodes, struct space *s) {

  /* Geometric grid partitioning. */
  if (initial_partition->type == INITPART_GRID) {
    int j, k;
    int ind[3];
    struct cell *c;

    /* If we've got the wrong number of nodes, fail. */
    if (nr_nodes != initial_partition->grid[0] * initial_partition->grid[1] *
                        initial_partition->grid[2])
      error("Grid size does not match number of nodes.");

    /* Run through the cells and set their nodeID. */
    // message("s->dim = [%e,%e,%e]", s->dim[0], s->dim[1], s->dim[2]);
    for (k = 0; k < s->nr_cells; k++) {
      c = &s->cells[k];
      for (j = 0; j < 3; j++)
        ind[j] = c->loc[j] / s->dim[j] * initial_partition->grid[j];
      c->nodeID = ind[0] + initial_partition->grid[0] *
                               (ind[1] + initial_partition->grid[1] * ind[2]);
      // message("cell at [%e,%e,%e]: ind = [%i,%i,%i], nodeID = %i", c->loc[0],
      // c->loc[1], c->loc[2], ind[0], ind[1], ind[2], c->nodeID);
    }

    /* The grid technique can fail, so check for this before proceeding. */
    if (!check_complete(s, (nodeID == 0), nr_nodes)) {
      if (nodeID == 0)
        message("Grid initial partition failed, using a vectorised partition");
      initial_partition->type = INITPART_VECTORIZE;
      partition_initial_partition(initial_partition, nodeID, nr_nodes, s);
      return;
    }

  } else if (initial_partition->type == INITPART_METIS_WEIGHT ||
             initial_partition->type == INITPART_METIS_NOWEIGHT) {
#if defined(WITH_MPI) && defined(HAVE_METIS)
    /* Simple k-way partition selected by METIS using cell particle counts as
     * weights or not. Should be best when starting with a inhomogeneous dist.
     */

    /* Space for particles per cell counts, which will be used as weights or
     * not. */
    int *weights = NULL;
    if (initial_partition->type == INITPART_METIS_WEIGHT) {
      if ((weights = (int *)malloc(sizeof(int) * s->nr_cells)) == NULL)
        error("Failed to allocate weights buffer.");
      bzero(weights, sizeof(int) * s->nr_cells);

      /* Check each particle and accumilate the counts per cell. */
      struct part *parts = s->parts;
      int *cdim = s->cdim;
      double ih[3], dim[3];
      ih[0] = s->ih[0];
      ih[1] = s->ih[1];
      ih[2] = s->ih[2];
      dim[0] = s->dim[0];
      dim[1] = s->dim[1];
      dim[2] = s->dim[2];
      for (int k = 0; k < s->nr_parts; k++) {
        for (int j = 0; j < 3; j++) {
          if (parts[k].x[j] < 0.0)
            parts[k].x[j] += dim[j];
          else if (parts[k].x[j] >= dim[j])
            parts[k].x[j] -= dim[j];
        }
        const int cid =
            cell_getid(cdim, parts[k].x[0] * ih[0], parts[k].x[1] * ih[1],
                       parts[k].x[2] * ih[2]);
        weights[cid]++;
      }

      /* Get all the counts from all the nodes. */
      if (MPI_Allreduce(MPI_IN_PLACE, weights, s->nr_cells, MPI_INT, MPI_SUM,
                        MPI_COMM_WORLD) != MPI_SUCCESS)
        error("Failed to allreduce particle cell weights.");
    }

    /* Main node does the partition calculation. */
    int *celllist = (int *)malloc(sizeof(int) * s->nr_cells);
    if (celllist == NULL) error("Failed to allocate celllist");
    if (nodeID == 0) pick_metis(s, nr_nodes, weights, NULL, celllist);

    /* Distribute the celllist partition and apply. */
    int res = MPI_Bcast(celllist, s->nr_cells, MPI_INT, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS) mpi_error(res, "Failed to bcast the cell list");

    /* And apply to our cells */
    split_metis(s, nr_nodes, celllist);

    /* It's not known if this can fail, but check for this before
     * proceeding. */
    if (!check_complete(s, (nodeID == 0), nr_nodes)) {
      if (nodeID == 0)
        message("METIS initial partition failed, using a vectorised partition");
      initial_partition->type = INITPART_VECTORIZE;
      partition_initial_partition(initial_partition, nodeID, nr_nodes, s);
    }

    if (weights != NULL) free(weights);
    free(celllist);
#else
    error("SWIFT was not compiled with METIS support");
#endif

  } else if (initial_partition->type == INITPART_VECTORIZE) {

#if defined(WITH_MPI)
    /* Vectorised selection, guaranteed to work for samples less than the
     * number of cells, but not very clumpy in the selection of regions. */
    int *samplecells = (int *)malloc(sizeof(int) * nr_nodes * 3);
    if (samplecells == NULL) error("Failed to allocate samplecells");

    if (nodeID == 0) {
      pick_vector(s, nr_nodes, samplecells);
    }

    /* Share the samplecells around all the nodes. */
    int res = MPI_Bcast(samplecells, nr_nodes * 3, MPI_INT, 0, MPI_COMM_WORLD);
    if (res != MPI_SUCCESS)
      mpi_error(res, "Failed to bcast the partition sample cells.");

    /* And apply to our cells */
    split_vector(s, nr_nodes, samplecells);
    free(samplecells);
#else
    error("SWIFT was not compiled with MPI support");
#endif
  }
}

/*  General support */
/*  =============== */

/**
 * @brief Check if all regions have been assigned a node in the
 *        cells of a space.
 *
 * @param s the space containing the cells to check.
 * @param nregions number of regions expected.
 * @param verbose if true report the missing regions.
 * @return true if all regions have been found, false otherwise.
 */
static int check_complete(struct space *s, int verbose, int nregions) {

  int *present = (int *)malloc(sizeof(int) * nregions);
  if (present == NULL) error("Failed to allocate present array");

  int failed = 0;
  for (int i = 0; i < nregions; i++) present[i] = 0;
  for (int i = 0; i < s->nr_cells; i++) {
    present[s->cells[i].nodeID]++;
  }
  for (int i = 0; i < nregions; i++) {
    if (!present[i]) {
      failed = 1;
      if (verbose) message("Region %d is not present in partition", i);
    }
  }
  free(present);
  return (!failed);
}
