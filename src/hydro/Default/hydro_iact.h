/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_RUNNER_IACT_H
#define SWIFT_RUNNER_IACT_H

/* Includes. */
#include "const.h"
#include "kernel.h"
#include "part.h"
#include "vector.h"

/**
 * @brief SPH interaction functions following the Gadget-2 version of SPH.
 *
 * The interactions computed here are the ones presented in the Gadget-2 paper
 *and use the same
 * numerical coefficients as the Gadget-2 code. When used with the Spline-3
 *kernel, the results
 * should be equivalent to the ones obtained with Gadget-2 up to the rounding
 *errors and interactions
 * missed by the Gadget-2 tree-code neighbours search.
 *
 * The code uses internal energy instead of entropy as a thermodynamical
 *variable.
 */

/**
 * @brief Density loop
 */

__attribute__((always_inline)) INLINE static void runner_iact_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float r = sqrtf(r2), ri = 1.0f / r;                                           // 1xSQRT, 1xRECIP                     2
  float xi, xj;
  float h_inv;
  float wi, wj, wi_dx, wj_dx;
  float mi, mj;
  float dvdr;
  float dv[3], curlvr[3];
  int k;

  /* Get the masses. */
  mi = pi->mass;
  mj = pj->mass;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];                                                  // 1xSUB                               1
  dv[1] = pi->v[1] - pj->v[1];                                                  // 1xSUB                               1
  dv[2] = pi->v[2] - pj->v[2];                                                  // 1xSUB                               1
  dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];                         // 3xMUL, 2xADD                        5
  dvdr *= ri;                                                                   // 1xMUL                               1

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];                                    // 2xMUL, 1xSUB                        3
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];                                    // 2xMUL, 1xSUB                        3
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];                                    // 2xMUL, 1xSUB                        3
  for (k = 0; k < 3; k++) curlvr[k] *= ri;                                      // 3xMUL                               3

  /* Compute density of pi. */
  h_inv = 1.0 / hi;                                                             // 1xRECIP                             1
  xi = r * h_inv;                                                               // 1xMUL                               1
  kernel_deval(xi, &wi, &wi_dx);                                                // 11xFLOP                            11

  pi->rho += mj * wi;                                                           // 1xADD, 1xMUL                        2
  pi->rho_dh -= mj * (3.0 * wi + xi * wi_dx);                                   // 1xSUB, 3xMUL, 1xADD                 5
  pi->density.wcount += wi;                                                     // 1xADD                               1
  pi->density.wcount_dh -= xi * wi_dx;                                          // 1xSUB, 1xMUL                        2

  pi->density.div_v += mj * dvdr * wi_dx;                                       // 1xADD, 2xMUL                        3
  for (k = 0; k < 3; k++) pi->density.curl_v[k] += mj * curlvr[k] * wi_dx;      // 1xADD, 2xMUL                        3

  /* Compute density of pj. */
  h_inv = 1.0 / hj;                                                             // 1xRECIP                             1
  xj = r * h_inv;                                                               // 1xMUL                               1
  kernel_deval(xj, &wj, &wj_dx);                                                // 11xFLOP                            11

  pj->rho += mi * wj;                                                           // 1xADD, 1xMUL                        2
  pj->rho_dh -= mi * (3.0 * wj + xj * wj_dx);                                   // 1xSUB, 3xMUL, 1xADD                 5
  pj->density.wcount += wj;                                                     // 1xADD                               1
  pj->density.wcount_dh -= xj * wj_dx;                                          // 1xSUB, 1xMUL                        2

  pj->density.div_v += mi * dvdr * wj_dx;                                       // 1xADD, 2xMUL                        3
  for (k = 0; k < 3; k++) pj->density.curl_v[k] += mi * curlvr[k] * wj_dx;      // 3xADD, 3x2xMUL                      9
                                                                                //
                                                                                // 87 FLOPs
}

/**
 * @brief Density loop (Vectorized version)
 */
__attribute__((always_inline)) INLINE static void runner_iact_vec_density(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {

#ifdef VECTORIZE

  vector r, ri, r2, xi, xj, hi, hj, hi_inv, hj_inv, wi, wj, wi_dx, wj_dx;
  vector rhoi, rhoj, rhoi_dh, rhoj_dh, wcounti, wcountj, wcounti_dh, wcountj_dh;
  vector mi, mj;
  vector dx[3], dv[3];
  vector vi[3], vj[3];
  vector dvdr, div_vi, div_vj;
  vector curlvr[3], curl_vi[3], curl_vj[3];
  int k, j;

#if VEC_SIZE == 8
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass,
                 pi[4]->mass, pi[5]->mass, pi[6]->mass, pi[7]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass,
                 pj[4]->mass, pj[5]->mass, pj[6]->mass, pj[7]->mass);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k],
                      pi[4]->v[k], pi[5]->v[k], pi[6]->v[k], pi[7]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k],
                      pj[4]->v[k], pj[5]->v[k], pj[6]->v[k], pj[7]->v[k]);
  }
#elif VEC_SIZE == 4
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k]);
  }
#endif
  for (k = 0; k < 3; k++)
    dx[k].v = vec_load(Dx + VEC_SIZE * k);

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri.v = vec_rsqrt(r2.v);                                                       // 8xSQRT, 8xRECIP                    16
  ri.v = ri.v - vec_set1(0.5f) * ri.v * (r2.v * ri.v * ri.v - vec_set1(1.0f));  // 2x8xSUB, 4x8xMUL                   48
  r.v = r2.v * ri.v;                                                            // 8xMUL                               8

  hi.v = vec_load(Hi);
  hi_inv.v = vec_rcp(hi.v);                                                     // 8xRECIP                             8
  hi_inv.v = hi_inv.v - hi_inv.v * (hi_inv.v * hi.v - vec_set1(1.0f));          // 2x8xSUB, 2x8xMUL                   32
  xi.v = r.v * hi_inv.v;                                                        // 8xMUL                               8

  hj.v = vec_load(Hj);
  hj_inv.v = vec_rcp(hj.v);                                                     // 8xRECIP                             8
  hj_inv.v = hj_inv.v - hj_inv.v * (hj_inv.v * hj.v - vec_set1(1.0f));          // 2x8xSUB + 2x8xMUL                  32
  xj.v = r.v * hj_inv.v;                                                        // 8xMUL                               8

  kernel_deval_vec(&xi, &wi, &wi_dx);                                           // 88xFLOP                            88
  kernel_deval_vec(&xj, &wj, &wj_dx);                                           // 88xFLOP                            88

  /* Compute dv. */
  dv[0].v = vi[0].v - vj[0].v;                                                  // 8xSUB                               8
  dv[1].v = vi[1].v - vj[1].v;                                                  // 8xSUB                               8
  dv[2].v = vi[2].v - vj[2].v;                                                  // 8xSUB                               8

  /* Compute dv dot r */
  dvdr.v = (dv[0].v * dx[0].v) + (dv[1].v * dx[1].v) + (dv[2].v * dx[2].v);     // 3x8xMUL, 2x8xADD                   40
  dvdr.v = dvdr.v * ri.v;                                                       // 8xMUL                               8

  /* Compute dv cross r */
  curlvr[0].v = dv[1].v * dx[2].v - dv[2].v * dx[1].v;                          // 2x8xMUL, 8xSUB                     24
  curlvr[1].v = dv[2].v * dx[0].v - dv[0].v * dx[2].v;                          // 2x8xMUL, 8xSUB                     24
  curlvr[2].v = dv[0].v * dx[1].v - dv[1].v * dx[0].v;                          // 2x8xMUL, 8xSUB                     24
  for (k = 0; k < 3; k++) curlvr[k].v *= ri.v;                                  // 3x8xMUL                            24

  rhoi.v = mj.v * wi.v;                                                         // 8xMUL                               8
  rhoi_dh.v = mj.v * (vec_set1(3.0f) * wi.v + xi.v * wi_dx.v);                  // 3x8xMUL, 8xADD                     32
  wcounti.v = wi.v;
  wcounti_dh.v = xi.v * wi_dx.v;                                                // 8xMUL                               8
  div_vi.v = mj.v * dvdr.v * wi_dx.v;                                           // 2x8xMUL                            16
  for (k = 0; k < 3; k++) curl_vi[k].v = mj.v * curlvr[k].v * wi_dx.v;          // 3x2x8xMUL                          48

  rhoj.v = mi.v * wj.v;                                                         // 8xMUL                               8
  rhoj_dh.v = mi.v * (vec_set1(3.0f) * wj.v + xj.v * wj_dx.v);                  // 3x8xMUL, 8xADD                     32
  wcountj.v = wj.v;
  wcountj_dh.v = xj.v * wj_dx.v;                                                // 8xMUL                               8
  div_vj.v = mi.v * dvdr.v * wj_dx.v;                                           // 2x8xMUL                            16
  for (k = 0; k < 3; k++) curl_vj[k].v = mi.v * curlvr[k].v * wj_dx.v;          // 3x2x8xMUL                          48

  for (k = 0; k < VEC_SIZE; k++) {                                              // 8x...                             128
    pi[k]->rho += rhoi.f[k];                                                    //   1xADD
    pi[k]->rho_dh -= rhoi_dh.f[k];                                              //   1xSUB
    pi[k]->density.wcount += wcounti.f[k];                                      //   1xADD
    pi[k]->density.wcount_dh -= wcounti_dh.f[k];                                //   1xSUB
    pi[k]->density.div_v += div_vi.f[k];                                        //   1xADD
    for (j = 0; j < 3; j++) pi[k]->density.curl_v[j] += curl_vi[j].f[k];        //   3xADD
    pj[k]->rho += rhoj.f[k];                                                    //   1xADD
    pj[k]->rho_dh -= rhoj_dh.f[k];                                              //   1xSUB
    pj[k]->density.wcount += wcountj.f[k];                                      //   1xADD
    pj[k]->density.wcount_dh -= wcountj_dh.f[k];                                //   1xSUB
    pj[k]->density.div_v += div_vj.f[k];                                        //   1xADD
    for (j = 0; j < 3; j++) pj[k]->density.curl_v[j] += curl_vj[j].f[k];        //   3xADD
  }
                                                                                //
                                                                                // 776 FLOPs

#else

  for (int k = 0; k < VEC_SIZE; k++)
    runner_iact_density(R2[k], &Dx[3 * k], Hi[k], Hj[k], pi[k], pj[k]);

#endif
}

/**
 * @brief Density loop (non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float r, ri;
  float xi;
  float h_inv;
  float wi, wi_dx;
  float mj;
  float dvdr;
  float dv[3], curlvr[3];
  int k;

  /* Get the masses. */
  mj = pj->mass;

  /* Get r and r inverse. */
  r = sqrtf(r2);                                                                // 1xSQRT                              1
  ri = 1.0f / r;                                                                // 1xRECIP                             1

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];                                                  // 1xSUB                               1
  dv[1] = pi->v[1] - pj->v[1];                                                  // 1xSUB                               1
  dv[2] = pi->v[2] - pj->v[2];                                                  // 1xSUB                               1
  dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];                         // 3xMUL, 2xADD                        5
  dvdr *= ri;                                                                   // 1xMUL                               1

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];                                    // 2xMUL, 1xSUB                        3
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];                                    // 2xMUL, 1xSUB                        3
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];                                    // 2xMUL, 1xSUB                        3
  for (k = 0; k < 3; k++) curlvr[k] *= ri;                                      // 3xMUL                               3

  h_inv = 1.0 / hi;                                                             // 1xRECIP                             1
  xi = r * h_inv;                                                               // 1xMUL                               1
  kernel_deval(xi, &wi, &wi_dx);                                                // 11xFLOP                            11

  pi->rho += mj * wi;                                                           // 1xADD, 1xMUL                        2
  pi->rho_dh -= mj * (3.0 * wi + xi * wi_dx);                                   // 1xSUB, 3xMUL, 1xADD                 5
  pi->density.wcount += wi;                                                     // 1xADD                               1
  pi->density.wcount_dh -= xi * wi_dx;                                          // 1xSUB, 1xMUL                        2

  pi->density.div_v += mj * dvdr * wi_dx;                                       // 1xADD, 2xMUL                        3
  for (k = 0; k < 3; k++) pi->density.curl_v[k] += mj * curlvr[k] * wi_dx;      // 3xADD, 3x2xMUL                      9
                                                                                //
                                                                                // 58 FLOPs
}

/**
 * @brief Density loop (non-symmetric vectorized version)
 */

__attribute__((always_inline))
    INLINE static void runner_iact_nonsym_vec_density(float *R2, float *Dx,
                                                      float *Hi, float *Hj,
                                                      struct part **pi,
                                                      struct part **pj) {

#ifdef VECTORIZE

  vector r, ri, r2, xi, hi, hi_inv, wi, wi_dx;
  vector rhoi, rhoi_dh, wcounti, wcounti_dh, div_vi;
  vector mj;
  vector dx[3], dv[3];
  vector vi[3], vj[3];
  vector dvdr;
  vector curlvr[3], curl_vi[3];
  int k, j;

#if VEC_SIZE == 8
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass,
                 pj[4]->mass, pj[5]->mass, pj[6]->mass, pj[7]->mass);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k],
                      pi[4]->v[k], pi[5]->v[k], pi[6]->v[k], pi[7]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k],
                      pj[4]->v[k], pj[5]->v[k], pj[6]->v[k], pj[7]->v[k]);
  }
#elif VEC_SIZE == 4
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k]);
  }
#endif
  for (k = 0; k < 3; k++)
    dx[k].v = vec_load(Dx + VEC_SIZE * k);

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri.v = vec_rsqrt(r2.v);                                                       // 8xSQRT, 8xRECIP                    16
  ri.v = ri.v - vec_set1(0.5f) * ri.v * (r2.v * ri.v * ri.v - vec_set1(1.0f));  // 2x8xSUB, 4x8xMUL                   48
  r.v = r2.v * ri.v;                                                            // 8xMUL                               8

  hi.v = vec_load(Hi);
  hi_inv.v = vec_rcp(hi.v);                                                     // 8xRECIP                             8
  hi_inv.v = hi_inv.v - hi_inv.v * (hi_inv.v * hi.v - vec_set1(1.0f));          // 2x8xSUB, 2x8xMUL                   32
  xi.v = r.v * hi_inv.v;                                                        // 8xMUL                               8

  kernel_deval_vec(&xi, &wi, &wi_dx);                                           // 88xFLOP                            88

  /* Compute dv. */
  dv[0].v = vi[0].v - vj[0].v;                                                  // 8xSUB                               8
  dv[1].v = vi[1].v - vj[1].v;                                                  // 8xSUB                               8
  dv[2].v = vi[2].v - vj[2].v;                                                  // 8xSUB                               8

  /* Compute dv dot r */
  dvdr.v = (dv[0].v * dx[0].v) + (dv[1].v * dx[1].v) + (dv[2].v * dx[2].v);     // 3x8xMUL, 2x8xADD                   40
  dvdr.v = dvdr.v * ri.v;                                                       // 8xMUL                               8

  /* Compute dv cross r */
  curlvr[0].v = dv[1].v * dx[2].v - dv[2].v * dx[1].v;                          // 2x8xMUL, 8xSUB                     24
  curlvr[1].v = dv[2].v * dx[0].v - dv[0].v * dx[2].v;                          // 2x8xMUL, 8xSUB                     24
  curlvr[2].v = dv[0].v * dx[1].v - dv[1].v * dx[0].v;                          // 2x8xMUL, 8xSUB                     24
  for (k = 0; k < 3; k++) curlvr[k].v *= ri.v;                                  // 3x8xMUL                            24

  rhoi.v = mj.v * wi.v;                                                         // 8xMUL                               8
  rhoi_dh.v = mj.v * (vec_set1(3.0f) * wi.v + xi.v * wi_dx.v);                  // 3x8xMUL, 8xADD                     32
  wcounti.v = wi.v;
  wcounti_dh.v = xi.v * wi_dx.v;                                                // 8xMUL                               8
  div_vi.v = mj.v * dvdr.v * wi_dx.v;                                           // 2x8xMUL                            16
  for (k = 0; k < 3; k++) curl_vi[k].v = mj.v * curlvr[k].v * wi_dx.v;          // 3x2x8xMUL                          48

  for (k = 0; k < VEC_SIZE; k++) {                                              // 8x...                             232
    pi[k]->rho += rhoi.f[k];                                                    //   1xADD
    pi[k]->rho_dh -= rhoi_dh.f[k];                                              //   1xSUB
    pi[k]->density.wcount += wcounti.f[k];                                      //   1xADD
    pi[k]->density.wcount_dh -= wcounti_dh.f[k];                                //   1xSUB
    pi[k]->density.div_v += div_vi.f[k];                                        //   1xADD
    for (j = 0; j < 3; j++) pi[k]->density.curl_v[j] += curl_vi[j].f[k];        //   3x8xADD
  }
                                                                                //
                                                                                // 720 FLOPs

#else

  for (int k = 0; k < VEC_SIZE; k++)
    runner_iact_nonsym_density(R2[k], &Dx[3 * k], Hi[k], Hj[k], pi[k], pj[k]);

#endif
}

/**
 * @brief Force loop
 */

__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float r = sqrtf(r2), ri = 1.0f / r;
  float xi, xj;
  float hi_inv, hi2_inv;
  float hj_inv, hj2_inv;
  float wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, w, dvdr;
  float mi, mj, POrho2i, POrho2j, rhoi, rhoj;
  float v_sig, omega_ij, Pi_ij, alpha_ij, tc, v_sig_u;
  float f;
  int k;

  /* Get some values in local variables. */
  mi = pi->mass;
  mj = pj->mass;
  rhoi = pi->rho;
  rhoj = pj->rho;
  POrho2i = pi->force.POrho2;
  POrho2j = pj->force.POrho2;

  /* Get the kernel for hi. */
  hi_inv = 1.0f / hi;                                                           // 1xRECIP                             1
  hi2_inv = hi_inv * hi_inv;                                                    // 1xMUL                               1
  xi = r * hi_inv;                                                              // 1xMUL                               1
  kernel_deval(xi, &wi, &wi_dx);                                                // 11xFLOP                            11
  wi_dr = hi2_inv * hi2_inv * wi_dx;                                            // 2xMUL                               2

  /* Get the kernel for hj. */
  hj_inv = 1.0f / hj;                                                           // 1xRECIP                             1
  hj2_inv = hj_inv * hj_inv;                                                    // 1xMUL                               1
  xj = r * hj_inv;                                                              // 1xMUL                               1
  kernel_deval(xj, &wj, &wj_dx);                                                // 11xFLOP                            11
  wj_dr = hj2_inv * hj2_inv * wj_dx;                                            // 2xMUL                               2

  /* Compute dv dot r. */
  dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] +        // 2xSUB, 2xMUL, 2xADD                 6
         (pi->v[2] - pj->v[2]) * dx[2];                                         // 1xSUB, 1xMUL                        2
  dvdr *= ri;                                                                   // 1xMUL                               1

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij = fminf(dvdr, 0.f);                                                  // 1xFMIN                              1

  /* Compute signal velocity */
  v_sig = pi->force.c + pj->force.c - 2.0f * omega_ij;                          // 1xADD, 1xSUB, 1xMUL                 3

  /* Compute viscosity parameter */
  alpha_ij = -0.5f * (pi->alpha + pj->alpha);                                   // 1xMUL, 1xADD                        2

  /* Compute viscosity tensor */
  Pi_ij = alpha_ij * v_sig * omega_ij / (rhoi + rhoj);                          // 12xMUL, 1xDIV, 1xADD               14

  /* Apply balsara switch */
  Pi_ij *= (pi->force.balsara + pj->force.balsara);                             // 1xMUL, 1xADD                        2

  /* Thermal conductivity */
  v_sig_u = sqrtf(2.f * (const_hydro_gamma - 1.f) *                             // 1xSQRT, 2xMUL, 1xSUB                4
                  fabs(rhoi * pi->u - rhoj * pj->u) / (rhoi + rhoj));           // 1xFABS, 2xMUL, 1xSUB, 1xDIV, 1xADD  6
  tc = const_conductivity_alpha * v_sig_u / (rhoi + rhoj);                      // 1xMUL, 1xDIV, 1xADD                 3
  tc *= (wi_dr + wj_dr);                                                        // 1xMUL, 1xADD                        2

  /* Get the common factor out. */
  w = ri *                                                                      // 1xMUL                               1
      ((POrho2i * wi_dr + POrho2j * wj_dr) + 0.25f * Pi_ij * (wi_dr + wj_dr));  // 4xMUL, 3xADD                        7

  /* Use the force, Luke! */
  for (k = 0; k < 3; k++) {                                                     // 3x...                              15
    f = dx[k] * w;                                                              //   1xMUL
    pi->a_hydro[k] -= mj * f;                                                   //   1xSUB, 1xMUL
    pj->a_hydro[k] += mi * f;                                                   //   1xADD, 1xMUL
  }

  /* Get the time derivative for u. */
  pi->force.u_dt +=                                                             // 1xADD                               1
      mj * dvdr * (POrho2i * wi_dr + 0.125f * Pi_ij * (wi_dr + wj_dr));         // 5xMUL, 2xADD                        7
  pj->force.u_dt +=                                                             // 1xADD                               1
      mi * dvdr * (POrho2j * wj_dr + 0.125f * Pi_ij * (wi_dr + wj_dr));         // 5xMUL, 2xADD                        7

  /* Add the thermal conductivity */
  pi->force.u_dt += mj * tc * (pi->u - pj->u);                                  // 1xADD, 2xMUL, 1xSUB                 4
  pj->force.u_dt += mi * tc * (pj->u - pi->u);                                  // 1xADD, 2xMUL, 1xSUB                 4

  /* Get the time derivative for h. */
  pi->h_dt -= mj * dvdr / rhoj * wi_dr;                                         // 1xSUB, 2xMUL, 1xDIV                 4
  pj->h_dt -= mi * dvdr / rhoi * wj_dr;                                         // 1xSUB, 2xMUL, 1xDIV                 4

  /* Update the signal velocity. */
  pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig);                              // 1xFMAX                              1
  pj->force.v_sig = fmaxf(pj->force.v_sig, v_sig);                              // 1xFMAX                              1
                                                                                //
                                                                                // 135 FLOPs
}

/**
 * @brief Force loop (Vectorized version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_vec_force(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {

#ifdef VECTORIZE

  vector r, r2, ri;
  vector xi, xj;
  vector hi, hj, hi_inv, hj_inv;
  vector hi2_inv, hj2_inv;
  vector wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, dvdr;
  vector w;
  vector piPOrho2, pjPOrho2, pirho, pjrho, piu, pju;
  vector mi, mj;
  vector f;
  vector dx[3];
  vector vi[3], vj[3];
  vector pia[3], pja[3];
  vector piu_dt, pju_dt;
  vector pih_dt, pjh_dt;
  vector ci, cj, v_sig, vi_sig, vj_sig;
  vector omega_ij, Pi_ij, balsara;
  vector pialpha, pjalpha, alpha_ij, v_sig_u, tc;
  int j, k;

/* Load stuff. */
#if VEC_SIZE == 8
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass,
                 pi[4]->mass, pi[5]->mass, pi[6]->mass, pi[7]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass,
                 pj[4]->mass, pj[5]->mass, pj[6]->mass, pj[7]->mass);
  piPOrho2.v =
      vec_set(pi[0]->force.POrho2, pi[1]->force.POrho2, pi[2]->force.POrho2,
              pi[3]->force.POrho2, pi[4]->force.POrho2, pi[5]->force.POrho2,
              pi[6]->force.POrho2, pi[7]->force.POrho2);
  pjPOrho2.v =
      vec_set(pj[0]->force.POrho2, pj[1]->force.POrho2, pj[2]->force.POrho2,
              pj[3]->force.POrho2, pj[4]->force.POrho2, pj[5]->force.POrho2,
              pj[6]->force.POrho2, pj[7]->force.POrho2);
  pirho.v = vec_set(pi[0]->rho, pi[1]->rho, pi[2]->rho, pi[3]->rho, pi[4]->rho,
                    pi[5]->rho, pi[6]->rho, pi[7]->rho);
  pjrho.v = vec_set(pj[0]->rho, pj[1]->rho, pj[2]->rho, pj[3]->rho, pj[4]->rho,
                    pj[5]->rho, pj[6]->rho, pj[7]->rho);
  piu.v = vec_set(pi[0]->u, pi[1]->u, pi[2]->u, pi[3]->u, pi[4]->u, pi[5]->u,
                  pi[6]->u, pi[7]->u);
  pju.v = vec_set(pj[0]->u, pj[1]->u, pj[2]->u, pj[3]->u, pj[4]->u, pj[5]->u,
                  pj[6]->u, pj[7]->u);
  ci.v =
      vec_set(pi[0]->force.c, pi[1]->force.c, pi[2]->force.c, pi[3]->force.c,
              pi[4]->force.c, pi[5]->force.c, pi[6]->force.c, pi[7]->force.c);
  cj.v =
      vec_set(pj[0]->force.c, pj[1]->force.c, pj[2]->force.c, pj[3]->force.c,
              pj[4]->force.c, pj[5]->force.c, pj[6]->force.c, pj[7]->force.c);
  vi_sig.v = vec_set(pi[0]->force.v_sig, pi[1]->force.v_sig, pi[2]->force.v_sig,
                     pi[3]->force.v_sig, pi[4]->force.v_sig, pi[5]->force.v_sig,
                     pi[6]->force.v_sig, pi[7]->force.v_sig);
  vj_sig.v = vec_set(pj[0]->force.v_sig, pj[1]->force.v_sig, pj[2]->force.v_sig,
                     pj[3]->force.v_sig, pj[4]->force.v_sig, pj[5]->force.v_sig,
                     pj[6]->force.v_sig, pj[7]->force.v_sig);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k],
                      pi[4]->v[k], pi[5]->v[k], pi[6]->v[k], pi[7]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k],
                      pj[4]->v[k], pj[5]->v[k], pj[6]->v[k], pj[7]->v[k]);
  }
  balsara.v =
      vec_set(pi[0]->force.balsara, pi[1]->force.balsara, pi[2]->force.balsara,
              pi[3]->force.balsara, pi[4]->force.balsara, pi[5]->force.balsara,
              pi[6]->force.balsara, pi[7]->force.balsara) +
      vec_set(pj[0]->force.balsara, pj[1]->force.balsara, pj[2]->force.balsara,
              pj[3]->force.balsara, pj[4]->force.balsara, pj[5]->force.balsara,
              pj[6]->force.balsara, pj[7]->force.balsara);
  pialpha.v = vec_set(pi[0]->alpha, pi[1]->alpha, pi[2]->alpha, pi[3]->alpha,
                      pi[4]->alpha, pi[5]->alpha, pi[6]->alpha, pi[7]->alpha);
  pjalpha.v = vec_set(pj[0]->alpha, pj[1]->alpha, pj[2]->alpha, pj[3]->alpha,
                      pj[4]->alpha, pj[5]->alpha, pj[6]->alpha, pj[7]->alpha);
#elif VEC_SIZE == 4
  mi.v = vec_set(pi[0]->mass, pi[1]->mass, pi[2]->mass, pi[3]->mass);
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass);
  piPOrho2.v = vec_set(pi[0]->force.POrho2, pi[1]->force.POrho2,
                       pi[2]->force.POrho2, pi[3]->force.POrho2);
  pjPOrho2.v = vec_set(pj[0]->force.POrho2, pj[1]->force.POrho2,
                       pj[2]->force.POrho2, pj[3]->force.POrho2);
  pirho.v = vec_set(pi[0]->rho, pi[1]->rho, pi[2]->rho, pi[3]->rho);
  pjrho.v = vec_set(pj[0]->rho, pj[1]->rho, pj[2]->rho, pj[3]->rho);
  piu.v = vec_set(pi[0]->u, pi[1]->u, pi[2]->u, pi[3]->u);
  pju.v = vec_set(pj[0]->u, pj[1]->u, pj[2]->u, pj[3]->u);
  ci.v =
      vec_set(pi[0]->force.c, pi[1]->force.c, pi[2]->force.c, pi[3]->force.c);
  cj.v =
      vec_set(pj[0]->force.c, pj[1]->force.c, pj[2]->force.c, pj[3]->force.c);
  vi_sig.v = vec_set(pi[0]->force.v_sig, pi[1]->force.v_sig, pi[2]->force.v_sig,
                     pi[3]->force.v_sig);
  vj_sig.v = vec_set(pj[0]->force.v_sig, pj[1]->force.v_sig, pj[2]->force.v_sig,
                     pj[3]->force.v_sig);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k]);
  }
  balsara.v = vec_set(pi[0]->force.balsara, pi[1]->force.balsara,
                      pi[2]->force.balsara, pi[3]->force.balsara) +
              vec_set(pj[0]->force.balsara, pj[1]->force.balsara,
                      pj[2]->force.balsara, pj[3]->force.balsara);
  pialpha.v = vec_set(pi[0]->alpha, pi[1]->alpha, pi[2]->alpha, pi[3]->alpha);
  pjalpha.v = vec_set(pj[0]->alpha, pj[1]->alpha, pj[2]->alpha, pj[3]->alpha);
#else
#error
#endif
  for (k = 0; k < 3; k++)
    dx[k].v = vec_load(Dx + VEC_SIZE * k);

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri.v = vec_rsqrt(r2.v);                                                       // 8xSQRT, 8xRECIP                    16
  ri.v = ri.v - vec_set1(0.5f) * ri.v * (r2.v * ri.v * ri.v - vec_set1(1.0f));  // 2x8xSUB, 4x8xMUL                   48
  r.v = r2.v * ri.v;                                                            // 8xMUL                               8

  /* Get the kernel for hi. */
  hi.v = vec_load(Hi);
  hi_inv.v = vec_rcp(hi.v);                                                     // 8xRECIP                             8
  hi_inv.v = hi_inv.v - hi_inv.v * (hi.v * hi_inv.v - vec_set1(1.0f));          // 2x8xSUB, 2x8xMUL                   32
  hi2_inv.v = hi_inv.v * hi_inv.v;                                              // 8xMUL                               8
  xi.v = r.v * hi_inv.v;                                                        // 8xMUL                               8
  kernel_deval_vec(&xi, &wi, &wi_dx);                                           // 88xFLOP                            88
  wi_dr.v = hi2_inv.v * hi2_inv.v * wi_dx.v;                                    // 2x8xMUL                            16

  /* Get the kernel for hj. */
  hj.v = vec_load(Hj);
  hj_inv.v = vec_rcp(hj.v);                                                     // 8xRECIP                             8
  hj_inv.v = hj_inv.v - hj_inv.v * (hj.v * hj_inv.v - vec_set1(1.0f));          // 2x8xSUB, 2x8xMUL                   32
  hj2_inv.v = hj_inv.v * hj_inv.v;                                              // 8xMUL                               8
  xj.v = r.v * hj_inv.v;                                                        // 8xMUL                               8
  kernel_deval_vec(&xj, &wj, &wj_dx);                                           // 88xFLOP                            88
  wj_dr.v = hj2_inv.v * hj2_inv.v * wj_dx.v;                                    // 2x8xMUL                            16

  /* Compute dv dot r. */
  dvdr.v = ((vi[0].v - vj[0].v) * dx[0].v) + ((vi[1].v - vj[1].v) * dx[1].v) +  // 2x8xSUB, 2x8xADD, 2x8xMUL          48
           ((vi[2].v - vj[2].v) * dx[2].v);                                     // 8xSUB, 8xADD                       16
  dvdr.v = dvdr.v * ri.v;

  /* Get the time derivative for h. */
  pih_dt.v = mj.v / pjrho.v * dvdr.v * wi_dr.v;                                 // 8xDIV, 2x8xMUL                     24
  pjh_dt.v = mi.v / pirho.v * dvdr.v * wj_dr.v;                                 // 8xDIV, 2x8xMUL                     24

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij.v = vec_fmin(dvdr.v, vec_set1(0.0f));                                // 8xFMIN                              8

  /* Compute signal velocity */
  v_sig.v = ci.v + cj.v - vec_set1(2.0f) * omega_ij.v;                          // 8xADD, 8xSUB, 8xMUL                24

  /* Compute viscosity parameter */
  alpha_ij.v = vec_set1(-0.5f) * (pialpha.v + pjalpha.v);                       // 8xMUL, 8xADD                       16

  /* Compute viscosity tensor */
  Pi_ij.v = balsara.v * alpha_ij.v * v_sig.v * omega_ij.v / (pirho.v + pjrho.v);// 3x8xMUL, 8xDIV, 8xADD              40
  Pi_ij.v *= (wi_dr.v + wj_dr.v);                                               // 8xMUL, 8xADD                       16

  /* Thermal conductivity */
  v_sig_u.v = vec_sqrt(vec_set1(2.f * (const_hydro_gamma - 1.f)) *              // 8xSQRT, 1xMUL, 1xSUB               10
                       vec_fabs(pirho.v * piu.v - pjrho.v * pju.v) /            // 8xFABS, 2x8xMUL, 8xSUB, 8xDIV      40
                       (pirho.v + pjrho.v));                                    // 8xADD                               8
  tc.v = vec_set1(const_conductivity_alpha) * v_sig_u.v / (pirho.v + pjrho.v);  // 8xMUL, 8xDIV, 8xADD                24
  tc.v *= (wi_dr.v + wj_dr.v);                                                  // 8xMUL, 8xADD                       16

  /* Get the common factor out. */
  w.v = ri.v * ((piPOrho2.v * wi_dr.v + pjPOrho2.v * wj_dr.v) +                 // 3x8xMUL, 2x8xADD                   40
                vec_set1(0.25f) * Pi_ij.v);                                     // 8xMUL                               8

  /* Use the force, Luke! */
  for (k = 0; k < 3; k++) {                                                     // 3x...                              48
    f.v = dx[k].v * w.v;                                                        //   8xMUL
    pia[k].v = mj.v * f.v;                                                      //   8xMUL
    pja[k].v = mi.v * f.v;
  }

  /* Get the time derivative for u. */
  piu_dt.v =
      mj.v * dvdr.v * (piPOrho2.v * wi_dr.v + vec_set1(0.125f) * Pi_ij.v);      // 4x8xMUL, 8xADD                     40
  pju_dt.v =
      mi.v * dvdr.v * (pjPOrho2.v * wj_dr.v + vec_set1(0.125f) * Pi_ij.v);      // 4x8xMUL, 8xADD                     40

  /* Add the thermal conductivity */
  piu_dt.v += mj.v * tc.v * (piu.v - pju.v);                                    // 8xADD, 2x8xMUL, 8xSUB              32
  pju_dt.v += mi.v * tc.v * (pju.v - piu.v);                                    // 8xADD, 2x8xMUL, 8xSUB              32

  /* compute the signal velocity (this is always symmetrical). */
  vi_sig.v = vec_fmax(vi_sig.v, v_sig.v);                                       // 8xFMAX                              8
  vj_sig.v = vec_fmax(vj_sig.v, v_sig.v);                                       // 8xFMAX                              8

  /* Store the forces back on the particles. */
  for (k = 0; k < VEC_SIZE; k++) {                                              // 8x...                              80
    pi[k]->force.u_dt += piu_dt.f[k];                                           //   1xADD
    pj[k]->force.u_dt += pju_dt.f[k];                                           //   1xADD
    pi[k]->h_dt -= pih_dt.f[k];                                                 //   1xSUB
    pj[k]->h_dt -= pjh_dt.f[k];                                                 //   1xSUB
    pi[k]->force.v_sig = vi_sig.f[k];
    pj[k]->force.v_sig = vj_sig.f[k];
    for (j = 0; j < 3; j++) {                                                   //   3x...
      pi[k]->a_hydro[j] -= pia[j].f[k];                                         //     1xSUB
      pj[k]->a_hydro[j] += pja[j].f[k];                                         //     1xADD
    }
  }
                                                                                //
                                                                                // 978 FLOPs

#else

  for (int k = 0; k < VEC_SIZE; k++)
    runner_iact_force(R2[k], &Dx[3 * k], Hi[k], Hj[k], pi[k], pj[k]);

#endif
}

/**
 * @brief Force loop (non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float r = sqrtf(r2), ri = 1.0f / r;
  float xi, xj;
  float hi_inv, hi2_inv;
  float hj_inv, hj2_inv;
  float wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, w, dvdr;
  float /*mi,*/ mj, POrho2i, POrho2j, rhoi, rhoj;
  float v_sig, omega_ij, Pi_ij, alpha_ij, tc, v_sig_u;
  float f;
  int k;

  /* Get some values in local variables. */
  // mi = pi->mass;
  mj = pj->mass;
  rhoi = pi->rho;
  rhoj = pj->rho;
  POrho2i = pi->force.POrho2;
  POrho2j = pj->force.POrho2;

  /* Get the kernel for hi. */
  hi_inv = 1.0f / hi;                                                           // 1xRECIP                             1
  hi2_inv = hi_inv * hi_inv;                                                    // 1xMUL                               1
  xi = r * hi_inv;                                                              // 1xMUL                               1
  kernel_deval(xi, &wi, &wi_dx);                                                // 11xFLOP                            11
  wi_dr = hi2_inv * hi2_inv * wi_dx;                                            // 2xMUL                               2

  /* Get the kernel for hj. */
  hj_inv = 1.0f / hj;                                                           // 1xRECIP                             1
  hj2_inv = hj_inv * hj_inv;                                                    // 1xMUL                               1
  xj = r * hj_inv;                                                              // 1xMUL                               1
  kernel_deval(xj, &wj, &wj_dx);                                                // 11xFLOP                            11
  wj_dr = hj2_inv * hj2_inv * wj_dx;                                            // 2xMUL                               2

  /* Compute dv dot r. */
  dvdr = (pi->v[0] - pj->v[0]) * dx[0] + (pi->v[1] - pj->v[1]) * dx[1] +        // 2xSUB, 2xMUL, 2xADD                 6
         (pi->v[2] - pj->v[2]) * dx[2];                                         // 1xSUB, 1xMUL                        2
  dvdr *= ri;                                                                   // 1xMUL                               1

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij = fminf(dvdr, 0.f);                                                  // 1xFMIN                              1

  /* Compute signal velocity */
  v_sig = pi->force.c + pj->force.c - 2.0f * omega_ij;                          // 1xADD, 1xSUB, 1xMUL                 3

  /* Compute viscosity parameter */
  alpha_ij = -0.5f * (pi->alpha + pj->alpha);                                   // 1xMUL, 1xADD                        2

  /* Compute viscosity tensor */
  Pi_ij = alpha_ij * v_sig * omega_ij / (rhoi + rhoj);                          // 2xMUL, 1xDIV, 1xADD                 4

  /* Apply balsara switch */
  Pi_ij *= (pi->force.balsara + pj->force.balsara);                             // 1xMUL, 1xADD                        2

  /* Thermal conductivity */
  v_sig_u = sqrtf(2.f * (const_hydro_gamma - 1.f) *                             // 1xSQRT, 2xMUL, 1xSUB                4
                  fabs(rhoi * pi->u - rhoj * pj->u) / (rhoi + rhoj));           // 1xFABS, 2xMUL, 1xSUB, 1xDIV, 1xADD  6
  tc = const_conductivity_alpha * v_sig_u / (rhoi + rhoj);                      // 1xMUL, 1xDIV, 1xADD                 3
  tc *= (wi_dr + wj_dr);                                                        // 1xMUL, 1xADD                        2

  /* Get the common factor out. */
  w = ri *                                                                      // 1xMUL                               1
      ((POrho2i * wi_dr + POrho2j * wj_dr) + 0.25f * Pi_ij * (wi_dr + wj_dr));  // 4xMUL, 3xADD                        7

  /* Use the force, Luke! */
  for (k = 0; k < 3; k++) {                                                     // 3x...                               9
    f = dx[k] * w;                                                              //   1xMUL
    pi->a_hydro[k] -= mj * f;                                                   //   1xSUB, 1xMUL
  }

  /* Get the time derivative for u. */
  pi->force.u_dt +=                                                             // 1xADD                               1
      mj * dvdr * (POrho2i * wi_dr + 0.125f * Pi_ij * (wi_dr + wj_dr));         // 5xMUL, 2xADD                        7

  /* Add the thermal conductivity */
  pi->force.u_dt += mj * tc * (pi->u - pj->u);                                  // 1xADD, 2xMUL, 1xSUB                 4

  /* Get the time derivative for h. */
  pi->h_dt -= mj * dvdr / rhoj * wi_dr;                                         // 1xSUB, 2xMUL, 1xDIV                 4

  /* Update the signal velocity. */
  pi->force.v_sig = fmaxf(pi->force.v_sig, v_sig);                              // 1xFMAX                              1
  pj->force.v_sig = fmaxf(pj->force.v_sig, v_sig);                              // 1xFMAX                              1
                                                                                //
                                                                                // 103 FLOPs
}

/**
 * @brief Force loop (Vectorized non-symmetric version)
 */

static void runner_iact_nonsym_vec_force(
    float *R2, float *Dx, float *Hi, float *Hj, struct part **pi,
    struct part **pj) {

#ifdef VECTORIZE

  vector r, r2, ri;
  vector xi, xj;
  vector hi, hj, hi_inv, hj_inv;
  vector hi2_inv, hj2_inv;
  vector wi, wj, wi_dx, wj_dx, wi_dr, wj_dr, dvdr;
  vector w;
  vector piPOrho2, pjPOrho2, pirho, pjrho, piu, pju;
  vector mj;
  vector f;
  vector dx[3];
  vector vi[3], vj[3];
  vector pia[3];
  vector piu_dt;
  vector pih_dt;
  vector ci, cj, v_sig, vi_sig, vj_sig;
  vector omega_ij, Pi_ij, balsara;
  vector pialpha, pjalpha, alpha_ij, v_sig_u, tc;
  int j, k;

/* Load stuff. */
#if VEC_SIZE == 8
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass,
                 pj[4]->mass, pj[5]->mass, pj[6]->mass, pj[7]->mass);
  piPOrho2.v =
      vec_set(pi[0]->force.POrho2, pi[1]->force.POrho2, pi[2]->force.POrho2,
              pi[3]->force.POrho2, pi[4]->force.POrho2, pi[5]->force.POrho2,
              pi[6]->force.POrho2, pi[7]->force.POrho2);
  pjPOrho2.v =
      vec_set(pj[0]->force.POrho2, pj[1]->force.POrho2, pj[2]->force.POrho2,
              pj[3]->force.POrho2, pj[4]->force.POrho2, pj[5]->force.POrho2,
              pj[6]->force.POrho2, pj[7]->force.POrho2);
  pirho.v = vec_set(pi[0]->rho, pi[1]->rho, pi[2]->rho, pi[3]->rho, pi[4]->rho,
                    pi[5]->rho, pi[6]->rho, pi[7]->rho);
  pjrho.v = vec_set(pj[0]->rho, pj[1]->rho, pj[2]->rho, pj[3]->rho, pj[4]->rho,
                    pj[5]->rho, pj[6]->rho, pj[7]->rho);
  piu.v = vec_set(pi[0]->u, pi[1]->u, pi[2]->u, pi[3]->u, pi[4]->u, pi[5]->u,
                  pi[6]->u, pi[7]->u);
  pju.v = vec_set(pj[0]->u, pj[1]->u, pj[2]->u, pj[3]->u, pj[4]->u, pj[5]->u,
                  pj[6]->u, pj[7]->u);
  ci.v =
      vec_set(pi[0]->force.c, pi[1]->force.c, pi[2]->force.c, pi[3]->force.c,
              pi[4]->force.c, pi[5]->force.c, pi[6]->force.c, pi[7]->force.c);
  cj.v =
      vec_set(pj[0]->force.c, pj[1]->force.c, pj[2]->force.c, pj[3]->force.c,
              pj[4]->force.c, pj[5]->force.c, pj[6]->force.c, pj[7]->force.c);
  vi_sig.v = vec_set(pi[0]->force.v_sig, pi[1]->force.v_sig, pi[2]->force.v_sig,
                     pi[3]->force.v_sig, pi[4]->force.v_sig, pi[5]->force.v_sig,
                     pi[6]->force.v_sig, pi[7]->force.v_sig);
  vj_sig.v = vec_set(pj[0]->force.v_sig, pj[1]->force.v_sig, pj[2]->force.v_sig,
                     pj[3]->force.v_sig, pj[4]->force.v_sig, pj[5]->force.v_sig,
                     pj[6]->force.v_sig, pj[7]->force.v_sig);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k],
                      pi[4]->v[k], pi[5]->v[k], pi[6]->v[k], pi[7]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k],
                      pj[4]->v[k], pj[5]->v[k], pj[6]->v[k], pj[7]->v[k]);
  }
  balsara.v =
      vec_set(pi[0]->force.balsara, pi[1]->force.balsara, pi[2]->force.balsara,
              pi[3]->force.balsara, pi[4]->force.balsara, pi[5]->force.balsara,
              pi[6]->force.balsara, pi[7]->force.balsara) +
      vec_set(pj[0]->force.balsara, pj[1]->force.balsara, pj[2]->force.balsara,
              pj[3]->force.balsara, pj[4]->force.balsara, pj[5]->force.balsara,
              pj[6]->force.balsara, pj[7]->force.balsara);
  pialpha.v = vec_set(pi[0]->alpha, pi[1]->alpha, pi[2]->alpha, pi[3]->alpha,
                      pi[4]->alpha, pi[5]->alpha, pi[6]->alpha, pi[7]->alpha);
  pjalpha.v = vec_set(pj[0]->alpha, pj[1]->alpha, pj[2]->alpha, pj[3]->alpha,
                      pj[4]->alpha, pj[5]->alpha, pj[6]->alpha, pj[7]->alpha);
#elif VEC_SIZE == 4
  mj.v = vec_set(pj[0]->mass, pj[1]->mass, pj[2]->mass, pj[3]->mass);
  piPOrho2.v = vec_set(pi[0]->force.POrho2, pi[1]->force.POrho2,
                       pi[2]->force.POrho2, pi[3]->force.POrho2);
  pjPOrho2.v = vec_set(pj[0]->force.POrho2, pj[1]->force.POrho2,
                       pj[2]->force.POrho2, pj[3]->force.POrho2);
  pirho.v = vec_set(pi[0]->rho, pi[1]->rho, pi[2]->rho, pi[3]->rho);
  pjrho.v = vec_set(pj[0]->rho, pj[1]->rho, pj[2]->rho, pj[3]->rho);
  piu.v = vec_set(pi[0]->u, pi[1]->u, pi[2]->u, pi[3]->u);
  pju.v = vec_set(pj[0]->u, pj[1]->u, pj[2]->u, pj[3]->u);
  ci.v =
      vec_set(pi[0]->force.c, pi[1]->force.c, pi[2]->force.c, pi[3]->force.c);
  cj.v =
      vec_set(pj[0]->force.c, pj[1]->force.c, pj[2]->force.c, pj[3]->force.c);
  vi_sig.v = vec_set(pi[0]->force.v_sig, pi[1]->force.v_sig, pi[2]->force.v_sig,
                     pi[3]->force.v_sig);
  vj_sig.v = vec_set(pj[0]->force.v_sig, pj[1]->force.v_sig, pj[2]->force.v_sig,
                     pj[3]->force.v_sig);
  for (k = 0; k < 3; k++) {
    vi[k].v = vec_set(pi[0]->v[k], pi[1]->v[k], pi[2]->v[k], pi[3]->v[k]);
    vj[k].v = vec_set(pj[0]->v[k], pj[1]->v[k], pj[2]->v[k], pj[3]->v[k]);
  }
  balsara.v = vec_set(pi[0]->force.balsara, pi[1]->force.balsara,
                      pi[2]->force.balsara, pi[3]->force.balsara) +
              vec_set(pj[0]->force.balsara, pj[1]->force.balsara,
                      pj[2]->force.balsara, pj[3]->force.balsara);
  pialpha.v = vec_set(pi[0]->alpha, pi[1]->alpha, pi[2]->alpha, pi[3]->alpha);
  pjalpha.v = vec_set(pj[0]->alpha, pj[1]->alpha, pj[2]->alpha, pj[3]->alpha);
#else
#error
#endif
  for (k = 0; k < 3; k++)
    dx[k].v = vec_load(Dx + VEC_SIZE * k);

  /* Get the radius and inverse radius. */
  r2.v = vec_load(R2);
  ri.v = vec_rsqrt(r2.v);                                                       // 8xSQRT, 8xRECIP                    16
  ri.v = ri.v - vec_set1(0.5f) * ri.v * (r2.v * ri.v * ri.v - vec_set1(1.0f));  // 2x8xSUB, 4x8xMUL                   48
  r.v = r2.v * ri.v;                                                            // 8xMUL                               8

  /* Get the kernel for hi. */
  hi.v = vec_load(Hi);
  hi_inv.v = vec_rcp(hi.v);                                                     // 8xRECIP                             8
  hi_inv.v = hi_inv.v - hi_inv.v * (hi.v * hi_inv.v - vec_set1(1.0f));          // 2x8xSUB, 2x8xMUL                   32
  hi2_inv.v = hi_inv.v * hi_inv.v;                                              // 8xMUL                               8
  xi.v = r.v * hi_inv.v;                                                        // 8xMUL                               8
  kernel_deval_vec(&xi, &wi, &wi_dx);                                           // 88xFLOP                            88
  wi_dr.v = hi2_inv.v * hi2_inv.v * wi_dx.v;                                    // 2x8xMUL                            16

  /* Get the kernel for hj. */
  hj.v = vec_load(Hj);
  hj_inv.v = vec_rcp(hj.v);                                                     // 8xRECIP                             8
  hj_inv.v = hj_inv.v - hj_inv.v * (hj.v * hj_inv.v - vec_set1(1.0f));          // 2x8xSUB, 2x8xMUL                   32
  hj2_inv.v = hj_inv.v * hj_inv.v;                                              // 8xMUL                               8
  xj.v = r.v * hj_inv.v;                                                        // 8xMUL                               8
  kernel_deval_vec(&xj, &wj, &wj_dx);                                           // 88xFLOP                            88
  wj_dr.v = hj2_inv.v * hj2_inv.v * wj_dx.v;                                    // 2x8xMUL                            16

  /* Compute dv dot r. */
  dvdr.v = ((vi[0].v - vj[0].v) * dx[0].v) + ((vi[1].v - vj[1].v) * dx[1].v) +  // 2x8xSUB, 2x8xMUL, 2x8xADD          48
           ((vi[2].v - vj[2].v) * dx[2].v);                                     // 8xSUB, 8xMUL                       16
  dvdr.v = dvdr.v * ri.v;                                                       // 8xMUL                               8

  /* Get the time derivative for h. */
  pih_dt.v = mj.v / pjrho.v * dvdr.v * wi_dr.v;                                 // 8xDIV, 2x8xMUL                     24

  /* Compute the relative velocity. (This is 0 if the particles move away from
   * each other and negative otherwise) */
  omega_ij.v = vec_fmin(dvdr.v, vec_set1(0.0f));                                // 8xFMIN                              8

  /* Compute signal velocity */
  v_sig.v = ci.v + cj.v - vec_set1(2.0f) * omega_ij.v;                          // 8xADD, 8xSUB, 8xMUL                24

  /* Compute viscosity parameter */
  alpha_ij.v = vec_set1(-0.5f) * (pialpha.v + pjalpha.v);                       // 8xMUL, 8xADD                       16

  /* Compute viscosity tensor */
  Pi_ij.v = balsara.v * alpha_ij.v * v_sig.v * omega_ij.v / (pirho.v + pjrho.v);// 3x8xMUL, 8xDIV, 8xADD              40
  Pi_ij.v *= (wi_dr.v + wj_dr.v);                                               // 8xMUL, 8xADD                       16

  /* Thermal conductivity */
  v_sig_u.v = vec_sqrt(vec_set1(2.f * (const_hydro_gamma - 1.f)) *              // 8xSQRT, 9xMUL, 1xSUB               18
                       vec_fabs(pirho.v * piu.v - pjrho.v * pju.v) /            // 8xFABS, 2x8xMUL, 8xSUB, 8xDIV      40
                       (pirho.v + pjrho.v));                                    // 8xADD                               8
  tc.v = vec_set1(const_conductivity_alpha) * v_sig_u.v / (pirho.v + pjrho.v);  // 8xMUL, 8xDIV, 8xADD                24
  tc.v *= (wi_dr.v + wj_dr.v);                                                  // 8xMUL, 8xADD                       16

  /* Get the common factor out. */
  w.v = ri.v * ((piPOrho2.v * wi_dr.v + pjPOrho2.v * wj_dr.v) +                 // 3x8xMUL, 2x8xADD                   40
                vec_set1(0.25f) * Pi_ij.v);                                     // 8xMUL                               8

  /* Use the force, Luke! */
  for (k = 0; k < 3; k++) {                                                     // 3x...                              48
    f.v = dx[k].v * w.v;                                                        //   8xMUL
    pia[k].v = mj.v * f.v;                                                      //   8xMUL
  }

  /* Get the time derivative for u. */
  piu_dt.v =
      mj.v * dvdr.v * (piPOrho2.v * wi_dr.v + vec_set1(0.125f) * Pi_ij.v);      // 4x8xMUL, 8xADD                     40

  /* Add the thermal conductivity */
  piu_dt.v += mj.v * tc.v * (piu.v - pju.v);                                    // 8xADD, 2x8xMUL, 8xSUB              32

  /* compute the signal velocity (this is always symmetrical). */
  vi_sig.v = vec_fmax(vi_sig.v, v_sig.v);                                       // 8xFMAX                              8
  vj_sig.v = vec_fmax(vj_sig.v, v_sig.v);                                       // 8xFMAX                              8

  /* Store the forces back on the particles. */
  for (k = 0; k < VEC_SIZE; k++) {                                              // 8x...                              24
    pi[k]->force.u_dt += piu_dt.f[k];                                           //   1xADD
    pi[k]->h_dt -= pih_dt.f[k];                                                 //   1xSUB
    pi[k]->force.v_sig = vi_sig.f[k];
    pj[k]->force.v_sig = vj_sig.f[k];
    for (j = 0; j < 3; j++) pi[k]->a_hydro[j] -= pia[j].f[k];                   //   1xSUB
  }
                                                                                //
                                                                                // 906 FLOPs

#else

  for (int k = 0; k < VEC_SIZE; k++)
    runner_iact_nonsym_force(R2[k], &Dx[3 * k], Hi[k], Hj[k], pi[k], pj[k]);

#endif
}

#endif /* SWIFT_RUNNER_IACT_H */
