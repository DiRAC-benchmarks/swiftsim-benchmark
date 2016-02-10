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
#ifndef SWIFT_RUNNER_IACT_LEGACY_H
#define SWIFT_RUNNER_IACT_LEGACY_H

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

  float wi, wi_dx;
  float wj, wj_dx;
  float dv[3], curlvr[3];

  /* Get the masses. */
  const float mi = pj->mass;
  const float mj = pj->mass;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Compute the kernel function for pi */
  const float hi_inv = 1.f / hi;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);

  /* Compute contribution to the density */
  pi->rho += mj * wi;
  pi->rho_dh -= mj * kernel_igamma * (3.f * wi + ui * wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= ui * wi_dx;

  /* Compute the kernel function for pj */
  const float hj_inv = 1.f / hj;
  const float uj = r * hj_inv;
  kernel_deval(uj, &wj, &wj_dx);

  /* Compute contribution to the density */
  pj->rho += mi * wj;
  pj->rho_dh -= mi * kernel_igamma * (3.f * wj + uj * wj_dx);

  /* Compute contribution to the number of neighbours */
  pj->density.wcount += wj;
  pj->density.wcount_dh -= uj * wj_dx;

  const float faci = mj * wi_dx * r_inv;
  const float facj = mi * wj_dx * r_inv;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];

  pi->div_v += faci * dvdr;
  pj->div_v += facj * dvdr;

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->rot_v[0] += faci * curlvr[0];
  pi->rot_v[1] += faci * curlvr[1];
  pi->rot_v[2] += faci * curlvr[2];

  pj->rot_v[0] += facj * curlvr[0];
  pj->rot_v[1] += facj * curlvr[1];
  pj->rot_v[2] += facj * curlvr[2];
}

/**
 * @brief Density loop (non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_density(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wi_dx;
  float dv[3], curlvr[3];

  /* Get the masses. */
  const float mj = pj->mass;

  /* Get r and r inverse. */
  const float r = sqrtf(r2);
  const float ri = 1.0f / r;

  /* Compute the kernel function */
  const float h_inv = 1.0f / hi;
  const float u = r * h_inv;
  kernel_deval(u, &wi, &wi_dx);

  /* Compute contribution to the density */
  pi->rho += mj * wi;
  pi->rho_dh -= mj * kernel_igamma * (3.f * wi + u * wi_dx);

  /* Compute contribution to the number of neighbours */
  pi->density.wcount += wi;
  pi->density.wcount_dh -= u * wi_dx;

  /* const float ih3 = h_inv * h_inv * h_inv; */
  /* const float ih4 = h_inv * h_inv * h_inv * h_inv; */

  const float fac = mj * wi_dx * ri;

  /* Compute dv dot r */
  dv[0] = pi->v[0] - pj->v[0];
  dv[1] = pi->v[1] - pj->v[1];
  dv[2] = pi->v[2] - pj->v[2];
  const float dvdr = dv[0] * dx[0] + dv[1] * dx[1] + dv[2] * dx[2];
  pi->div_v -= fac * dvdr;

  /* if(pi->id == 515050 && pj->id == 504849) */
  /*   message("Interacting with %lld. r=%e hi=%e u=%e W=%e dW/dx=%e dh_drho1=%e
   * dh_drho2=%e\n fac=%e dvdr=%e pj->v=[%.3e,%.3e,%.3e]", */
  /* 	    pj->id, */
  /* 	    r, */
  /* 	    hi, */
  /* 	    u, */
  /* 	    wi * ih3, */
  /* 	    wi_dx * ih4, */
  /* 	    -mj * (3.f * kernel_igamma * wi) * ih4, */
  /* 	    -mj * u * wi_dx * kernel_igamma * ih4, */
  /* 	    fac * ih4, */
  /* 	    dvdr, */
  /* 	    pj->v[0], */
  /* 	    pj->v[1], */
  /* 	    pj->v[2] */
  /* 	    ); */

  /* Compute dv cross r */
  curlvr[0] = dv[1] * dx[2] - dv[2] * dx[1];
  curlvr[1] = dv[2] * dx[0] - dv[0] * dx[2];
  curlvr[2] = dv[0] * dx[1] - dv[1] * dx[0];

  pi->rot_v[0] += fac * curlvr[0];
  pi->rot_v[1] += fac * curlvr[1];
  pi->rot_v[2] += fac * curlvr[2];
}

/**
 * @brief Force loop
 */

__attribute__((always_inline)) INLINE static void runner_iact_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wj, wi_dx, wj_dx;

  const float fac_mu = 1.f; /* Will change with cosmological integration */

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Get some values in local variables. */
  // const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float pressurei = pi->pressure;
  const float pressurej = pj->pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hi2_inv = hi_inv * hi_inv;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);
  const float wi_dr = hi2_inv * hi2_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hj2_inv = hj_inv * hj_inv;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hj2_inv * hj2_inv * wj_dx;

  /* Compute gradient terms */
  const float P_over_rho_i = pressurei / (rhoi * rhoi) * pi->rho_dh;
  const float P_over_rho_j = pressurej / (rhoj * rhoj) * pj->rho_dh;

  /* Compute sound speeds */
  const float ci = sqrtf(const_hydro_gamma * pressurei / rhoi);
  const float cj = sqrtf(const_hydro_gamma * pressurej / rhoj);
  float v_sig = ci + cj;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Artificial viscosity term */
  float visc = 0.f;
  if (dvdr < 0.f) {
    const float mu_ij = fac_mu * dvdr * r_inv;
    v_sig -= 3.f * mu_ij;
    const float rho_ij = 0.5f * (rhoi + rhoj);
    const float balsara_i = fabsf(pi->div_v) / (fabsf(pi->div_v) + pi->curl_v +
                                                0.0001 * ci / fac_mu / hi);
    const float balsara_j = fabsf(pj->div_v) / (fabsf(pj->div_v) + pj->curl_v +
                                                0.0001 * cj / fac_mu / hj);
    visc = -0.25f * const_viscosity_alpha * v_sig * mu_ij / rho_ij *
           (balsara_i + balsara_j);
  }

  /* Now, convolve with the kernel */
  const float visc_term = 0.5f * mj * visc * (wi_dr + wj_dr) * r_inv;
  const float sph_term =
      mj * (P_over_rho_i * wi_dr + P_over_rho_j * wj_dr) * r_inv;

  /* Eventually got the acceleration */
  const float acc = visc_term + sph_term;

  /* //if(pi->id == 1000 && pj->id == 1100) */
  /* if(pi->id == 515050 && pj->id == 504849) */
  /*   message("Interacting with %lld. r=%e hi=%e hj=%e dWi/dx=%e dWj/dx=%3e
   * dvdr=%e visc=%e sph=%e", */
  /* 	    pj->id, */
  /* 	    r, */
  /* 	    2*hi, */
  /* 	    2*hj, */
  /* 	    wi_dr, */
  /* 	    wj_dr, */
  /* 	    dvdr, */
  /* 	    visc_term, */
  /* 	    sph_term */
  /* 	    ); */
  /* if(pi->id == 1100 && pj->id == 1000) */
  /*   message("oO"); */

  /* Use the force Luke ! */
  pi->a[0] -= acc * dx[0];
  pi->a[1] -= acc * dx[1];
  pi->a[2] -= acc * dx[2];

  pj->a[0] += acc * dx[0];
  pj->a[1] += acc * dx[1];
  pj->a[2] += acc * dx[2];

  /* Update the signal velocity. */
  pi->v_sig = fmaxf(pi->v_sig, v_sig);
  pj->v_sig = fmaxf(pj->v_sig, v_sig);

  /* Change in entropy */
  pi->entropy_dt += 0.5f * visc_term * dvdr;
  pj->entropy_dt -= 0.5f * visc_term * dvdr;
}

/**
 * @brief Force loop (non-symmetric version)
 */

__attribute__((always_inline)) INLINE static void runner_iact_nonsym_force(
    float r2, float *dx, float hi, float hj, struct part *pi, struct part *pj) {

  float wi, wj, wi_dx, wj_dx;

  const float fac_mu = 1.f; /* Will change with cosmological integration */

  const float r = sqrtf(r2);
  const float r_inv = 1.0f / r;

  /* Get some values in local variables. */
  // const float mi = pi->mass;
  const float mj = pj->mass;
  const float rhoi = pi->rho;
  const float rhoj = pj->rho;
  const float pressurei = pi->pressure;
  const float pressurej = pj->pressure;

  /* Get the kernel for hi. */
  const float hi_inv = 1.0f / hi;
  const float hi2_inv = hi_inv * hi_inv;
  const float ui = r * hi_inv;
  kernel_deval(ui, &wi, &wi_dx);
  const float wi_dr = hi2_inv * hi2_inv * wi_dx;

  /* Get the kernel for hj. */
  const float hj_inv = 1.0f / hj;
  const float hj2_inv = hj_inv * hj_inv;
  const float xj = r * hj_inv;
  kernel_deval(xj, &wj, &wj_dx);
  const float wj_dr = hj2_inv * hj2_inv * wj_dx;

  /* Compute gradient terms */
  const float P_over_rho_i = pressurei / (rhoi * rhoi) * pi->rho_dh;
  const float P_over_rho_j = pressurej / (rhoj * rhoj) * pj->rho_dh;

  /* Compute sound speeds */
  const float ci = sqrtf(const_hydro_gamma * pressurei / rhoi);
  const float cj = sqrtf(const_hydro_gamma * pressurej / rhoj);
  float v_sig = ci + cj;

  /* Compute dv dot r. */
  const float dvdr = (pi->v[0] - pj->v[0]) * dx[0] +
                     (pi->v[1] - pj->v[1]) * dx[1] +
                     (pi->v[2] - pj->v[2]) * dx[2];

  /* Artificial viscosity term */
  float visc = 0.f;
  if (dvdr < 0.f) {
    const float mu_ij = fac_mu * dvdr * r_inv;
    v_sig -= 3.f * mu_ij;
    const float rho_ij = 0.5f * (rhoi + rhoj);
    const float balsara_i = fabsf(pi->div_v) / (fabsf(pi->div_v) + pi->curl_v +
                                                0.0001 * ci / fac_mu / hi);
    const float balsara_j = fabsf(pj->div_v) / (fabsf(pj->div_v) + pj->curl_v +
                                                0.0001 * cj / fac_mu / hj);
    visc = -0.25f * const_viscosity_alpha * v_sig * mu_ij / rho_ij *
           (balsara_i + balsara_j);
  }

  /* Now, convolve with the kernel */
  const float visc_term = 0.5f * mj * visc * (wi_dr + wj_dr) * r_inv;
  const float sph_term =
      mj * (P_over_rho_i * wi_dr + P_over_rho_j * wj_dr) * r_inv;

  /* Eventually got the acceleration */
  const float acc = visc_term + sph_term;

  /* Use the force Luke ! */
  pi->a[0] -= acc * dx[0];
  pi->a[1] -= acc * dx[1];
  pi->a[2] -= acc * dx[2];

  /* Update the signal velocity. */
  pi->v_sig = fmaxf(pi->v_sig, v_sig);

  /* Change in entropy */
  pi->entropy_dt += 0.5f * visc_term * dvdr;
}

#endif /* SWIFT_RUNNER_IACT_LEGACY_H */