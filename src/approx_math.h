/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_APPROX_MATH_H
#define SWIFT_APPROX_MATH_H

/**
 * @brief Approximate version of expf(x) using a 4th order Taylor expansion
 *
 * The absolute error is of order 10^-6 for -0.2 < x < 0.2.
 *
 * @param x The number to take the exponential of.
 */
__attribute__((always_inline)) INLINE static float approx_expf(float x) {
  return 1.f + x * (1.f + x * (0.5f + x * (1.f / 6.0f + 1.f / 24.0f * x))); // 4xADD, 4xMUL, 2xDIV
                                                                            //
                                                                            // 10 FLOPs
}

#endif /* SWIFT_APPROX_MATH_H */
