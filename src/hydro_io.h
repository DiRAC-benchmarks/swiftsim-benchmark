/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_HYDRO_IO_H
#define SWIFT_HYDRO_IO_H

#include "./const.h"

/* Import the right functions */
#if defined(MINIMAL_SPH)
#include "./hydro/Minimal/hydro_io.h"
#elif defined(GADGET2_SPH)
#include "./hydro/Gadget2/hydro_io.h"
#elif defined(DEFAULT_SPH)
#include "./hydro/Default/hydro_io.h"
#else
#error "Invalid choice of SPH variant"
#endif

#endif /* SWIFT_HYDRO_IO_H */
