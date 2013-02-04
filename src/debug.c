/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2012 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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


#include <stdio.h>

#include "part.h"

void printParticle(struct part *parts, int i)
{
  printf("## Particle[%d]: id= %lld x=( %f, %f, %f) v=( %f, %f, %f) h= %f m= %f rho= %f u= %f dt= %f\n",
	 i,
	 parts[i].id,
	 parts[i].x[0], parts[i].x[1], parts[i].x[2],
	 parts[i].v[0], parts[i].v[1], parts[i].v[2],
	 parts[i].h,
	 parts[i].mass,
	 parts[i].rho,
	 parts[i].u,
	 parts[i].dt
	 );
}
