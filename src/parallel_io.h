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


#if defined(HAVE_HDF5) && defined(WITH_MPI)

void read_ic_parallel ( char* fileName, double dim[3], struct part **parts,  int* N, int* periodic, int mpi_rank, int mpi_size, MPI_Comm comm, MPI_Info info);

void write_output_parallel ( struct engine* e, int mpi_rank, int mpi_size, MPI_Comm comm, MPI_Info info);

#endif
