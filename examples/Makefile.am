
# This file is part of SWIFT.
# Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
#                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# Common flags
MYFLAGS = -DTIMER

# Add the source directory and debug to CFLAGS
AM_CFLAGS = -I../src $(HDF5_CPPFLAGS)

AM_LDFLAGS =

MPI_THREAD_LIBS = @MPI_THREAD_LIBS@
MPI_LIBS = $(METIS_LIBS) $(MPI_THREAD_LIBS)
MPI_FLAGS = -DWITH_MPI $(METIS_INCS)

# Set-up the library
bin_PROGRAMS = swift swift_fixdt

# Build MPI versions as well?
if HAVEMPI
bin_PROGRAMS += swift_mpi swift_fixdt_mpi
endif

# engine_policy_setaffinity is available?
if HAVESETAFFINITY
ENGINE_POLICY_SETAFFINITY=| engine_policy_setaffinity
else
ENGINE_POLICY_SETAFFINITY=
endif

# Sources for swift
swift_SOURCES = main.c
swift_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) -DENGINE_POLICY="engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
swift_LDADD =  ../src/.libs/libswiftsim.a $(HDF5_LDFLAGS) $(HDF5_LIBS)

swift_fixdt_SOURCES = main.c
swift_fixdt_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) -DENGINE_POLICY="engine_policy_fixdt | engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
swift_fixdt_LDADD =  ../src/.libs/libswiftsim.a $(HDF5_LDFLAGS) $(HDF5_LIBS)

# Sources for swift_mpi, do we need an affinity policy for MPI?
swift_mpi_SOURCES = main.c
swift_mpi_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) $(MPI_FLAGS) -DENGINE_POLICY="engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
swift_mpi_LDADD =  ../src/.libs/libswiftsim_mpi.a $(HDF5_LDFLAGS) $(HDF5_LIBS) $(MPI_LIBS)

swift_fixdt_mpi_SOURCES = main.c
swift_fixdt_mpi_CFLAGS = $(MYFLAGS) $(AM_CFLAGS) $(MPI_FLAGS) -DENGINE_POLICY="engine_policy_fixdt | engine_policy_keep $(ENGINE_POLICY_SETAFFINITY)"
swift_fixdt_mpi_LDADD =  ../src/.libs/libswiftsim_mpi.a $(HDF5_LDFLAGS) $(HDF5_LIBS) $(MPI_LIBS)

# Scripts to generate ICs
EXTRA_DIST = UniformBox/makeIC.py \
	     UniformDMBox/makeIC.py \
	     PerturbedBox/makeIC.py \
	     SedovBlast/makeIC.py SedovBlast/makeIC_fcc.py SedovBlast/solution.py \
	     SodShock/makeIC.py SodShock/solution.py SodShock/glass_001.hdf5 SodShock/glass_002.hdf5 SodShock/rhox.py \
	     CosmoVolume/getIC.sh \
	     BigCosmoVolume/makeIC.py \
	     BigPerturbedBox/makeIC_fcc.py \
             GreshoVortex/makeIC.py GreshoVortex/solution.py \
             MultiTypes/makeIC.py


# Scripts to plot task graphs
EXTRA_DIST += plot_tasks_MPI.py plot_tasks.py \
	      process_plot_tasks_MPI process_plot_tasks

# Simple run scripts
EXTRA_DIST += runs.sh
