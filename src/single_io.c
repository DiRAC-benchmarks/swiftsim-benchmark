/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk).
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

/* Config parameters. */
#include "../config.h"

#if defined(HAVE_HDF5) && !defined(WITH_MPI)

/* Some standard headers. */
#include <hdf5.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* This object's header. */
#include "single_io.h"

/* Local includes. */
#include "const.h"
#include "common_io.h"
#include "error.h"

/*-----------------------------------------------------------------------------
 * Routines reading an IC file
 *-----------------------------------------------------------------------------*/

/**
 * @brief Reads a data array from a given HDF5 group.
 *
 * @param grp The group from which to read.
 * @param name The name of the array to read.
 * @param type The #DATA_TYPE of the attribute.
 * @param N The number of particles.
 * @param dim The dimension of the data (1 for scalar, 3 for vector)
 * @param part_c A (char*) pointer on the first occurrence of the field of
 *interest in the parts array
 * @param partSize The size in bytes of the particle structure.
 * @param importance If COMPULSORY, the data must be present in the IC file. If
 *OPTIONAL, the array will be zeroed when the data is not present.
 *
 * @todo A better version using HDF5 hyper-slabs to read the file directly into
 *the part array
 * will be written once the structures have been stabilized.
 */
void readArrayBackEnd(hid_t grp, char* name, enum DATA_TYPE type, int N,
                      int dim, char* part_c, size_t partSize,
                      enum DATA_IMPORTANCE importance) {
  hid_t h_data = 0, h_err = 0, h_type = 0;
  htri_t exist = 0;
  void* temp;
  int i = 0;
  const size_t typeSize = sizeOfType(type);
  const size_t copySize = typeSize * dim;
  char* temp_c = 0;

  /* Check whether the dataspace exists or not */
  exist = H5Lexists(grp, name, 0);
  if (exist < 0) {
    error("Error while checking the existence of data set '%s'.", name);
  } else if (exist == 0) {
    if (importance == COMPULSORY) {
      error("Compulsory data set '%s' not present in the file.", name);
    } else {
      /* message("Optional data set '%s' not present. Zeroing this particle
       * field...", name);	   */

      for (i = 0; i < N; ++i) memset(part_c + i * partSize, 0, copySize);

      return;
    }
  }

  /* message( "Reading %s '%s' array...", importance == COMPULSORY ?
   * "compulsory": "optional  ", name); */

  /* Open data space */
  h_data = H5Dopen(grp, name, H5P_DEFAULT);
  if (h_data < 0) {
    error("Error while opening data space '%s'.", name);
  }

  /* Check data type */
  h_type = H5Dget_type(h_data);
  if (h_type < 0) error("Unable to retrieve data type from the file");
  // if (!H5Tequal(h_type, hdf5Type(type)))
  //  error("Non-matching types between the code and the file");

  /* Allocate temporary buffer */
  temp = malloc(N * dim * typeSize);
  if (temp == NULL) error("Unable to allocate memory for temporary buffer");

  /* Read HDF5 dataspace in temporary buffer */
  /* Dirty version that happens to work for vectors but should be improved */
  /* Using HDF5 dataspaces would be better */
  h_err = H5Dread(h_data, hdf5Type(type), H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
  if (h_err < 0) {
    error("Error while reading data array '%s'.", name);
  }

  /* Copy temporary buffer to particle data */
  temp_c = temp;
  for (i = 0; i < N; ++i)
    memcpy(part_c + i * partSize, &temp_c[i * copySize], copySize);

  /* Free and close everything */
  free(temp);
  H5Tclose(h_type);
  H5Dclose(h_data);
}

/*-----------------------------------------------------------------------------
 * Routines writing an output file
 *-----------------------------------------------------------------------------*/

/**
 * @brief Writes a data array in given HDF5 group.
 *
 * @param grp The group in which to write.
 * @param fileName The name of the file in which the data is written
 * @param xmfFile The FILE used to write the XMF description
 * @param partTypeGroupName The name of the group containing the particles in
 *the HDF5 file.
 * @param name The name of the array to write.
 * @param type The #DATA_TYPE of the array.
 * @param N The number of particles to write.
 * @param dim The dimension of the data (1 for scalar, 3 for vector)
 * @param part_c A (char*) pointer on the first occurrence of the field of
 *interest in the parts array.
 * @param partSize The size in bytes of the particle structure.
 * @param us The UnitSystem currently in use
 * @param convFactor The UnitConversionFactor for this array
 *
 * @todo A better version using HDF5 hyper-slabs to write the file directly from
 *the part array
 * will be written once the structures have been stabilized.
 */
void writeArrayBackEnd(hid_t grp, char* fileName, FILE* xmfFile,
                       char* partTypeGroupName, char* name, enum DATA_TYPE type,
                       int N, int dim, char* part_c, size_t partSize,
                       struct UnitSystem* us,
                       enum UnitConversionFactor convFactor) {
  hid_t h_data = 0, h_err = 0, h_space = 0, h_prop = 0;
  void* temp = 0;
  int i = 0, rank = 0;
  const size_t typeSize = sizeOfType(type);
  const size_t copySize = typeSize * dim;
  char* temp_c = 0;
  hsize_t shape[2];
  hsize_t chunk_shape[2];
  char buffer[FILENAME_BUFFER_SIZE];

  /* message("Writing '%s' array...", name); */

  /* Allocate temporary buffer */
  temp = malloc(N * dim * sizeOfType(type));
  if (temp == NULL) error("Unable to allocate memory for temporary buffer");

  /* Copy particle data to temporary buffer */
  temp_c = temp;
  for (i = 0; i < N; ++i)
    memcpy(&temp_c[i * copySize], part_c + i * partSize, copySize);

  /* Create data space */
  h_space = H5Screate(H5S_SIMPLE);
  if (h_space < 0) {
    error("Error while creating data space for field '%s'.", name);
  }

  if (dim > 1) {
    rank = 2;
    shape[0] = N;
    shape[1] = dim;
    chunk_shape[0] = 1 << 16; /* Just a guess...*/
    chunk_shape[1] = dim;
  } else {
    rank = 1;
    shape[0] = N;
    shape[1] = 0;
    chunk_shape[0] = 1 << 16; /* Just a guess...*/
    chunk_shape[1] = 0;
  }

  /* Make sure the chunks are not larger than the dataset */
  if (chunk_shape[0] > N) chunk_shape[0] = N;

  /* Change shape of data space */
  h_err = H5Sset_extent_simple(h_space, rank, shape, NULL);
  if (h_err < 0) {
    error("Error while changing data space shape for field '%s'.", name);
  }

  /* Dataset properties */
  h_prop = H5Pcreate(H5P_DATASET_CREATE);

  /* Set chunk size */
  h_err = H5Pset_chunk(h_prop, rank, chunk_shape);
  if (h_err < 0) {
    error("Error while setting chunk size (%lld, %lld) for field '%s'.",
          chunk_shape[0], chunk_shape[1], name);
  }

  /* Impose data compression */
  h_err = H5Pset_deflate(h_prop, 4);
  if (h_err < 0) {
    error("Error while setting compression options for field '%s'.", name);
  }

  /* Create dataset */
  h_data = H5Dcreate(grp, name, hdf5Type(type), h_space, H5P_DEFAULT, h_prop,
                     H5P_DEFAULT);
  if (h_data < 0) {
    error("Error while creating dataspace '%s'.", name);
  }

  /* Write temporary buffer to HDF5 dataspace */
  h_err = H5Dwrite(h_data, hdf5Type(type), h_space, H5S_ALL, H5P_DEFAULT, temp);
  if (h_err < 0) {
    error("Error while writing data array '%s'.", name);
  }

  /* Write XMF description for this data set */
  writeXMFline(xmfFile, fileName, partTypeGroupName, name, N, dim, type);

  /* Write unit conversion factors for this data set */
  conversionString(buffer, us, convFactor);
  writeAttribute_d(h_data, "CGS conversion factor",
                   conversionFactor(us, convFactor));
  writeAttribute_f(h_data, "h-scale exponent", hFactor(us, convFactor));
  writeAttribute_f(h_data, "a-scale exponent", aFactor(us, convFactor));
  writeAttribute_s(h_data, "Conversion factor", buffer);

  /* Free and close everything */
  free(temp);
  H5Pclose(h_prop);
  H5Dclose(h_data);
  H5Sclose(h_space);
}

/**
 * @brief A helper macro to call the readArrayBackEnd function more easily.
 *
 * @param grp The group from which to read.
 * @param name The name of the array to read.
 * @param type The #DATA_TYPE of the attribute.
 * @param N The number of particles.
 * @param dim The dimension of the data (1 for scalar, 3 for vector)
 * @param part The array of particles to fill
 * @param N_total Unused parameter in non-MPI mode
 * @param offset Unused parameter in non-MPI mode
 * @param field The name of the field (C code name as defined in part.h) to fill
 * @param importance Is the data compulsory or not
 *
 */
#define readArray(grp, name, type, N, dim, part, N_total, offset, field, \
                  importance)                                            \
  readArrayBackEnd(grp, name, type, N, dim, (char*)(&(part[0]).field),   \
                   sizeof(part[0]), importance)

/**
 * @brief A helper macro to call the readArrayBackEnd function more easily.
 *
 * @param grp The group in which to write.
 * @param fileName The name of the file in which the data is written
 * @param xmfFile The FILE used to write the XMF description
 * @param name The name of the array to write.
 * @param partTypeGroupName The name of the group containing the particles in
 *the HDF5 file.
 * @param type The #DATA_TYPE of the array.
 * @param N The number of particles to write.
 * @param dim The dimension of the data (1 for scalar, 3 for vector)
 * @param part A (char*) pointer on the first occurrence of the field of
 * interest in the parts array
 * @param N_total Unused parameter in non-MPI mode
 * @param mpi_rank Unused parameter in non-MPI mode
 * @param offset Unused parameter in non-MPI mode
 * @param field The name (code name) of the field to read from.
 * @param us The UnitSystem currently in use
 * @param convFactor The UnitConversionFactor for this array
 *
 */
#define writeArray(grp, fileName, xmfFile, partTypeGroupName, name, type, N,  \
                   dim, part, N_total, mpi_rank, offset, field, us,           \
                   convFactor)                                                \
  writeArrayBackEnd(grp, fileName, xmfFile, partTypeGroupName, name, type, N, \
                    dim, (char*)(&(part[0]).field), sizeof(part[0]), us,      \
                    convFactor)

/* Import the right hydro definition */
#include "hydro_io.h"
/* Import the right gravity definition */
#include "gravity_io.h"

/**
 * @brief Reads an HDF5 initial condition file (GADGET-3 type)
 *
 * @param fileName The file to read.
 * @param dim (output) The dimension of the volume.
 * @param parts (output) Array of Gas particles.
 * @param gparts (output) Array of #gpart particles.
 * @param Ngas (output) number of Gas particles read.
 * @param Ngparts (output) The number of #gpart read.
 * @param periodic (output) 1 if the volume is periodic, 0 if not.
 *
 * Opens the HDF5 file fileName and reads the particles contained
 * in the parts array. N is the returned number of particles found
 * in the file.
 *
 * @warning Can not read snapshot distributed over more than 1 file !!!
 * @todo Read snapshots distributed in more than one file.
 *
 * Calls #error() if an error occurs.
 *
 */
void read_ic_single(char* fileName, double dim[3], struct part** parts,
                    struct gpart** gparts, size_t* Ngas, size_t* Ngparts,
                    int* periodic) {
  hid_t h_file = 0, h_grp = 0;
  /* GADGET has only cubic boxes (in cosmological mode) */
  double boxSize[3] = {0.0, -1.0, -1.0};
  /* GADGET has 6 particle types. We only keep the type 0 & 1 for now...*/
  int numParticles[NUM_PARTICLE_TYPES] = {0};
  int numParticles_highWord[NUM_PARTICLE_TYPES] = {0};
  size_t N[NUM_PARTICLE_TYPES] = {0};
  size_t Ndm;

  /* Open file */
  /* message("Opening file '%s' as IC.", fileName); */
  h_file = H5Fopen(fileName, H5F_ACC_RDONLY, H5P_DEFAULT);
  if (h_file < 0) {
    error("Error while opening file '%s'.", fileName);
  }

  /* Open header to read simulation properties */
  /* message("Reading runtime parameters..."); */
  h_grp = H5Gopen(h_file, "/RuntimePars", H5P_DEFAULT);
  if (h_grp < 0) error("Error while opening runtime parameters\n");

  /* Read the relevant information */
  readAttribute(h_grp, "PeriodicBoundariesOn", INT, periodic);

  /* Close runtime parameters */
  H5Gclose(h_grp);

  /* Open header to read simulation properties */
  /* message("Reading file header..."); */
  h_grp = H5Gopen(h_file, "/Header", H5P_DEFAULT);
  if (h_grp < 0) error("Error while opening file header\n");

  /* Read the relevant information and print status */
  readAttribute(h_grp, "BoxSize", DOUBLE, boxSize);
  readAttribute(h_grp, "NumPart_Total", UINT, numParticles);
  readAttribute(h_grp, "NumPart_Total_HighWord", UINT, numParticles_highWord);

  for (int ptype = 0; ptype < NUM_PARTICLE_TYPES; ++ptype)
    N[ptype] = ((long long)numParticles[ptype]) +
               ((long long)numParticles_highWord[ptype] << 32);

  dim[0] = boxSize[0];
  dim[1] = (boxSize[1] < 0) ? boxSize[0] : boxSize[1];
  dim[2] = (boxSize[2] < 0) ? boxSize[0] : boxSize[2];

  /* message("Found %d particles in a %speriodic box of size [%f %f %f].",  */
  /* 	  *N, (periodic ? "": "non-"), dim[0], dim[1], dim[2]);  */

  /* Close header */
  H5Gclose(h_grp);

  /* Allocate memory to store SPH particles */
  *Ngas = N[0];
  if (posix_memalign((void*)parts, part_align, *Ngas * sizeof(struct part)) !=
      0)
    error("Error while allocating memory for SPH particles");
  bzero(*parts, *Ngas * sizeof(struct part));

  /* Allocate memory to store all particles */
  Ndm = N[1];
  *Ngparts = N[1] + N[0];
  if (posix_memalign((void*)gparts, gpart_align,
                     *Ngparts * sizeof(struct gpart)) != 0)
    error("Error while allocating memory for gravity particles");
  bzero(*gparts, *Ngparts * sizeof(struct gpart));

  /* message("Allocated %8.2f MB for particles.", *N * sizeof(struct part) /
   * (1024.*1024.)); */

  /* message("BoxSize = %lf", dim[0]); */
  /* message("NumPart = [%zd, %zd] Total = %zd", *Ngas, Ndm, *Ngparts); */

  /* Loop over all particle types */
  for (int ptype = 0; ptype < NUM_PARTICLE_TYPES; ptype++) {

    /* Don't do anything if no particle of this kind */
    if (N[ptype] == 0) continue;

    /* Open the particle group in the file */
    char partTypeGroupName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(partTypeGroupName, PARTICLE_GROUP_BUFFER_SIZE, "/PartType%d",
             ptype);
    h_grp = H5Gopen(h_file, partTypeGroupName, H5P_DEFAULT);
    if (h_grp < 0) {
      error("Error while opening particle group %s.", partTypeGroupName);
    }

    /* message("Group %s found - reading...", partTypeGroupName); */

    /* Read particle fields into the particle structure */
    switch (ptype) {

      case GAS:
        hydro_read_particles(h_grp, *Ngas, *Ngas, 0, *parts);
        break;

      case DM:
        darkmatter_read_particles(h_grp, Ndm, Ndm, 0, *gparts);
        break;

      default:
        error("Particle Type %d not yet supported. Aborting", ptype);
    }

    /* Close particle group */
    H5Gclose(h_grp);
  }

  /* Prepare the DM particles */
  prepare_dm_gparts(*gparts, Ndm);

  /* Now duplicate the hydro particle into gparts */
  duplicate_hydro_gparts(*parts, *gparts, *Ngas, Ndm);

  /* message("Done Reading particles..."); */

  /* Close file */
  H5Fclose(h_file);
}

/**
 * @brief Writes an HDF5 output file (GADGET-3 type) with its XMF descriptor
 *
 * @param e The engine containing all the system.
 * @param us The UnitSystem used for the conversion of units in the output
 *
 * Creates an HDF5 output file and writes the particles contained
 * in the engine. If such a file already exists, it is erased and replaced
 * by the new one.
 * The companion XMF file is also updated accordingly.
 *
 * Calls #error() if an error occurs.
 *
 */
void write_output_single(struct engine* e, struct UnitSystem* us) {

  hid_t h_file = 0, h_grp = 0;
  const size_t Ngas = e->s->nr_parts;
  const size_t Ntot = e->s->nr_gparts;
  int periodic = e->s->periodic;
  int numFiles = 1;
  struct part* parts = e->s->parts;
  struct gpart* gparts = e->s->gparts;
  struct gpart* dmparts = NULL;
  static int outputCount = 0;

  /* Number of particles of each type */
  // const size_t Ndm = Ntot - Ngas;

  /* MATTHIEU: Temporary fix to preserve master */
  const size_t Ndm = Ntot > 0 ? Ntot - Ngas : 0;
  /* MATTHIEU: End temporary fix */

  long long N_total[NUM_PARTICLE_TYPES] = {Ngas, Ndm, 0};

  /* File name */
  char fileName[FILENAME_BUFFER_SIZE];
  snprintf(fileName, FILENAME_BUFFER_SIZE, "output_%03i.hdf5", outputCount);

  /* First time, we need to create the XMF file */
  if (outputCount == 0) createXMFfile();

  /* Prepare the XMF file for the new entry */
  FILE* xmfFile = 0;
  xmfFile = prepareXMFfile();

  /* Write the part corresponding to this specific output */
  writeXMFoutputheader(xmfFile, fileName, e->time);

  /* Open file */
  /* message("Opening file '%s'.", fileName); */
  h_file = H5Fcreate(fileName, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  if (h_file < 0) {
    error("Error while opening file '%s'.", fileName);
  }

  /* Open header to write simulation properties */
  /* message("Writing runtime parameters..."); */
  h_grp =
      H5Gcreate(h_file, "/RuntimePars", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating runtime parameters group\n");

  /* Write the relevant information */
  writeAttribute(h_grp, "PeriodicBoundariesOn", INT, &periodic, 1);

  /* Close runtime parameters */
  H5Gclose(h_grp);

  /* Open header to write simulation properties */
  /* message("Writing file header..."); */
  h_grp = H5Gcreate(h_file, "/Header", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating file header\n");

  /* Print the relevant information and print status */
  writeAttribute(h_grp, "BoxSize", DOUBLE, e->s->dim, 3);
  double dblTime = e->time;
  writeAttribute(h_grp, "Time", DOUBLE, &dblTime, 1);

  /* GADGET-2 legacy values */
  /* Number of particles of each type */
  unsigned int numParticles[NUM_PARTICLE_TYPES] = {0};
  unsigned int numParticlesHighWord[NUM_PARTICLE_TYPES] = {0};
  for (int ptype = 0; ptype < NUM_PARTICLE_TYPES; ++ptype) {
    numParticles[ptype] = (unsigned int)N_total[ptype];
    numParticlesHighWord[ptype] = (unsigned int)(N_total[ptype] >> 32);
  }
  writeAttribute(h_grp, "NumPart_ThisFile", LONGLONG, N_total,
                 NUM_PARTICLE_TYPES);
  writeAttribute(h_grp, "NumPart_Total", UINT, numParticles,
                 NUM_PARTICLE_TYPES);
  writeAttribute(h_grp, "NumPart_Total_HighWord", UINT, numParticlesHighWord,
                 NUM_PARTICLE_TYPES);
  double MassTable[NUM_PARTICLE_TYPES] = {0};
  writeAttribute(h_grp, "MassTable", DOUBLE, MassTable, NUM_PARTICLE_TYPES);
  unsigned int flagEntropy[NUM_PARTICLE_TYPES] = {0};
  writeAttribute(h_grp, "Flag_Entropy_ICs", UINT, flagEntropy,
                 NUM_PARTICLE_TYPES);
  writeAttribute(h_grp, "NumFilesPerSnapshot", INT, &numFiles, 1);

  /* Close header */
  H5Gclose(h_grp);

  /* Print the code version */
  writeCodeDescription(h_file);

  /* Print the SPH parameters */
  h_grp = H5Gcreate(h_file, "/SPH", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  if (h_grp < 0) error("Error while creating SPH group");
  writeSPHflavour(h_grp);
  H5Gclose(h_grp);

  /* Print the system of Units */
  writeUnitSystem(h_file, us);

  /* Loop over all particle types */
  for (int ptype = 0; ptype < NUM_PARTICLE_TYPES; ptype++) {

    /* Don't do anything if no particle of this kind */
    if (numParticles[ptype] == 0) continue;

    /* Add the global information for that particle type to the XMF meta-file */
    writeXMFgroupheader(xmfFile, fileName, numParticles[ptype], ptype);

    /* Open the particle group in the file */
    char partTypeGroupName[PARTICLE_GROUP_BUFFER_SIZE];
    snprintf(partTypeGroupName, PARTICLE_GROUP_BUFFER_SIZE, "/PartType%d",
             ptype);
    h_grp = H5Gcreate(h_file, partTypeGroupName, H5P_DEFAULT, H5P_DEFAULT,
                      H5P_DEFAULT);
    if (h_grp < 0) {
      error("Error while creating particle group.\n");
    }

    /* message("Writing particle arrays..."); */

    /* Write particle fields from the particle structure */
    switch (ptype) {

      case GAS:
        hydro_write_particles(h_grp, fileName, partTypeGroupName, xmfFile, Ngas,
                              Ngas, 0, 0, parts, us);
        break;

      case DM:
        /* Allocate temporary array */
        if (posix_memalign((void*)&dmparts, gpart_align,
                           Ndm * sizeof(struct gpart)) != 0)
          error("Error while allocating temporart memory for DM particles");
        bzero(dmparts, Ndm * sizeof(struct gpart));

        /* Collect the DM particles from gpart */
        collect_dm_gparts(gparts, Ntot, dmparts, Ndm);

        /* Write DM particles */
        darkmatter_write_particles(h_grp, fileName, partTypeGroupName, xmfFile,
                                   Ndm, Ndm, 0, 0, dmparts, us);

        /* Free temporary array */
        free(dmparts);
        break;

      default:
        error("Particle Type %d not yet supported. Aborting", ptype);
    }

    /* Close particle group */
    H5Gclose(h_grp);

    /* Close this particle group in the XMF file as well */
    writeXMFgroupfooter(xmfFile, ptype);
  }

  /* Write LXMF file descriptor */
  writeXMFoutputfooter(xmfFile, outputCount, e->time);

  /* message("Done writing particles..."); */

  /* Close file */
  H5Fclose(h_file);

  ++outputCount;
}

#endif /* HAVE_HDF5 */
