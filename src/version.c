/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (C) 2012 Matthieu Schaller (matthieu.schaller@durham.ac.uk).
 * Copyright (C) 2015 Peter W. Draper (p.w.draper@durham.ac.uk).
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

/* Some standard headers. */
#include <stdio.h>

/* This object's header. */
#include "version.h"

/**
 * @brief Return the source code git revision
 *
 * @details The SHA of the code checked out when the library was last built.
 * Will include -dirty if they are local modifications.
 */
const char *git_revision(void) {
  static const char *revision = GIT_REVISION;
  return revision;
}

/**
 * @brief Return the source code git branch
 *
 * @details The name of the current branch when the code was last built.
 */
const char *git_branch(void) {
  static const char *branch = GIT_BRANCH;
  return branch;
}


/**
 * @brief The version of SWIFT
 */
const char *package_version(void) {
  static const char *version = PACKAGE_VERSION;
  return version;
}

/**
 * @brief A description of the package version and code status.
 */
const char *package_description(void) {
  static char buf[256];
  static int initialised = 0;
  if (!initialised) {
    sprintf(buf, "SWIFT version: %s, at revision: %s, branch: %s", 
            PACKAGE_VERSION, GIT_REVISION, GIT_BRANCH);
    initialised = 1;
  }
  return buf;
}

/**
 * @brief Prints a greeting message to the standard output containing code
 * version and revision number
 */
void greetings(void) {

  printf(" Welcome to the cosmological code\n");
  printf("    ______       _________________\n");
  printf("   / ___/ |     / /  _/ ___/_  __/\n");
  printf("   \\__ \\| | /| / // // /_   / /   \n");
  printf("  ___/ /| |/ |/ // // __/  / /    \n");
  printf(" /____/ |__/|__/___/_/    /_/     \n");
  printf(" SPH With Inter-dependent Fine-grained Tasking\n\n");

  printf(" Version : %s\n", package_version());
  printf(" Revision: %s, branch: %s\n", git_revision(), git_branch());
  printf(" Webpage : www.swiftsim.com\n\n");
}
