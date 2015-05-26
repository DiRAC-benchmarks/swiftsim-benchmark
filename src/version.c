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

#include <stdio.h>
#include "version.h"

/**
 * @brief Return the source code git revision
 *
 * @details The SHA of the code checked out when the library was last built. 
 * Will include -dirty if they are local modifications.
 */
const char *git_revision( void )
{
    static const char *revision = GIT_REVISION;
    return revision;
}

/**
 * @brief The version of SWIFT
 */
const char *package_version( void )
{
    static const char *version = PACKAGE_VERSION;
    return version;
}

/**
 * @brief A description of the package version and code status.
 */
const char *package_description( void )
{
    static char buf[256];
    static int initialised = 0;
    if ( ! initialised ) {
        sprintf( buf, "SWIFT version: %s, at revision: %s", 
                 PACKAGE_VERSION, GIT_REVISION );
        initialised = 1;
    }
    return buf;
}