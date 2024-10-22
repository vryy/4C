# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# Finder for CLN
# Exports cln::cln as an imported target
# Note: The CLN_ROOT variable is automatically considered by the find_ calls below.

find_path(CLN_INCLUDE_DIR cln/cln.h)

find_library(CLN_LIBRARY NAMES cln)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CLN DEFAULT_MSG CLN_LIBRARY CLN_INCLUDE_DIR)

if(CLN_FOUND AND NOT TARGET cln::cln)
  add_library(cln::cln UNKNOWN IMPORTED)
  set_target_properties(
    cln::cln
    PROPERTIES IMPORTED_LOCATION "${CLN_LIBRARY}"
               INTERFACE_INCLUDE_DIRECTORIES "${CLN_INCLUDE_DIR}"
    )
endif()
