# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

# Finder for FFTW
# Exports fftw::fftw as an imported target
# Note: The FFTW_ROOT variable is automatically considered by the find_ calls below.

find_path(FFTW_INCLUDE_DIR fftw3.h)

find_library(FFTW_LIBRARY fftw3)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW DEFAULT_MSG FFTW_LIBRARY FFTW_INCLUDE_DIR)

if(FFTW_FOUND AND NOT TARGET fftw::fftw)
  add_library(fftw::fftw UNKNOWN IMPORTED)
  set_target_properties(
    fftw::fftw
    PROPERTIES IMPORTED_LOCATION "${FFTW_LIBRARY}"
               INTERFACE_INCLUDE_DIRECTORIES "${FFTW_INCLUDE_DIR}"
    )
endif()
