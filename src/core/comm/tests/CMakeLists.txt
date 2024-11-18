# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(TESTNAME unittests_comm)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_comm_mpi_utils_test.cpp
    4C_comm_pack_unpack_test.cpp
    4C_comm_utils_test.cpp
    )

four_c_add_google_test_executable(
  ${TESTNAME}
  NP
  5
  SOURCE
  ${SOURCE_LIST}
  )
