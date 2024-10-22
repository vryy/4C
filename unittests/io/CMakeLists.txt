# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(TESTNAME unittests_io)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_discretization_nodal_coordinates_test.cpp
    4C_domainreader_test.cpp
    4C_gridgenerator_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

set(TESTNAME unittests_io_np_3)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_discretization_nodal_coordinates_np_3_test.cpp
    4C_gridgenerator_np_3_test.cpp
    )

four_c_add_google_test_executable(
  ${TESTNAME}
  NP
  3
  SOURCE
  ${SOURCE_LIST}
  )
