# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(TESTNAME unittests_contact_constitutivelaw)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_brokenrational_contactconstitutivelaw_test.cpp
    4C_cubic_contactconstitutivelaw_test.cpp
    4C_linear_contactconstitutivelaw_test.cpp
    4C_mirco_force_contactconstitutivelaw_test.cpp
    4C_mirco_pressure_contactconstitutivelaw_test.cpp
    4C_power_contactconstitutivelaw_test.cpp
    )

four_c_add_google_test_executable(
  ${TESTNAME}
  THREADS
  2
  SOURCE
  ${SOURCE_LIST}
  )
