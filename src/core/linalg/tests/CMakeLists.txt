# This file is part of 4C multiphysics licensed under the
# GNU Lesser General Public License v3.0 or later.
#
# See the LICENSE.md file in the top-level for license information.
#
# SPDX-License-Identifier: LGPL-3.0-or-later

set(TESTNAME unittests_linalg)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_linalg_utils_densematrix_determinant_test.cpp
    4C_linalg_utils_densematrix_eigen_test.cpp
    4C_linalg_utils_densematrix_inverse_test.cpp
    4C_linalg_utils_densematrix_svd_test.cpp
    4C_linalg_utils_fixedsizematrix_test.cpp
    4C_linalg_utils_tensor_transformation_test.cpp
    4C_linalg_utils_voigt_notation_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

set(TESTNAME unittests_linalg_parallel)
set(SOURCE_LIST
    # cmake-format: sortable
    4C_linalg_utils_sparse_algebra_manipulation_test.cpp
    4C_linalg_utils_sparse_algebra_math_test.cpp
    4C_linalg_vector_test.cpp
    )

file(GLOB_RECURSE SUPPORT_FILES test_matrices/*.mm)

four_c_add_google_test_executable(
  ${TESTNAME}
  NP
  2
  SOURCE
  ${SOURCE_LIST}
  SUPPORT_FILES
  ${SUPPORT_FILES}
  )
