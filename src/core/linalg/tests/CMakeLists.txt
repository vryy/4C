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
    4C_linalg_utils_sparse_algebra_math_test.cpp
    4C_linalg_vector_test.cpp
    )

four_c_add_google_test_executable(
  ${TESTNAME}
  NP
  2
  SOURCE
  ${SOURCE_LIST}
  )

# copy test matrices to the unit-test build
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_matrices DESTINATION ${CMAKE_BINARY_DIR}/tests)
