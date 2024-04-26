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
