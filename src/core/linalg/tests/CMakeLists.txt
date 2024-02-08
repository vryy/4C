set(TESTNAME unittests_linalg)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_linalg_utils_densematrix_determinant_test.cpp
    baci_linalg_utils_densematrix_eigen_test.cpp
    baci_linalg_utils_densematrix_inverse_test.cpp
    baci_linalg_utils_densematrix_svd_test.cpp
    baci_linalg_utils_fixedsizematrix_test.cpp
    baci_linalg_utils_tensor_transformation_test.cpp
    baci_linalg_utils_voigt_notation_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
