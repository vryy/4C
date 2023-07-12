set(TESTNAME unittests_linalg)

set(SOURCE_LIST
    # cmake-format: sortable
    linalg_utils_densematrix_determinant_test.cpp
    linalg_utils_densematrix_eigen_test.cpp
    linalg_utils_densematrix_inverse_test.cpp
    linalg_utils_densematrix_svd_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME} core_linalg)
