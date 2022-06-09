set(TESTNAME unittests_linalg)

set(SOURCE_LIST
    unit_linalg_utils_densematrix_determinant.cpp
    unit_linalg_utils_densematrix_eigen.cpp
    unit_linalg_utils_densematrix_inverse.cpp
    unit_linalg_utils_densematrix_svd.cpp
    )

baci_add_google_test_executable(${TESTNAME} ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME} linalg)
