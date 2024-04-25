set(TESTNAME unittests_numerics)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_cubic_spline_interpolation_test.cpp
    4C_utils_local_integration_test.cpp
    4C_utils_local_newton_test.cpp
    4C_utils_local_numeric_methods_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
