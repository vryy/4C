set(TESTNAME unittests_numerics)

set(SOURCE_LIST
    # cmake-format: sortable
    utils_local_newton_test.cpp
    utils_local_numeric_methods_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME} core_utils_numerics)
