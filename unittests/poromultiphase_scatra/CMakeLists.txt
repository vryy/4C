set(TESTNAME unittests_poromultiphase_scatra)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_poromultiphase_scatra_function_test.cpp
    4C_poromultiphase_scatra_utils_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
