set(TESTNAME unittests_mortar_interface_utils)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_mortar_interface_utils_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
