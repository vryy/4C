set(TESTNAME unittests_mortar_interface_utils)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_mortar_interface_utils_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
