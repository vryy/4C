set(TESTNAME unittests_geometry)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_geometry_element_volume_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
