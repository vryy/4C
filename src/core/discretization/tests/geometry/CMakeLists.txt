set(TESTNAME unittests_geometry1)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_geometry_element_volume_test.cpp
    4C_geometry_reference_configuration_utils_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
