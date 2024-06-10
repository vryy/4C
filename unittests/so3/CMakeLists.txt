set(TESTNAME unittests_so3)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_so_contact_element_reference_configuration_test.cpp
    4C_so_element_service_test.cpp
    4C_so_hex8_determinant_analysis_test.cpp
    4C_so_hex8_test.cpp
    4C_so_tet4_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
