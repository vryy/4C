set(TESTNAME unittests_so3)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_so_element_service_test.cpp
    baci_so_hex8_determinant_analysis_test.cpp
    baci_so_hex8_test.cpp
    baci_so_tet4_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
