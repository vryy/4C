set(TESTNAME unittests_so3)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_so_element_service.cpp
    unit_so_hex8.cpp
    unit_so_hex8_determinant_analysis.cpp
    unit_so_tet4.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME})

baci_target_requires_cyclic_legacy_linking(${TESTNAME})
