set(TESTNAME unittests_geometry)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_geometry_element_volume.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME})
baci_target_requires_cyclic_legacy_linking(${TESTNAME})
