set(TESTNAME unittests_inpar)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_inpar_structure_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME})
baci_target_requires_cyclic_legacy_linking(${TESTNAME})
