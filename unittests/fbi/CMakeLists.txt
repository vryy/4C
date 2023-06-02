set(TESTNAME unittests_fbi)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_beam_to_fluid_meshtying_pair_gauss_point.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME})

baci_target_requires_cyclic_legacy_linking(${TESTNAME})
