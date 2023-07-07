set(TESTNAME unittests_mat)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_elasthyper_service.cpp
    unit_electrode.cpp
    unit_inelastic_defgrad_factors.cpp
    unit_material_service.cpp
    unit_multiplicative_split_defgrad_elasthyper.cpp
    unit_multiplicative_split_defgrad_elasthyper_service.cpp
    unit_muscle_utils.cpp
    unit_stvenantkirchhoff.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME})
baci_target_requires_cyclic_legacy_linking(${TESTNAME})
