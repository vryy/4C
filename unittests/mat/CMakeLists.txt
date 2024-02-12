set(TESTNAME unittests_mat)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_druckerprager_test.cpp
    baci_elasthyper_service_test.cpp
    baci_electrode_test.cpp
    baci_inelastic_defgrad_factors_test.cpp
    baci_material_service_test.cpp
    baci_multiplicative_split_defgrad_elasthyper_service_test.cpp
    baci_multiplicative_split_defgrad_elasthyper_test.cpp
    baci_muscle_utils_test.cpp
    baci_stvenantkirchhoff_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_target_requires_cyclic_legacy_linking(${TESTNAME})
