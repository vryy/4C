set(TESTNAME unittests_mat)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_druckerprager_test.cpp
    4C_elasthyper_service_test.cpp
    4C_electrode_test.cpp
    4C_inelastic_defgrad_factors_test.cpp
    4C_material_service_test.cpp
    4C_multiplicative_split_defgrad_elasthyper_service_test.cpp
    4C_multiplicative_split_defgrad_elasthyper_test.cpp
    4C_muscle_utils_test.cpp
    4C_stvenantkirchhoff_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
