set(TESTNAME unittests_poromultiphase_scatra)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_poromultiphase_scatra_function_test.cpp
    baci_poromultiphase_scatra_utils_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_target_requires_cyclic_legacy_linking(${TESTNAME})
