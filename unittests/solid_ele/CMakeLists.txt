set(TESTNAME unittests_solid_ele)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_solid_ele_calc_lib_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_target_requires_cyclic_legacy_linking(${TESTNAME})
