set(TESTNAME unittests_solid_ele)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_solid_3D_ele_calc_lib_test.cpp
    4C_solid_3D_ele_utils_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
