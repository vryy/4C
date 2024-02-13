set(TESTNAME unittests_fbi)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_beam_to_fluid_meshtying_pair_gauss_point_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
