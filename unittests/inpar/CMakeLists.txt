set(TESTNAME unittests_inpar)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_inpar_structure_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
