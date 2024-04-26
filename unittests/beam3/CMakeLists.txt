set(TESTNAME unittests_beam3)

set(SOURCE_LIST
    #cmake - format : sortable
    4C_beam3_euler_bernoulli_test.cpp 4C_beam3_kirchhoff_test.cpp 4C_beam3_reissner_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
