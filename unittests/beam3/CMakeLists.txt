set(TESTNAME unittests_beam3)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_beam3_euler_bernoulli.cpp
    unit_beam3_kirchhoff.cpp
    unit_beam3_reissner.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME})
baci_link_executable_with_all(${TESTNAME})
