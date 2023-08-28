set(TESTNAME unittests_beam3)

set(SOURCE_LIST
    #cmake - format : sortable
    baci_beam3_euler_bernoulli_test.cpp baci_beam3_kirchhoff_test.cpp baci_beam3_reissner_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME})
baci_target_requires_cyclic_legacy_linking(${TESTNAME})
