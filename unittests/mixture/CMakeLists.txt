set(TESTNAME unittests_mixture)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_remodelfiber.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME})
baci_link_executable_with_all(${TESTNAME})
