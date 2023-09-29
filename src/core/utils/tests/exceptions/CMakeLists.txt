set(TESTNAME unittests_exceptions)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_utils_exceptions_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME} core_utils_exceptions)
