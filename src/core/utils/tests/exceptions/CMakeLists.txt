set(TESTNAME unittests_exceptions)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_utils_exceptions_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
