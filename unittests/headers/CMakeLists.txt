set(TESTNAME unittests_headers)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_singleton_owner_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
