set(TESTNAME unittests_headers)

set(SOURCE_LIST unit_singleton_owner.cpp)

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME})
