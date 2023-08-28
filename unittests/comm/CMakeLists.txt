set(TESTNAME unittests_comm)

set(SOURCE_LIST baci_comm_utils_test.cpp)

baci_add_google_test_executable(
  ${TESTNAME}
  NP
  5
  SOURCE
  ${SOURCE_LIST}
  )

baci_link_google_test_necessary_libraries(${TESTNAME} comm io)
