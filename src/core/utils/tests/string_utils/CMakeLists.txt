set(TESTNAME unittests_string_utils)

set(SOURCE_LIST 4C_utils_string_test.cpp)

four_c_add_google_test_executable(
  ${TESTNAME}
  NP
  1
  SOURCE
  ${SOURCE_LIST}
  )
