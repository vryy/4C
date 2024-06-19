set(TESTNAME unittests_comm)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_comm_broadcast_utils_test.cpp
    4C_comm_utils_test.cpp
    )

four_c_add_google_test_executable(
  ${TESTNAME}
  NP
  5
  SOURCE
  ${SOURCE_LIST}
  )
