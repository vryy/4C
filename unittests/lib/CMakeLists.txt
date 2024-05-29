set(TESTNAME unittests_lib)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_demangle_test.cpp
    4C_inputreader_test.cpp
    4C_linedefinition_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

set(TESTNAME unittests_lib_np_3)

four_c_add_google_test_executable(
  ${TESTNAME}
  NP
  3
  SOURCE
  ${SOURCE_LIST}
  )
