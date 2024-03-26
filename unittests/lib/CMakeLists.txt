set(TESTNAME unittests_lib)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_demangle_test.cpp
    baci_inputreader_test.cpp
    baci_linedefinition_test.cpp
    baci_utils_reference_configuration_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

set(TESTNAME unittests_lib_np_3)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_gid_vector_test.cpp
    )

baci_add_google_test_executable(
  ${TESTNAME}
  NP
  3
  SOURCE
  ${SOURCE_LIST}
  )
