set(TESTNAME unittests_io)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_discretization_nodal_coordinates_test.cpp
    4C_domainreader_test.cpp
    4C_file_reader_test.cpp
    4C_gridgenerator_test.cpp
    4C_pstream_test.cpp
    4C_string_converter_test.cpp
    4C_utils_reader_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

set(TESTNAME unittests_io_np_3)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_discretization_nodal_coordinates_np_3_test.cpp
    4C_gridgenerator_np_3_test.cpp
    )

baci_add_google_test_executable(
  ${TESTNAME}
  NP
  3
  SOURCE
  ${SOURCE_LIST}
  )
