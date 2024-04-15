set(TESTNAME unittests_io)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_discretization_nodal_coordinates_test.cpp
    baci_domainreader_test.cpp
    baci_file_reader_test.cpp
    baci_gridgenerator_test.cpp
    baci_pstream_test.cpp
    baci_string_converter_test.cpp
    baci_utils_reader_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

set(TESTNAME unittests_io_np_3)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_discretization_nodal_coordinates_np_3_test.cpp
    baci_gridgenerator_np_3_test.cpp
    )

baci_add_google_test_executable(
  ${TESTNAME}
  NP
  3
  SOURCE
  ${SOURCE_LIST}
  )
