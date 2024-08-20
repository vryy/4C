set(TESTNAME unittests_core_io)

set(SOURCE_LIST
    4C_io_file_reader_test.cpp
    4C_io_inputreader_test.cpp
    4C_io_linedefinition_test.cpp
    4C_io_pstream_test.cpp
    4C_io_string_converter_test.cpp
    )

four_c_add_google_test_executable(
  ${TESTNAME}
  NP
  1
  SOURCE
  ${SOURCE_LIST}
  )
