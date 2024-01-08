set(TESTNAME unittests_io)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_csv_reader_test.cpp
    baci_pstream_test.cpp
    baci_string_converter_test.cpp
    baci_utils_reader_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME} io)
