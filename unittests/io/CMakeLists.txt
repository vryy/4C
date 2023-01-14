set(TESTNAME unittests_io)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_csv_reader.cpp
    unit_pstream.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME} lib_io)
