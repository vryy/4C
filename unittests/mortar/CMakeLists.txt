set(TESTNAME unittests_mortar_interface_utils)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_mortar_interface_utils.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME} mortar)
