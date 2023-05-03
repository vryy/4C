set(TESTNAME unittests_utils)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_utils_local_numeric_methods.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(
  ${TESTNAME} -Wl,--start-group ${BACI_LIBRARIES} -Wl,--end-group
  )
