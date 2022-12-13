set(TESTNAME unittests_drt_poromultiphase_scatra)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_poromultiphase_scatra_function.cpp
    unit_poromultiphase_scatra_utils.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(
  ${TESTNAME} -Wl,--start-group ${BACI_LIBRARIES} -Wl,--end-group
  )
