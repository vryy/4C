set(TESTNAME unittests_drt_poromultiphase_scatra)

set(SOURCE_LIST unit_poromultiphase_scatra_function.cpp)

BACI_ADD_GOOGLE_TEST_EXECUTABLE(${TESTNAME} ${SOURCE_LIST})

BACI_LINK_GOOGLE_TEST_NECESSARY_LIBRARIES(${TESTNAME} -Wl,--start-group ${BACI_LIBRARIES} -Wl,--end-group ${LIBRARIES})
