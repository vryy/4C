set(TESTNAME unittests_drt_so3)

set(SOURCE_LIST unit_so_hex8.cpp)

baci_add_google_test_executable(${TESTNAME} ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(
  ${TESTNAME}
  -Wl,--start-group
  ${BACI_LIBRARIES}
  -Wl,--end-group
  ${LIBRARIES}
  )
