set(TESTNAME unittests_drt_fbi)

set(SOURCE_LIST unit_beam_to_fluid_meshtying_pair_gauss_point.cpp)

baci_add_google_test_executable(${TESTNAME} ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(
  ${TESTNAME}
  -Wl,--start-group
  ${BACI_LIBRARIES}
  -Wl,--end-group
  ${LIBRARIES}
  )
