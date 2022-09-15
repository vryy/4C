set(TESTNAME unittests_drt_beaminteraction)

set(SOURCE_LIST unit_beaminteraction_beam_to_solid_volume_meshtying_mortar.cpp)

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(
  ${TESTNAME}
  -Wl,--start-group
  ${BACI_LIBRARIES}
  -Wl,--end-group
  ${LIBRARIES}
  )
