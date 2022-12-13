set(TESTNAME unittests_drt_matelast)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_elast_coupanisoexpo_elementfibers.cpp
    unit_elast_coupanisoexpo_gausspointfibers.cpp
    unit_elast_coupanisoexpobase.cpp
    unit_elast_coupanisoexposhear_elementfibers.cpp
    unit_elast_coupanisoexposhear_gausspointfibers.cpp
    unit_elast_isoogden.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(
  ${TESTNAME} -Wl,--start-group ${BACI_LIBRARIES} -Wl,--end-group
  )
