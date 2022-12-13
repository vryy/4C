set(TESTNAME unittests_drt_contact_constitutivelaw)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_brokenrational_contactconstitutivelaw.cpp
    unit_cubic_contactconstitutivelaw.cpp
    unit_linear_contactconstitutivelaw.cpp
    unit_power_contactconstitutivelaw.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(
  ${TESTNAME} -Wl,--start-group ${BACI_LIBRARIES} -Wl,--end-group
  )
