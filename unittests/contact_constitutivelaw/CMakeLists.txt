set(TESTNAME unittests_contact_constitutivelaw)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_brokenrational_contactconstitutivelaw.cpp
    unit_cubic_contactconstitutivelaw.cpp
    unit_linear_contactconstitutivelaw.cpp
    unit_mirco_contactconstitutivelaw.cpp
    unit_power_contactconstitutivelaw.cpp
    )

baci_add_google_test_executable(
  ${TESTNAME}
  THREADS
  2
  SOURCE
  ${SOURCE_LIST}
  )

baci_link_google_test_necessary_libraries(${TESTNAME})
baci_link_executable_with_all(${TESTNAME})
