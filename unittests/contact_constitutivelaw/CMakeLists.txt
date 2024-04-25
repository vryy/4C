set(TESTNAME unittests_contact_constitutivelaw)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_brokenrational_contactconstitutivelaw_test.cpp
    4C_cubic_contactconstitutivelaw_test.cpp
    4C_linear_contactconstitutivelaw_test.cpp
    4C_mirco_force_contactconstitutivelaw_test.cpp
    4C_mirco_pressure_contactconstitutivelaw_test.cpp
    4C_power_contactconstitutivelaw_test.cpp
    )

baci_add_google_test_executable(
  ${TESTNAME}
  THREADS
  2
  SOURCE
  ${SOURCE_LIST}
  )
