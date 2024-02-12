set(TESTNAME unittests_contact_constitutivelaw)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_brokenrational_contactconstitutivelaw_test.cpp
    baci_cubic_contactconstitutivelaw_test.cpp
    baci_linear_contactconstitutivelaw_test.cpp
    baci_mirco_contactconstitutivelaw_test.cpp
    baci_power_contactconstitutivelaw_test.cpp
    )

baci_add_google_test_executable(
  ${TESTNAME}
  THREADS
  2
  SOURCE
  ${SOURCE_LIST}
  )

baci_target_requires_cyclic_legacy_linking(${TESTNAME})
