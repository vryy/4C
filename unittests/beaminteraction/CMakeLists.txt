set(TESTNAME unittests_beaminteraction)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_beaminteraction_beam_to_solid_volume_meshtying_mortar_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
