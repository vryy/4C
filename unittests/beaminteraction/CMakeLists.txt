set(TESTNAME unittests_beaminteraction)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_beaminteraction_beam_to_solid_volume_meshtying_mortar_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
