set(TESTNAME unittests_particle_engine)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_particle_container_bundle_test.cpp
    4C_particle_container_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
