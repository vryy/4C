set(TESTNAME unittests_particle_rigidbody)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_particle_rigidbody_utils_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
