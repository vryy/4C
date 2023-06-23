set(TESTNAME unittests_particle_engine)

set(SOURCE_LIST
    # cmake-format: sortable
    unit_particle_container.cpp
    unit_particle_container_bundle.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME} particle_engine)
