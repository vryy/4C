set(TESTNAME unittests_drt_particle_rigidbody)

set(SOURCE_LIST unit_particle_rigidbody_utils.cpp)

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME} Trilinos::all_selected_libs)
