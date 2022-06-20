set(TESTNAME unittests_drt_inpar)

set(SOURCE_LIST unit_inpar_structure.cpp)

baci_add_google_test_executable(${TESTNAME} ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME} drt_lib_inpar)
