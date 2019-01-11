set_property(GLOBAL PROPERTY geometry_pair_unit_files) # use this instead of lists, because it operates at a global scope
add_subdirectory ( src )
get_property(geometry_pair_unit_list GLOBAL PROPERTY geometry_pair_unit_files)

add_subdirectory ( stubs )

# add a template to cxxtest, which includes a call to MPI_Init and MPI_Finalize
set( CXXTEST_TESTGEN_ARGS --error-printer --template ${CMAKE_CURRENT_SOURCE_DIR}/runner.tpl )
cxxtest_add_test ( geometry_pair_tests geometry_pair_unittest_all.cpp ${geometry_pair_unit_list} )

target_link_libraries( geometry_pair_tests geometry_pair_unittest_stubs -Wl,--start-group ${BACI_LIBRARIES} -Wl,--end-group geometry_pair_dserror_stub ${LIBRARIES})
target_compile_options( geometry_pair_tests PRIVATE -g -O0 )
target_include_directories( geometry_pair_tests PRIVATE ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )

add_dependencies(full geometry_pair_tests)
set_tests_properties( geometry_pair_tests PROPERTIES TIMEOUT 1000 LABELS minimal)
