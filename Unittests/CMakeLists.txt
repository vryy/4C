set_property(GLOBAL PROPERTY unit_files) # use this instead of lists, because it operates at a global scope
add_subdirectory ( src )
get_property(unit_list GLOBAL PROPERTY unit_files)

add_subdirectory ( stubs )

# add a template to cxxtest, which includes a call to MPI_Init and MPI_Finalize
set( CXXTEST_TESTGEN_ARGS --error-printer --template ${CMAKE_CURRENT_SOURCE_DIR}/runner.tpl )
cxxtest_add_test ( unittest unittest_all.cpp ${unit_list} )

target_link_libraries( unittest unittest_stubs -Wl,--start-group ${BACI_LIBRARIES} -Wl,--end-group dserror_stub ${LIBRARIES})
target_compile_options( unittest PRIVATE -g -O0 )
target_include_directories( unittest PRIVATE ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )

add_dependencies(full unittest)
set_tests_properties( unittest PROPERTIES TIMEOUT 10 LABELS minimal)
