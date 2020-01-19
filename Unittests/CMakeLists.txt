# Determine timeout for the unit tests. Use default one, if it is not passed from the outside.
if (DEFINED ENV{UNITTEST_TIMEOUT})
  set(UNITTEST_TIMEOUT $ENV{UNITTEST_TIMEOUT})
  message("The unit test timeout is $ENV{UNITTEST_TIMEOUT} s.")
else ()
  # default unit test timeout, if not passed as an environment variable
  set(UNITTEST_TIMEOUT 10) # Default unit test timeout

  message("The unit test timeout is not passed as an environment variable. It is set to the default value ${UNITTEST_TIMEOUT} s.")
endif ()

set_property(GLOBAL PROPERTY unit_files) # use this instead of lists, because it operates at a global scope
add_subdirectory ( src )
get_property(unit_list GLOBAL PROPERTY unit_files)

add_subdirectory ( stubs )

# add a template to cxxtest, which includes a call to MPI_Init and MPI_Finalize
set( CXXTEST_TESTGEN_ARGS --error-printer --template ${CMAKE_CURRENT_SOURCE_DIR}/runner.tpl )
cxxtest_add_test ( unittests unittest_all.cpp ${unit_list} )

target_link_libraries( unittests unittest_stubs -Wl,--start-group ${BACI_LIBRARIES} -Wl,--end-group dserror_stub ${LIBRARIES})
target_compile_options( unittests PRIVATE -g -O0 )
target_include_directories( unittests PRIVATE ${PROJECT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR} )

add_dependencies(full unittests)
set_tests_properties( unittests PROPERTIES TIMEOUT ${UNITTEST_TIMEOUT} LABELS minimal)
