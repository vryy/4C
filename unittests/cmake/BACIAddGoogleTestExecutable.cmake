#! Add a unit test executable
#
# This executable may be run in serial or with a given number of processes as the NP argument
# Usage: baci_add_google_test_executable(<name> [NP <number of processes>] SOURCE source1 [source2 ...])
#
function(baci_add_google_test_executable TESTNAME)
  set(options "")
  set(oneValueArgs NP)
  set(multiValueArgs SOURCE)
  cmake_parse_arguments(
    BACI_ADD_GOOGLE_TEST_EXECUTABLE
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )

  if(DEFINED BACI_ADD_GOOGLE_TEST_EXECUTABLE_UNPARSED_ARGUMENTS)
    message(
      SEND_ERROR
        "There are unparsed arguments: ${BACI_ADD_GOOGLE_TEST_EXECUTABLE_UNPARSED_ARGUMENTS}"
      )
  endif()

  if(NOT DEFINED BACI_ADD_GOOGLE_TEST_EXECUTABLE_SOURCE)
    message(SEND_ERROR "Need to specify at least one source file.")
  endif()

  if(NOT DEFINED BACI_ADD_GOOGLE_TEST_EXECUTABLE_NP)
    set(BACI_ADD_GOOGLE_TEST_EXECUTABLE_NP 1)
  endif()

  set(assert_mpi_file
      ${CMAKE_CURRENT_BINARY_DIR}/assert_mpi_${TESTNAME}_${BACI_ADD_GOOGLE_TEST_EXECUTABLE_NP}.cpp
      )
  configure_file(${PROJECT_SOURCE_DIR}/unittests/assert_mpi.cpp.in ${assert_mpi_file})

  add_executable(
    ${TESTNAME}
    ${PROJECT_SOURCE_DIR}/unittests/gtest_main_mpi.cpp
    ${assert_mpi_file}
    ${BACI_ADD_GOOGLE_TEST_EXECUTABLE_SOURCE}
    )
  target_include_directories(${TESTNAME} PRIVATE ${PROJECT_SOURCE_DIR})

  # the first process will write a unit test report
  set(mpi_arguments
      -np 1 $<TARGET_FILE:${TESTNAME}> --gtest_output=xml:unittest_reports/${TESTNAME}_report.xml
      )
  # if there is more than one process, spawn the remaining ones without a report
  if(BACI_ADD_GOOGLE_TEST_EXECUTABLE_NP GREATER "1")
    math(EXPR remaining_procs "${BACI_ADD_GOOGLE_TEST_EXECUTABLE_NP}-1")
    list(
      APPEND
      mpi_arguments
      :
      -np
      ${remaining_procs}
      $<TARGET_FILE:${TESTNAME}>
      )
  endif()

  add_test(NAME ${TESTNAME} COMMAND "mpirun" ${mpi_arguments})
  set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT ${UNITTEST_TIMEOUT} LABELS minimal)
  set_tests_properties(${TESTNAME} PROPERTIES PROCESSORS ${BACI_ADD_GOOGLE_TEST_EXECUTABLE_NP})

  add_dependencies(unittests ${TESTNAME})
endfunction()
