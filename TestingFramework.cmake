###------------------------------------------------------------------ Helper Functions

# define this test (name_of_test) as a "setup fixture" (named name_of_fixture). Other tests may
# add a dependency on such a setup fixture through the require_fixture() function. CTest will
# then ensure that a) required fixtures are included if necessary and b) tests are executed in
# the correct order.
#
# In the 4C test suite we need such dependencies between a test e.g. for an initial simulation
# and a test for a restart, or between an initial simulation and a post-processing test.
function(define_setup_fixture name_of_test name_of_fixture)
  set_tests_properties(${name_of_test} PROPERTIES FIXTURES_SETUP ${name_of_fixture})
endfunction()

# add a required test (name_of_required_test) to this test (name_of_test). The required
# test must be defined through define_setup_fixture() as "setup fixture". For more details on why
# these functions are needed, have a look at the documentation of define_setup_fixture().
function(require_fixture name_of_test name_of_required_test)
  set_tests_properties(${name_of_test} PROPERTIES FIXTURES_REQUIRED "${name_of_required_test}")
endfunction()

function(set_environment name_of_test)
  set_tests_properties(${name_of_test} PROPERTIES ENVIRONMENT "PATH=$ENV{PATH}")
endfunction()

# set fail expressions to this test (name_of_test)
function(set_fail_expression name_of_test)
  set_tests_properties(${name_of_test} PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR:; ERROR ;Error ")
endfunction()

# set label to this test (name_of_test)
function(set_label name_of_test label)
  set_tests_properties(${name_of_test} PROPERTIES LABELS ${label})
endfunction()

# set number of processors (num_proc) to this test (name_of_test)
function(set_processors name_of_test num_proc)
  set_tests_properties(${name_of_test} PROPERTIES PROCESSORS ${num_proc})
endfunction()

# set this test(name_of_test) to run in serial
function(set_run_serial name_of_test)
  set_tests_properties(${name_of_test} PROPERTIES RUN_SERIAL TRUE)
endfunction()

# set timeout to this test (name_of_test). Optional, set timeout value
function(set_timeout name_of_test)
  if("${ARGN}" STREQUAL "")
    set_tests_properties(${name_of_test} PROPERTIES TIMEOUT ${FOUR_C_TEST_GLOBAL_TIMEOUT})
  else()
    set_tests_properties(${name_of_test} PROPERTIES TIMEOUT ${ARGN})
  endif()
endfunction()

# add test with options
function(_add_test_with_options)

  set(options "")
  set(oneValueArgs NAME_OF_TEST ADDITIONAL_FIXTURE NP TIMEOUT)
  set(multiValueArgs TEST_COMMAND LABELS)
  cmake_parse_arguments(
    _parsed
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )

  if(DEFINED _parsed_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "There are unparsed arguments: ${_parsed_UNPARSED_ARGUMENTS}!")
  endif()

  if(NOT DEFINED _parsed_NAME_OF_TEST)
    message(FATAL_ERROR "Name of test is a necessary input argument!")
  endif()

  if(NOT DEFINED _parsed_TEST_COMMAND)
    message(FATAL_ERROR "Test command is a necessary input argument!")
  endif()

  if(NOT DEFINED _parsed_ADDITIONAL_FIXTURE)
    set(_parsed_ADDITIONAL_FIXTURE "")
  endif()

  if(NOT DEFINED _parsed_NP)
    set(_parsed_NP 1)
  endif()

  if(NOT DEFINED _parsed_TIMEOUT)
    set(_parsed_TIMEOUT "")
  endif()

  if(NOT DEFINED _parsed_LABELS)
    set(_parsed_LABELS "")
  endif()

  add_test(NAME ${_parsed_NAME_OF_TEST} COMMAND bash -c "${_parsed_TEST_COMMAND}")

  require_fixture(${_parsed_NAME_OF_TEST} "${_parsed_ADDITIONAL_FIXTURE};test_cleanup")
  set_processors(${_parsed_NAME_OF_TEST} ${_parsed_NP})
  define_setup_fixture(${_parsed_NAME_OF_TEST} ${_parsed_NAME_OF_TEST})
  set_timeout(${_parsed_NAME_OF_TEST} ${_parsed_TIMEOUT})

  if(NOT ${_parsed_LABELS} STREQUAL "")
    set_label(${_parsed_NAME_OF_TEST} ${_parsed_LABELS})
  endif()

endfunction()

###------------------------------------------------------------------ 4C Test
# Run simulation with .dat file
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test(<input_file> optional: NP <> RESTART_STEP <> TIMEOUT <> OMP_THREADS <> POST_ENSIGHT_STRUCTURE <> LABEL <>)"

# TEST_FILE:              must equal the name of a .dat file in directory tests/input_files; without ".dat".
#                         If two files are provided the second input file is restarted based on the results of the first input file.

# optional:
# NP:                     Number of processors the test should use. Fallback to 1 if not specified.
#                         For two input files two NP's are required.
# RESTART_STEP:           Number of restart step; not defined indicates no restart
# TIMEOUT:                Manually defined duration for test timeout; defaults to global timeout if not specified
# OMP_THREADS:            Number of OpenMP threads per proccessor the test should use; defaults to deactivated
# POST_ENSIGHT_STRUCTURE: Test post_ensight options in serial and parallel (for structure simulation only!)
# LABELS:                 Add labels to the test

function(four_c_test)

  set(options "")
  set(oneValueArgs RESTART_STEP TIMEOUT OMP_THREADS POST_ENSIGHT_STRUCTURE)
  set(multiValueArgs TEST_FILE NP LABELS)
  cmake_parse_arguments(
    _parsed
    "${options}"
    "${oneValueArgs}"
    "${multiValueArgs}"
    ${ARGN}
    )

  # validate input arguments
  if(DEFINED _parsed_UNPARSED_ARGUMENTS)
    message(FATAL_ERROR "There are unparsed arguments: ${_parsed_UNPARSED_ARGUMENTS}!")
  endif()

  if(NOT DEFINED _parsed_TEST_FILE)
    message(FATAL_ERROR "Test file is required for test!")
  endif()

  if(NOT DEFINED _parsed_NP)
    set(_parsed_NP 1)
  endif()

  if(NOT DEFINED _parsed_RESTART_STEP)
    set(_parsed_RESTART_STEP "")
  endif()

  if(NOT DEFINED _parsed_TIMEOUT)
    set(_parsed_TIMEOUT "")
  endif()

  if(NOT DEFINED _parsed_OMP_THREADS)
    set(_parsed_OMP_THREADS 0)
  endif()

  if(NOT DEFINED _parsed_POST_ENSIGHT_STRUCTURE)
    set(_parsed_POST_ENSIGHT_STRUCTURE OFF)
  endif()

  if(NOT DEFINED _parsed_LABELS)
    set(_parsed_LABELS "")
  endif()

  list(LENGTH _parsed_TEST_FILE num_TEST_FILE)
  list(LENGTH _parsed_NP num_NP)
  if(num_TEST_FILE GREATER 2 OR NOT num_NP EQUAL num_TEST_FILE)
    message(
      FATAL_ERROR
        "You provided more than two test files or the provided number of processors do not match your provided test files!"
      )
  endif()

  # Assert that NP <= 3
  foreach(_np IN LISTS _parsed_NP)
    if(_np GREATER 3)
      message(FATAL_ERROR "Number of processors must be less than or equal to 3!")
    endif()
  endforeach()

  # check if source files exist
  set(source_file "")
  foreach(string IN LISTS _parsed_TEST_FILE)
    set(file_name "${PROJECT_SOURCE_DIR}/tests/input_files/${string}.dat")
    if(NOT EXISTS ${file_name})
      message(FATAL_ERROR "Test source file ${file_name} does not exist!")
    endif()
    list(APPEND source_file ${file_name})
  endforeach()

  # set base test name and directory
  if(num_TEST_FILE EQUAL 1)
    set(base_test_file ${source_file})
    set(base_NP ${_parsed_NP})
    set(name_of_test ${_parsed_TEST_FILE}-p${_parsed_NP})
    set(test_directory ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_test})
  elseif(num_TEST_FILE EQUAL 2)
    list(GET source_file 0 base_test_file)
    list(GET source_file 1 restart_test_file)
    list(GET _parsed_TEST_FILE 0 base_test)
    list(GET _parsed_TEST_FILE 1 restart_test)
    list(GET _parsed_NP 0 base_NP)
    list(GET _parsed_NP 1 restart_NP)
    set(name_of_test
        ${base_test}-p${base_NP}_for_${restart_test}-p${restart_NP}-restart_step_${_parsed_RESTART_STEP}
        )
    set(test_directory ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_test})
  endif()

  set(test_command
      "mkdir -p ${test_directory} && ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_OPTS_FOR_TESTING} -np ${base_NP} $<TARGET_FILE:${FOUR_C_EXECUTABLE_NAME}> ${base_test_file} ${test_directory}/xxx"
      )

  # Optional timeout
  if(NOT "${_parsed_TIMEOUT}" STREQUAL "")
    # scale testtimeout with the global test timeout scale
    math(EXPR _parsed_TIMEOUT "${FOUR_C_TEST_TIMEOUT_SCALE} * ${_parsed_TIMEOUT}")
  endif()

  # Optional OpenMP threads per processor
  set(total_procs ${base_NP})
  if(${_parsed_OMP_THREADS})
    set(name_of_test ${name_of_test}-OMP${_parsed_OMP_THREADS})
    set(test_command
        "export OMP_NUM_THREADS=${_parsed_OMP_THREADS}; ${test_command}; unset OMP_NUM_THREADS"
        )
    math(EXPR total_procs "${_parsed_NP}*${_parsed_OMP_THREADS}")
  endif()

  _add_test_with_options(
    NAME_OF_TEST
    ${name_of_test}
    TEST_COMMAND
    ${test_command}
    NP
    ${total_procs}
    TIMEOUT
    "${_parsed_TIMEOUT}"
    LABELS
    "${_parsed_LABELS}"
    )

  # set additional fixture for restart or post_ensight
  set(additional_fixture "${name_of_test}")

  # restart option
  if(NOT _parsed_RESTART_STEP STREQUAL "" AND num_TEST_FILE EQUAL 1)
    # restart with same input file
    set(name_of_test "${name_of_test}-restart_from_same_input")
    set(test_command "${test_command} restart=${_parsed_RESTART_STEP}")
    _add_test_with_options(
      NAME_OF_TEST
      ${name_of_test}
      TEST_COMMAND
      ${test_command}
      ADDITIONAL_FIXTURE
      ${additional_fixture}
      NP
      ${total_procs}
      TIMEOUT
      "${_parsed_TIMEOUT}"
      LABELS
      "${_parsed_LABELS}"
      )

  elseif(NOT _parsed_RESTART_STEP STREQUAL "" AND num_TEST_FILE EQUAL 2)
    # restart with different input file
    set(name_of_test
        ${restart_test}-p${restart_NP}_from_${base_test}-p${base_NP}-restart_step_${_parsed_RESTART_STEP}
        )
    set(restart_test_directory ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_test})
    set(test_command
        "mkdir -p ${restart_test_directory} && ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_OPTS_FOR_TESTING} -np ${restart_NP} $<TARGET_FILE:${FOUR_C_EXECUTABLE_NAME}> ${restart_test_file} ${restart_test_directory}/xxx restartfrom=${test_directory}/xxx restart=${_parsed_RESTART_STEP}"
        )

    # Optional OpenMP threads per processor
    set(total_procs ${restart_NP})
    if(${_parsed_OMP_THREADS})
      set(name_of_test ${name_of_test}-OMP${_parsed_OMP_THREADS})
      set(test_command
          "export OMP_NUM_THREADS=${_parsed_OMP_THREADS}; ${test_command}; unset OMP_NUM_THREADS"
          )
      math(EXPR total_procs "${_parsed_NP}*${_parsed_OMP_THREADS}")
    endif()

    _add_test_with_options(
      NAME_OF_TEST
      ${name_of_test}
      TEST_COMMAND
      ${test_command}
      ADDITIONAL_FIXTURE
      ${additional_fixture}
      NP
      ${total_procs}
      TIMEOUT
      "${_parsed_TIMEOUT}"
      LABELS
      "${_parsed_LABELS}"
      )
    set_run_serial(${name_of_test})

  endif()

  # post_ensight_structure test in serial and parallel
  if(${_parsed_POST_ENSIGHT_STRUCTURE})
    # serial run
    set(name_of_ensight_test "${name_of_test}-post_ensight_serial")
    set(ensight_command
        "${FOUR_C_ENABLE_ADDRESS_SANITIZER_TEST_OPTIONS}\ ./post_ensight\ --file=${test_directory}/xxx\ --output=${test_directory}/xxx_serial\ --outputtype=bin\ --stress=ndxyz && ${PROJECT_SOURCE_DIR}/utilities/python-venv/bin/python3 ${PROJECT_SOURCE_DIR}/tests/post_processing_test/ensight_comparison.py ${source_file} ${test_directory}/xxx_serial_structure.case"
        )
    _add_test_with_options(
      NAME_OF_TEST
      ${name_of_ensight_test}
      TEST_COMMAND
      ${ensight_command}
      ADDITIONAL_FIXTURE
      ${additional_fixture}
      TIMEOUT
      "${_parsed_TIMEOUT}"
      LABELS
      "${_parsed_LABELS}"
      )

    # parallel run
    set(name_of_ensight_test "${name_of_test}-post_ensight_parallel")
    set(ensight_command
        "${MPIEXEC_EXECUTABLE}\ ${MPIEXEC_EXTRA_OPTS_FOR_TESTING}\ -np\ ${_parsed_NP}\ ./post_ensight\ --file=${test_directory}/xxx\ --output=${test_directory}/xxx_parallel --outputtype=bin\ --stress=ndxyz && ${PROJECT_SOURCE_DIR}/utilities/python-venv/bin/python3 ${PROJECT_SOURCE_DIR}/tests/post_processing_test/ensight_comparison.py ${source_file} ${test_directory}/xxx_parallel_structure.case"
        )
    _add_test_with_options(
      NAME_OF_TEST
      ${name_of_ensight_test}
      TEST_COMMAND
      ${ensight_command}
      ADDITIONAL_FIXTURE
      ${additional_fixture}
      NP
      _parsed_NP
      TIMEOUT
      "${_parsed_TIMEOUT}"
      LABELS
      "${_parsed_LABELS}"
      )
  endif()

endfunction()

###------------------------------------------------------------------ Nested Parallelism
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test_nested_parallelism(<name_of_input_file_1> <name_of_input_file_2> <restart_step>)"
# <name_of_input_file_1>: must equal the name of a .dat file in directory tests/input_files for the first test; without ".dat". This test will be executed using 1 process.
# <name_of_input_file_2>: must equal the name of a .dat file in directory tests/input_files for the second test; without ".dat". This test will be executed using 2 processes.
# <restart_step>: number of restart step; <""> indicates no restart
function(four_c_test_nested_parallelism name_of_input_file_1 name_of_input_file_2 restart_step)
  set(test_directory ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file_1})

  add_test(
    NAME ${name_of_input_file_1}-nestedPar
    COMMAND
      bash -c
      "mkdir -p ${test_directory} &&  ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_OPTS_FOR_TESTING} -np 3 $<TARGET_FILE:${FOUR_C_EXECUTABLE_NAME}> -ngroup=2 -glayout=1,2 -nptype=separateDatFiles ${PROJECT_SOURCE_DIR}/tests/input_files/${name_of_input_file_1}.dat ${test_directory}/xxx ${PROJECT_SOURCE_DIR}/tests/input_files/${name_of_input_file_2}.dat ${test_directory}/xxxAdditional"
    )

  require_fixture(${name_of_input_file_1}-nestedPar test_cleanup)
  set_processors(${name_of_input_file_1}-nestedPar 3)
  define_setup_fixture(${name_of_input_file_1}-nestedPar ${name_of_input_file_1}-nestedPar-p3)
  set_timeout(${name_of_input_file_1}-nestedPar)

  if(${restart_step})
    add_test(
      NAME ${name_of_input_file_1}-nestedPar-restart
      COMMAND
        bash -c
        "${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_OPTS_FOR_TESTING} -np 3 $<TARGET_FILE:${FOUR_C_EXECUTABLE_NAME}> -ngroup=2 -glayout=1,2 -nptype=separateDatFiles ${PROJECT_SOURCE_DIR}/tests/input_files/${name_of_input_file_1}.dat ${test_directory}/xxx restart=${restart_step} ${PROJECT_SOURCE_DIR}/tests/input_files/${name_of_input_file_2}.dat ${test_directory}/xxxAdditional restart=${restart_step}"
      )

    require_fixture(
      ${name_of_input_file_1}-nestedPar-restart "${name_of_input_file_1}-nestedPar-p3;test_cleanup"
      )
    set_processors(${name_of_input_file_1}-nestedPar-restart 3)
  endif()
endfunction()

###------------------------------------------------------------------ Framework Tests
# Testing the whole framework: pre_exodus, 4C, and post-filter
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test_framework(<name_of_input_file> <num_proc> <xml_filename>)"
# <name_of_input_file>: must equal the name of a .e/.bc/.head file in directory tests/framework-test
# <num_proc>: number of processors the test should use
# <xml_filename>: copy any xml-file to the build directory. May also be ""
function(four_c_test_framework name_of_input_file num_proc xml_filename)
  set(name_of_test ${name_of_input_file}-p${num_proc}-fw)
  set(test_directory framework_test_output/${name_of_input_file})

  set(RUNPREEXODUS
      ${FOUR_C_ENABLE_ADDRESS_SANITIZER_TEST_OPTIONS}\ ./pre_exodus\ --exo=${PROJECT_SOURCE_DIR}/tests/framework-test/${name_of_input_file}.e\ --bc=${PROJECT_SOURCE_DIR}/tests/framework-test/${name_of_input_file}.bc\ --head=${PROJECT_SOURCE_DIR}/tests/framework-test/${name_of_input_file}.head\ --dat=${test_directory}/xxx.dat
      ) # pre_exodus is run to generate a Dat file

  if(NOT ${xml_filename} STREQUAL "")
    # if a XML file name is given, it is copied from the 4C input directory to the build directory
    set(RUNCOPYXML
        "cp ${PROJECT_SOURCE_DIR}/tests/input_files/${xml_filename} ./${test_directory}/"
        )
  else()
    # no-op command to do nothing
    set(RUNCOPYXML :)
  endif()

  set(RUNFOURC
      ${MPIEXEC_EXECUTABLE}\ ${MPIEXEC_EXTRA_OPTS_FOR_TESTING}\ -np\ ${num_proc}\ $<TARGET_FILE:${FOUR_C_EXECUTABLE_NAME}>\ ${test_directory}/xxx.dat\ ${test_directory}/xxx
      ) # 4C is run using the generated dat file
  set(RUNPOSTFILTER
      ${MPIEXEC_EXECUTABLE}\ ${MPIEXEC_EXTRA_OPTS_FOR_TESTING}\ -np\ ${num_proc}\ ./post_ensight\ --file=${test_directory}/xxx
      ) # post_ensight is run for the resulting output

  add_test(
    NAME ${name_of_test}
    COMMAND
      bash -c
      "mkdir -p ${PROJECT_BINARY_DIR}/${test_directory} && ${RUNCOPYXML} && ${RUNPREEXODUS} && ${RUNFOURC} && ${RUNPOSTFILTER}"
    )

  require_fixture(${name_of_test} test_cleanup)
  set_environment(${name_of_test})
  set_fail_expression(${name_of_test})
  set_processors(${name_of_test} ${num_proc})
  set_timeout(${name_of_test})
endfunction()

###------------------------------------------------------------------ Cut Tests
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test_cut_test(<num_proc>)"
# <num_proc>: number of processors the test should use
function(four_c_test_cut_test num_proc)
  set(name_of_test test-p${num_proc}-cut)
  set(test_directory ${PROJECT_BINARY_DIR}/framework_test_output/cut_test_p${num_proc})

  set(RUNTESTS
      # Run all the cuttests with num_proc except from alex53
      ${MPIEXEC_EXECUTABLE}\ ${MPIEXEC_EXTRA_OPTS_FOR_TESTING}\ -np\ ${num_proc}\ ${PROJECT_BINARY_DIR}/cut_test\ --ignore_test=alex53
      # Run alex53 serially
      ${FOUR_C_ENABLE_ADDRESS_SANITIZER_TEST_OPTIONS}\ ${PROJECT_BINARY_DIR}/cut_test\ --test=alex53
      )

  add_test(
    NAME ${name_of_test}
    COMMAND bash -c "mkdir -p ${test_directory} && cd ${test_directory} && ${RUNTESTS}"
    )

  require_fixture(${name_of_test} test_cleanup)
  set_fail_expression(${name_of_test})
  set_processors(${name_of_test} ${num_proc})
  set_timeout(${name_of_test})
endfunction()

###------------------------------------------------------------------ Preprocessing Test
# Generate default header file and test pre_exo with it
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test_pre_processing(<name_of_input_file> <num_proc>)"
# <name_of_input_file>: must equal the name of .e/.bc/.head file in tests/pre_processing_test
# <num_proc>: number of processors the test should use
function(four_c_test_pre_processing name_of_input_file num_proc)
  set(name_of_test ${name_of_input_file}-p${num_proc}-pre_processing)
  set(test_directory ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}-p${num_proc})

  set(RUNPREEXODUS_NOHEAD
      ${FOUR_C_ENABLE_ADDRESS_SANITIZER_TEST_OPTIONS}\ ./pre_exodus\ --exo=${PROJECT_SOURCE_DIR}/tests/pre_processing_test/${name_of_input_file}.e
      ) # run pre_exodus to generate default head and bc file
  set(RUNPREEXODUS_DEFAULTHEAD
      ${FOUR_C_ENABLE_ADDRESS_SANITIZER_TEST_OPTIONS}\ ./pre_exodus\ --exo=${PROJECT_SOURCE_DIR}/tests/pre_processing_test/${name_of_input_file}.e\ --bc=${PROJECT_SOURCE_DIR}/tests/pre_processing_test/${name_of_input_file}.bc\ --head=default.head\ --dat=${test_directory}/xxx.dat
      ) # run pre_exodus to generate dat file using the default head file

  add_test(
    NAME ${name_of_test}
    COMMAND
      bash -c "mkdir -p ${test_directory} && ${RUNPREEXODUS_NOHEAD} && ${RUNPREEXODUS_DEFAULTHEAD}"
    )

  require_fixture(${name_of_test} test_cleanup)
  set_environment(${name_of_test})
  set_fail_expression(${name_of_test})
  set_processors(${name_of_test} ${num_proc})
  set_timeout(${name_of_test})
endfunction()

###------------------------------------------------------------------ Postprocessing Test
# Run ensight postprocessor on previous test
# CAUTION: This tests bases on results of a previous simulation/test
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test_post_processing(<name_of_input_file> <num_proc> <stresstype> <straintype> <startstep> <optional: identifier> <optional: field>"
# <name_of_input_file>: must equal the name of a .dat file from a previous tests
# <num_proc>: number of processors the test should use
# <num_proc_base_run>: number of processors of precursor base run
# <stresstype>: use post processor with this stresstype
# <straintype>: use post processot with this straintype
# <startstep>: start post processing at this step
# <optional: identifier>: additional identifier that can be added to the test name
# <optional: field>: additional field name that can be added to the test name
function(
  four_c_test_post_processing
  name_of_input_file
  num_proc
  num_proc_base_run
  stresstype
  straintype
  startstep
  )
  set(test_directory
      ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}-p${num_proc_base_run}
      )

  # set additional output prefix identifier to empty string "" in default case or to specific string if specified as optional input argument
  if(${ARGC} GREATER 7)
    set(IDENTIFIER ${ARGV7})
  else()
    set(IDENTIFIER "")
  endif()

  # set field name to empty string "" in default case or to specific string if specified as optional input argument
  if(${ARGC} GREATER 6)
    set(FIELD ${ARGV6})
  else()
    set(FIELD "")
  endif()

  set(name_of_test "${name_of_input_file}${IDENTIFIER}${FIELD}-p${num_proc}-pp")
  # define macros for serial and parallel runs
  set(RUNPOSTFILTER_SER
      ${FOUR_C_ENABLE_ADDRESS_SANITIZER_TEST_OPTIONS}\ ./post_ensight\ --file=${test_directory}/xxx${IDENTIFIER}\ --output=${test_directory}/xxx${IDENTIFIER}_SER_${name_of_input_file}\ --stress=${stresstype}\ --strain=${straintype}\ --start=${startstep}
      )
  set(RUNPOSTFILTER_PAR
      ${MPIEXEC_EXECUTABLE}\ ${MPIEXEC_EXTRA_OPTS_FOR_TESTING}\ -np\ ${num_proc}\ ./post_ensight\ --file=${test_directory}/xxx${IDENTIFIER}\ --output=${test_directory}/xxx${IDENTIFIER}_PAR_${name_of_input_file}\ --stress=${stresstype}\ --strain=${straintype}\ --start=${startstep}
      )

  # specify test case
  add_test(
    NAME "${name_of_test}"
    COMMAND
      sh -c
      " ${RUNPOSTFILTER_PAR} && ${RUNPOSTFILTER_SER} && ${FOUR_C_PVPYTHON} ${PROJECT_SOURCE_DIR}/tests/post_processing_test/comparison.py ${test_directory}/xxx${IDENTIFIER}_PAR_${name_of_input_file}${FIELD}*.case ${test_directory}/xxx${IDENTIFIER}_SER_${name_of_input_file}${FIELD}*.case ${PROJECT_SOURCE_DIR}/tests/input_files/${name_of_input_file}${IDENTIFIER}${FIELD}.csv ${test_directory}"
    )

  require_fixture("${name_of_test}" "${name_of_input_file}-p${num_proc_base_run};test_cleanup")
  set_environment(${name_of_test})
  set_processors(${name_of_test} ${num_proc})
  set_timeout(${name_of_test})

  # Set "RUN_SERIAL TRUE" because result files can only be read by one process.
  set_run_serial(${name_of_test})
endfunction()

###------------------------------------------------------------------ Compare absolute
# Compare arbitrary result csv file to corresponding reference file with absolute difference
# Implementation can be found in 'utilities/diff_with_tolerance.py'
# CAUTION: This tests bases on results of a previous simulation/test
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test_result_file_abs(<name_of_input_file> <num_proc> <filetag> <resultfilename> <referencefilename> <tolerance>)"
# <name_of_input_file>: must equal the name of a .dat file in directory tests/input_files; without ".dat"
# <num_proc>: number of processors the test should use
# <num_proc_base_run>: number of processors of precursor base run
# <filetag>: add tag to test name
# <resultfilename>: file that should be compared
# <referencefilename>: file to compare with
# <tolerance>: difference the values in <resultfilename> may have to the values in <referencefilename>
function(
  four_c_test_result_file_abs
  name_of_input_file
  num_proc
  num_proc_base_run
  filetag
  resultfilename
  referencefilename
  tolerance
  )
  set(name_of_test ${name_of_input_file}-p${num_proc}-${filetag})
  set(test_directory
      ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}-p${num_proc_base_run}
      )

  # add test to testing framework
  add_test(
    NAME ${name_of_test}
    COMMAND
      ${PROJECT_SOURCE_DIR}/utilities/python-venv/bin/python3
      ${PROJECT_SOURCE_DIR}/utilities/diff_with_tolerance.py ${tolerance}
      ${test_directory}/${resultfilename}
      ${PROJECT_SOURCE_DIR}/tests/input_files/${referencefilename} abs_tol 0.0
    )

  require_fixture(
    ${name_of_input_file}-p${num_proc}-${filetag}
    "${name_of_input_file}-p${num_proc_base_run};test_cleanup"
    )
  set_processors(${name_of_test} 1)
  set_timeout(${name_of_test})
endfunction()

###------------------------------------------------------------------ Compare Relative
# Compare arbitrary result csv file to corresponding reference file with relative difference
# Implementation can be found in 'utilities/diff_with_tolerance.py'
# CAUTION: This tests bases on results of a previous simulation/test
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test_result_file_rel(<name_of_input_file> <num_proc> <filetag> <resultfilename> <referencefilename> <tolerance>)"
# <name_of_input_file>: must equal the name of a .dat file in directory tests/input_files; without ".dat"
# <num_proc>: number of processors the test should use
# <num_proc_base_run>: number of processors of precursor base run
# <filetag>: add tag to test name
# <resultfilename>: file that should be compared
# <referencefilename>: file to compare with
# <tolerance>: difference the values in <resultfilename> may have to the values in <referencefilename>
# <min_val>: minimum value of denominator for calculation of relative difference
function(
  four_c_test_result_file_rel
  name_of_input_file
  num_proc
  num_proc_base_run
  filetag
  resultfilename
  referencefilename
  tolerance
  min_val
  )
  set(name_of_test ${name_of_input_file}-p${num_proc}-${filetag})
  set(test_directory
      ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}-p${num_proc_base_run}
      )

  # add test to testing framework
  add_test(
    NAME ${name_of_test}
    COMMAND
      ${PROJECT_SOURCE_DIR}/utilities/python-venv/bin/python3
      ${PROJECT_SOURCE_DIR}/utilities/diff_with_tolerance.py ${tolerance}
      ${test_directory}/${resultfilename}
      ${PROJECT_SOURCE_DIR}/tests/input_files/${referencefilename} rel_tol ${min_val}
    )

  require_fixture(${name_of_test} "${name_of_input_file}-p${num_proc_base_run};test_cleanup")
  set_processors(${name_of_test} 1)
  set_timeout(${name_of_test})
endfunction()

###------------------------------------------------------------------ Compare VTK
# Compare XML formatted .vtk result data set referenced by .pvd files to corresponding reference files
# CAUTION: This tests bases on results of a previous simulation/test
# Implementation can be found in '/tests/output_test/vtk_compare.py'
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test_vtk(<name_of_input_file> <num_proc> <filetag> <pvd_referencefilename> <tolerance> <optional: time_steps>)"
# <name_of_test>: name of this test
# <name_of_input_file>: must equal the name of a .dat file from a previous test
# <num_proc_base_run>: number of processors of precursor base run
# <pvd_referencefilename>: file to compare with
# <tolerance>: difference the values may have
# <optional: time_steps>: time steps when to compare
function(
  four_c_test_vtk
  name_of_test
  name_of_input_file
  num_proc_base_run
  pvd_resultfilename
  pvd_referencefilename
  tolerance
  )
  set(test_directory
      ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}-p${num_proc_base_run}/${pvd_resultfilename}
      )

  # this test takes a list of times as extra arguments to check results at those timesteps
  # if no extra arguments are given test checks every timestep
  set(extra_macro_args ${ARGN})

  # Did we get any optional args?
  list(LENGTH extra_macro_args num_extra_args)

  if(${num_extra_args} GREATER 0)
    list(GET extra_macro_args 0 optional_arg)
  endif()

  # add test to testing framework
  add_test(
    NAME "${name_of_test}-p${num_proc_base_run}"
    COMMAND
      ${PROJECT_SOURCE_DIR}/utilities/python-venv/bin/python3
      ${PROJECT_SOURCE_DIR}/tests/output_test/vtk_compare.py ${test_directory}
      ${PROJECT_SOURCE_DIR}/tests/input_files/${pvd_referencefilename} ${tolerance}
      ${num_extra_args} ${extra_macro_args}
    )

  require_fixture(
    ${name_of_test}-p${num_proc_base_run}
    "${name_of_input_file}-p${num_proc_base_run};${name_of_input_file}-p${num_proc_base_run}-restart;test_cleanup"
    )
  set_processors(${name_of_test}-p${num_proc_base_run} 1)
  set_timeout(${name_of_test}-p${num_proc_base_run})
endfunction()

###------------------------------------------------------------------ List of tests
include(TestingFrameworkListOfTests.cmake)

###------------------------------------------------------------------ Final cleanup
# remove any output files from our tests
# autogenerated core files (generated by kernel)
# special files generated by trilinos (such as amesos-failure.dat from ML)
add_test(
  NAME test_cleanup
  COMMAND
    sh -c
    "if [ -f *_CUTFAIL.pos ]; then mkdir -p ../cut-debug ; cp *_CUTFAIL.pos ../cut-debug/ ; fi ; rm -fr xxx* framework_test_output* core.* amesos-failure.dat default.bc default.head"
  )
set_processors(test_cleanup 1)
set_tests_properties(test_cleanup PROPERTIES FIXTURES_CLEANUP test_cleanup)
