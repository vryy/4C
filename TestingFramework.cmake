###------------------------------------------------------------------ Test definitions

# define this test (name_of_test) as a "setup fixture" (named name_of_fixture). Other tests may
# add a dependency on such a setup fixture through the require_fixture() function. CTest will
# then ensure that a) required fixtures are included if necessary and b) tests are executed in
# the correct order.
#
# In the 4C test suite we need such dependencies between a test e.g. for an initial simulation
# and a test for a restart, or between an initial simulation and a post-processing test.
function(define_setup_fixture name_of_test name_of_fixture)
  set_tests_properties(${name_of_test} PROPERTIES FIXTURES_SETUP ${name_of_fixture})
endfunction(define_setup_fixture)

# add a required test (name_of_required_test) to this test (name_of_test). The required
# test must be defined through define_setup_fixture() as "setup fixture". For more details on why
# these functions are needed, have a look at the documentation of define_setup_fixture().
function(require_fixture name_of_test name_of_required_test)
  set_tests_properties(${name_of_test} PROPERTIES FIXTURES_REQUIRED "${name_of_required_test}")
endfunction(require_fixture)

function(set_environment name_of_test)
  set_tests_properties(${name_of_test} PROPERTIES ENVIRONMENT "PATH=$ENV{PATH}")
endfunction(set_environment)

# set fail expressions to this test (name_of_test)
function(set_fail_expression name_of_test)
  set_tests_properties(${name_of_test} PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR:; ERROR ;Error ")
endfunction(set_fail_expression)

# set label to this test (name_of_test)
function(set_label name_of_test label)
  set_tests_properties(${name_of_test} PROPERTIES LABELS ${label})
endfunction(set_label)

# set number of processors (num_proc) to this test (name_of_test)
function(set_processors name_of_test num_proc)
  set_tests_properties(${name_of_test} PROPERTIES PROCESSORS ${num_proc})
endfunction(set_processors)

# set this test(name_of_test) to run in serial
function(set_run_serial name_of_test)
  set_tests_properties(${name_of_test} PROPERTIES RUN_SERIAL TRUE)
endfunction(set_run_serial)

# set timeout to this test (name_of_test). Optional, set timeout value
function(set_timeout name_of_test)
  if("${ARGN}" STREQUAL "")
    set_tests_properties(${name_of_test} PROPERTIES TIMEOUT ${FOUR_C_TEST_GLOBAL_TIMEOUT})
  else()
    set_tests_properties(${name_of_test} PROPERTIES TIMEOUT ${ARGN})
  endif()
endfunction(set_timeout)

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

# 4C Test - run simulation with .dat file
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test(<input_file> optional: NP <> RESTART_STEP <> TIMEOUT <> OMP_THREADS <> POST_ENSIGHT_STRUCTURE <> LABEL <>)"

# <input_file>:   must equal the name of a .dat file in directory tests/input_files; without ".dat"

# optional:
# NP:                     Number of processors the test should use; 1 if not specified
# RESTART_STEP:           Number of restart step; 0 or not defined indicates no restart
# TIMEOUT:                Manually defined duration for test timeout; default global timeout if not specified
# OMP_THREADS:            Number of OpenMP threads per proccessor the test should use; default deactivated
# POST_ENSIGHT_STRUCTURE: Test post_ensight options in serial and parallel (for structure simulation only!)
# LABELS:                 Add labels to the test

function(four_c_test input_file)

  set(options "")
  set(oneValueArgs
      NP
      RESTART_STEP
      TIMEOUT
      OMP_THREADS
      POST_ENSIGHT_STRUCTURE
      )
  set(multiValueArgs LABELS)
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

  if(NOT DEFINED _parsed_NP)
    set(_parsed_NP 1)
  endif()

  if(NOT DEFINED _parsed_RESTART_STEP)
    set(_parsed_RESTART_STEP 0)
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

  set(name_of_test ${input_file}-p${_parsed_NP})
  set(source_file ${PROJECT_SOURCE_DIR}/tests/input_files/${input_file}.dat)
  # check if .dat file exists
  if(NOT EXISTS ${source_file})
    message(FATAL_ERROR "Test source file ${source_file} does not exist!")
  endif()

  set(additional_fixture "")
  set(test_directory ${PROJECT_BINARY_DIR}/framework_test_output/${input_file}_p${_parsed_NP})

  set(test_command
      "mkdir -p ${test_directory} && ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_OPTS_FOR_TESTING} -np ${_parsed_NP} $<TARGET_FILE:${FOUR_C_EXECUTABLE_NAME}> ${source_file} ${test_directory}/xxx"
      )

  # Optional timeout
  if(NOT "${_parsed_TIMEOUT}" STREQUAL "")
    # scale testtimeout with the global test timeout scale
    math(EXPR _parsed_TIMEOUT "${FOUR_C_TEST_TIMEOUT_SCALE} * ${_parsed_TIMEOUT}")
  endif()

  # Optional OpenMP threads per processor
  set(total_procs ${_parsed_NP})
  if(${_parsed_OMP_THREADS})
    set(name_of_test ${name_of_test}-t${_parsed_OMP_THREADS})
    set(test_directory ${test_directory}_t${_parsed_OMP_THREADS})
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

  # restart option
  if(${_parsed_RESTART_STEP})
    set(additional_fixture "${name_of_test};${additional_fixture}")
    set(name_of_test "${name_of_test}-restart")
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
  endif()

  # post_ensight test in serial and parallel
  if(${_parsed_POST_ENSIGHT_STRUCTURE})
    set(additional_fixture ${name_of_test})
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

# The macros defined in this section can be used in the file 'TestingFrameworkListOfTests.cmake' to define tests

###########
# RESTART SIMULATION
# CAUTION: This tests bases on results of a previous simulation/test
# Usage in TestingFrameworkListOfTests.cmake: "four_c_restart_only(<name_of_input_file> <num_proc> <restart_step> <optional: identifier>)"
# <name_of_test>: give it a name, typically equal to a test file name
# <name_of_input_file>: must equal the name of a .dat file in directory tests/input_files; without ".dat"
# <name_of_input_file_restart>: the name of an input file the current restart can base on
# <num_proc>: number of processors the test should use
# <num_proc_base_run>: number of processors of precursor base run
# <restart_step>: number of restart step; <""> indicates no restart
# <optional: identifier>: add an identifier to the file results are read from
macro(
  four_c_restart_only_with_name
  name_of_restart_test
  name_of_input_file
  name_of_input_file_restart
  num_proc
  num_proc_base_run
  restart_step
  )
  # set additional output prefix identifier to empty string "" in default case or to specific string if specified as optional input argument
  if(${ARGC} GREATER 6)
    set(IDENTIFIER ${ARGV6})
  else()
    set(IDENTIFIER "")
  endif()

  set(name_of_test ${name_of_restart_test}-p${num_proc}-restart-${restart_step})
  set(source_file ${PROJECT_SOURCE_DIR}/tests/input_files/${name_of_input_file}.dat)
  set(test_directory ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}_p${num_proc})

  add_test(
    NAME ${name_of_test}
    COMMAND
      bash -c
      "mkdir -p ${test_directory} && ${MPIEXEC_EXECUTABLE} ${MPIEXEC_EXTRA_OPTS_FOR_TESTING} -np ${num_proc} $<TARGET_FILE:${FOUR_C_EXECUTABLE_NAME}> ${source_file} framework_test_output/${name_of_input_file}_p${num_proc}/xxx${IDENTIFIER} restartfrom=framework_test_output/${name_of_input_file_restart}_p${num_proc_base_run}/xxx${IDENTIFIER} restart=${restart_step}"
    )

  require_fixture(
    ${name_of_test} "${name_of_input_file_restart}-p${num_proc_base_run};test_cleanup"
    )
  set_processors(${name_of_test} ${num_proc})
  set_timeout(${name_of_test})

  # Set "RUN_SERIAL TRUE" because result files can only be read by one process.
  set_run_serial(${name_of_test})
endmacro(four_c_restart_only_with_name)

# RESTART SIMULATION
# CAUTION: This tests bases on results of a previous simulation/test
# Usage in TestingFrameworkListOfTests.cmake: "four_c_restart_only(<name_of_input_file> <num_proc> <restart_step> <optional: identifier>)"
# <name_of_input_file>: must equal the name of a .dat file in directory tests/input_files; without ".dat"
# <name_of_input_file_restart>: the name of an input file the current restart can base on
# <num_proc>: number of processors the test should use
# <num_proc_base_run>: number of processors of precursor base run
# <restart_step>: number of restart step; <""> indicates no restart
# <optional: identifier>: add an identifier to the file results are read from
macro(
  four_c_restart_only
  name_of_input_file
  name_of_input_file_restart
  num_proc
  num_proc_base_run
  restart_step
  )
  set(name_of_test ${name_of_input_file})
  four_c_restart_only_with_name(
    ${name_of_test}
    ${name_of_input_file}
    ${name_of_input_file_restart}
    ${num_proc}
    ${num_proc_base_run}
    ${restart_step}
    )
endmacro(four_c_restart_only)

###########
# NESTED PARALLELISM
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test_nested_parallelism(<name_of_input_file_1> <name_of_input_file_2> <restart_step>)"
# <name_of_input_file_1>: must equal the name of a .dat file in directory tests/input_files for the first test; without ".dat". This test will be executed using 1 process.
# <name_of_input_file_2>: must equal the name of a .dat file in directory tests/input_files for the second test; without ".dat". This test will be executed using 2 processes.
# <restart_step>: number of restart step; <""> indicates no restart
macro(four_c_test_nested_parallelism name_of_input_file_1 name_of_input_file_2 restart_step)
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
  endif(${restart_step})
endmacro(four_c_test_nested_parallelism)

###########
# FRAMEWORK TESTS - testing the whole framework: pre_exodus, 4C, and post-filter
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test_framework(<name_of_input_file> <num_proc> <xml_filename>)"
# <name_of_input_file>: must equal the name of a .e/.bc/.head file in directory tests/framework-test
# <num_proc>: number of processors the test should use
# <xml_filename>: copy any xml-file to the build directory. May also be ""
macro(four_c_test_framework name_of_input_file num_proc xml_filename)
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
  endif(NOT ${xml_filename} STREQUAL "")

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
endmacro(four_c_test_framework)

###########
# CUT TESTS
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test_cut_test(<num_proc>)"
# <num_proc>: number of processors the test should use
macro(four_c_test_cut_test num_proc)
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
endmacro(four_c_test_cut_test)

###########
# PREPROCESSING TEST - generate default header file and test pre_exo with it
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test_pre_processing(<name_of_input_file> <num_proc>)"
# <name_of_input_file>: must equal the name of .e/.bc/.head file in tests/pre_processing_test
# <num_proc>: number of processors the test should use
macro(four_c_test_pre_processing name_of_input_file num_proc)
  set(name_of_test ${name_of_input_file}-p${num_proc}-pre_processing)
  set(test_directory ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}_p${num_proc})

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
endmacro(four_c_test_pre_processing)

###########
# POSTPROCESSING TEST - run ensight postprocessor on previous test
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
macro(
  four_c_test_post_processing
  name_of_input_file
  num_proc
  num_proc_base_run
  stresstype
  straintype
  startstep
  )
  set(test_directory
      ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}_p${num_proc_base_run}
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
endmacro(four_c_test_post_processing)

###########
# COMPARE ABSOULTE - compare arbitrary result csv file to corresponding reference file with absolute difference
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
macro(
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
      ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}_p${num_proc_base_run}
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
endmacro(four_c_test_result_file_abs)

###########
# COMPARE RELATIVE - compare arbitrary result csv file to corresponding reference file with relative difference
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
macro(
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
      ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}_p${num_proc_base_run}
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
endmacro(four_c_test_result_file_rel)

###########
# COMPARE VTK - compare XML formatted .vtk result data set referenced by .pvd files to corresponding reference files
# CAUTION: This tests bases on results of a previous simulation/test
# Implementation can be found in '/tests/output_test/vtk_compare.py'
# Usage in TestingFrameworkListOfTests.cmake: "four_c_test_vtk(<name_of_input_file> <num_proc> <filetag> <pvd_referencefilename> <tolerance> <optional: time_steps>)"
# <name_of_test>: name of this test
# <name_of_input_file>: must equal the name of a .dat file from a previous test
# <num_proc_base_run>: number of processors of precursor base run
# <pvd_referencefilename>: file to compare with
# <tolerance>: difference the values may have
# <optional: time_steps>: time steps when to compare
macro(
  four_c_test_vtk
  name_of_test
  name_of_input_file
  num_proc_base_run
  pvd_resultfilename
  pvd_referencefilename
  tolerance
  )
  set(test_directory
      ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}_p${num_proc_base_run}/${pvd_resultfilename}
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
endmacro(four_c_test_vtk)

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
