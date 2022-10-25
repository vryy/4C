###------------------------------------------------------------------ Test definitions

# Determine timeout for each test. Use default one, if it is not passed from the outside.
if(DEFINED ENV{GLOBAL_TEST_TIMEOUT})
  set(GLOBAL_TEST_TIMEOUT $ENV{GLOBAL_TEST_TIMEOUT})
  message("Global test timeout is $ENV{GLOBAL_TEST_TIMEOUT} s (before scaling).")
else()
  # default test timeout, if not passed as an environment variable
  set(GLOBAL_TEST_TIMEOUT 260) # Default timeout

  message(
    "Global test timeout is not passed as an environment variable. It is set to the default ${GLOBAL_TEST_TIMEOUT} s (before scaling)."
    )
endif()

# Determine timeout scale factor. Use default one, if it is not passed from the outside
if(DEFINED ENV{GLOBAL_TEST_TIMEOUT_SCALE})
  set(GLOBAL_TEST_TIMEOUT_SCALE $ENV{GLOBAL_TEST_TIMEOUT_SCALE})
  message("Global test timeout scale is $ENV{GLOBAL_TEST_TIMEOUT_SCALE}.")
else()
  # default test timeout scale, if not passed as an environment variable
  if("${CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")
    set(GLOBAL_TEST_TIMEOUT_SCALE 4) # Default timeout scale for debug configuration
  else()
    set(GLOBAL_TEST_TIMEOUT_SCALE 1) # Default timeout scale
  endif()

  message(
    "Global test timeout scale is not passed as an environment variable. It is set to the default ${GLOBAL_TEST_TIMEOUT_SCALE} for this kind of build."
    )
endif()

# Determine scaled global test timeout
math(EXPR GLOBAL_TEST_TIMEOUT_SCALED "${GLOBAL_TEST_TIMEOUT}*${GLOBAL_TEST_TIMEOUT_SCALE}")
message("The scaled global test timeout is ${GLOBAL_TEST_TIMEOUT_SCALED} s.")

####################################################################
################        Definition of macros       #################
####################################################################

# define this test (name_of_test) as a "setup fixture" (named name_of_fixture). Other tests may
# add a dependency on such a setup fixture through the require_fixture() function. CTest will
# then ensure that a) required fixtures are included if necessary and b) tests are executed in
# the correct order.
#
# In the BACI test suite we need such dependencies between a test e.g. for an initial simulation
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
    set_tests_properties(${name_of_test} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
  else()
    set_tests_properties(${name_of_test} PROPERTIES TIMEOUT ${ARGN})
  endif()
endfunction(set_timeout)

# The macros defined in this section can be used in the file 'TestingFrameworkListOfTests.cmake' to define tests

# DEFAULT BACI TEST - run simulation with .dat file
# Usage in TestingFrameworkListOfTests.cmake: "baci_test(<name_of_input_file> <num_proc> <restart_step> optional: <label>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <num_proc>: number of processors the test should use
# <restart_step>: number of restart step; <""> indicates no restart
# optional: <label>: add a label to the test
macro(baci_test name_of_input_file num_proc restart_step)
  set(name_of_test ${name_of_input_file}-p${num_proc})
  set(test_directory ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}_p${num_proc})
  set(source_file ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat)

  add_test(
    NAME ${name_of_input_file}-p${num_proc}
    COMMAND
      bash -c
      "mkdir -p ${test_directory} && ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS} -np ${num_proc} $<TARGET_FILE:${baciname}> ${source_file} ${test_directory}/xxx"
    )

  require_fixture(${name_of_test} test_cleanup)
  set_processors(${name_of_test} ${num_proc})
  define_setup_fixture(${name_of_test} ${name_of_test})
  set_timeout(${name_of_test})

  if(NOT "${ARGN}" STREQUAL "")
    set_label(${name_of_test} ${ARGN})
  endif()

  if(${restart_step})
    add_test(
      NAME ${name_of_test}-restart
      COMMAND
        bash -c
        "${MPI_RUN} ${MPIEXEC_EXTRA_OPTS} -np ${num_proc} $<TARGET_FILE:${baciname}> ${source_file} ${test_directory}/xxx restart=${restart_step}"
      )

    require_fixture(${name_of_test}-restart "${name_of_test};test_cleanup")
    set_processors(${name_of_test}-restart ${num_proc})
    set_timeout(${name_of_test}-restart)
  endif(${restart_step})
endmacro(baci_test)

# BACI TEST TIMEOUT - run simulation with .dat file and manually defined time for timeout
# Usage in TestingFrameworkListOfTests.cmake: "baci_test_extended_timeout(<name_of_input_file> <num_proc> <restart_step> <testtimeout>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <num_proc>: number of processors the test should use
# <restart_step>: number of restart step; <""> indicates no restart
# <testtimeout>: manually defined duration for test timeout
macro(
  baci_test_extended_timeout
  name_of_input_file
  num_proc
  restart_step
  testtimeout
  )
  set(name_of_test ${name_of_input_file}-p${num_proc})
  set(test_directory ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}_p${num_proc})
  set(source_file ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat)

  # scale testtimeout with the global test timeout scale
  math(EXPR actualtesttimeout "${GLOBAL_TEST_TIMEOUT_SCALE} * ${testtimeout}")

  add_test(
    NAME ${name_of_test}
    COMMAND
      bash -c
      "mkdir -p ${test_directory} && ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS} -np ${num_proc} $<TARGET_FILE:${baciname}> ${source_file} ${test_directory}/xxx"
    )

  require_fixture(${name_of_test} test_cleanup)
  set_processors(${name_of_test} ${num_proc})
  define_setup_fixture(${name_of_test} ${name_of_test})
  set_timeout(${name_of_test} ${actualtesttimeout})

  if(${restart_step})
    add_test(
      NAME ${name_of_test}-restart
      COMMAND
        bash -c
        "${MPI_RUN} ${MPIEXEC_EXTRA_OPTS} -np ${num_proc} $<TARGET_FILE:${baciname}> ${source_file} ${test_directory}/xxx restart=${restart_step}"
      )

    require_fixture(${name_of_test}-restart "${name_of_test};test_cleanup")
    set_processors(${name_of_test}-restart ${num_proc})
    set_timeout(${name_of_test}-restart ${actualtesttimeout})
  endif(${restart_step})
endmacro(baci_test_extended_timeout)

# DEFAULT BACI TEST + POST ENSIGHT - run BACI test and subsequent post ensight test in serial and parallel
# Usage in TestingFrameworkListOfTests.cmake: "baci_test_and_post_ensight_test(<name_of_input_file> <num_proc> <restart_step> optional: <label>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <num_proc>: number of processors the test should use
# <restart_step>: number of restart step; <""> indicates no restart
# optional: <label>: add a label to the test
macro(baci_test_and_post_ensight_test name_of_input_file num_proc restart_step)
  set(test_directory framework_test_output/${name_of_input_file}_p${num_proc})
  set(source_file ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat)

  # run normal testing
  if("${ARGN}" STREQUAL "")
    baci_test(${name_of_input_file} ${num_proc} "${restart_step}")
  else()
    baci_test(${name_of_input_file} ${num_proc} "${restart_step}" ${ARGN})
  endif()

  # additionally run postprocessing in serial mode
  set(RUNPOSTFILTER_SER
      ./post_drt_ensight\ --file=${test_directory}/xxx\ --output=${test_directory}/xxx_SER\ --outputtype=bin\ --stress=ndxyz
      )

  add_test(
    NAME ${name_of_input_file}-p${num_proc}-post_ensight-ser
    COMMAND
      sh -c
      " ${RUNPOSTFILTER_SER} && ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3 ${PROJECT_SOURCE_DIR}/tests/post_processing_test/ensight_comparison.py ${source_file} ${test_directory}/xxx_SER_structure.case"
    )

  require_fixture(
    ${name_of_input_file}-p${num_proc}-post_ensight-ser
    "${name_of_input_file}-p${num_proc};test_cleanup"
    )
  set_processors(${name_of_input_file}-p${num_proc}-post_ensight-ser 1)
  set_timeout(${name_of_input_file}-p${num_proc}-post_ensight-ser)

  # additionally run postprocessing in parallel mode
  set(RUNPOSTFILTER_PAR
      ${MPI_RUN}\ ${MPIEXEC_EXTRA_OPTS}\ -np\ ${num_proc}\ ./post_drt_ensight\ --file=${test_directory}/xxx\ --output=${test_directory}/xxx_PAR\ --outputtype=bin\ --stress=ndxyz
      )

  add_test(
    NAME ${name_of_input_file}-p${num_proc}-post_ensight-par
    COMMAND
      sh -c
      " ${RUNPOSTFILTER_PAR} && ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3 ${PROJECT_SOURCE_DIR}/tests/post_processing_test/ensight_comparison.py ${source_file} ${test_directory}/xxx_PAR_structure.case"
    )

  require_fixture(
    ${name_of_input_file}-p${num_proc}-post_ensight-par
    "${name_of_input_file}-p${num_proc};test_cleanup"
    )
  set_processors(${name_of_input_file}-p${num_proc}-post_ensight-par ${num_proc})
  set_timeout(${name_of_input_file}-p${num_proc}-post_ensight-par)

endmacro(baci_test_and_post_ensight_test)

###########
# RESTART SIMULATION
# CAUTION: This tests bases on results of a previous simulation/test
# Usage in TestingFrameworkListOfTests.cmake: "baci_test_restartonly(<name_of_input_file> <num_proc> <restart_step> <optional: identifier>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <name_of_input_file_restart>: the name of an input file the current restart can base on
# <num_proc>: number of processors the test should use
# <num_proc_base_run>: number of processors of precursor base run
# <restart_step>: number of restart step; <""> indicates no restart
# <optional: identifier>: add an identifier to the file results are read from
macro(
  baci_test_restartonly
  name_of_input_file
  name_of_input_file_restart
  num_proc
  num_proc_base_run
  restart_step
  )
  # set additional output prefix identifier to empty string "" in default claase or to specific string if specified as optional input argument
  if(${ARGC} GREATER 5)
    set(IDENTIFIER ${ARGV5})
  else()
    set(IDENTIFIER "")
  endif()

  set(name_of_test ${name_of_input_file}-p${num_proc})
  set(source_file ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat)

  add_test(
    NAME ${name_of_test}-restart
    COMMAND
      bash -c
      "${MPI_RUN} ${MPIEXEC_EXTRA_OPTS} -np ${num_proc} $<TARGET_FILE:${baciname}> ${source_file} framework_test_output/${name_of_input_file_restart}_p${num_proc_base_run}/xxx${IDENTIFIER} restart=${restart_step}"
    )

  require_fixture(
    ${name_of_test}-restart "${name_of_input_file_restart}-p${num_proc_base_run};test_cleanup"
    )
  set_processors(${name_of_test}-restart ${num_proc})
  set_timeout(${name_of_input_file}-p${num_proc}-restart)

  # Set "RUN_SERIAL TRUE" because result files can only be read by one process.
  set_run_serial(${name_of_test}-restart)
endmacro(baci_test_restartonly)

###########
# NESTED PARALLELISM
# Usage in TestingFrameworkListOfTests.cmake: "baci_test_Nested_Par(<name_of_input_file_1> <name_of_input_file_2> <restart_step>)"
# <name_of_input_file_1>: must equal the name of a .dat file in directory Input for the first test; without ".dat". This test will be executed using 1 process.
# <name_of_input_file_2>: must equal the name of a .dat file in directory Input for the second test; without ".dat". This test will be executed using 2 processes.
# <restart_step>: number of restart step; <""> indicates no restart
macro(baci_test_Nested_Par name_of_input_file_1 name_of_input_file_2 restart_step)
  set(test_directory ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file_1})

  add_test(
    NAME ${name_of_input_file_1}-nestedPar
    COMMAND
      bash -c
      "mkdir -p ${test_directory} &&  ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS} -np 3 $<TARGET_FILE:${baciname}> -ngroup=2 -glayout=1,2 -nptype=separateDatFiles ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file_1}.dat ${test_directory}/xxx ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file_2}.dat ${test_directory}/xxxAdditional"
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
        "${MPI_RUN} ${MPIEXEC_EXTRA_OPTS} -np 3 $<TARGET_FILE:${baciname}> -ngroup=2 -glayout=1,2 -nptype=separateDatFiles ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file_1}.dat ${test_directory}/xxx restart=${restart_step} ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file_2}.dat ${test_directory}/xxxAdditional restart=${restart_step}"
      )

    require_fixture(
      ${name_of_input_file_1}-nestedPar-restart "${name_of_input_file_1}-nestedPar-p3;test_cleanup"
      )
    set_processors(${name_of_input_file_1}-nestedPar-restart 3)
  endif(${restart_step})
endmacro(baci_test_Nested_Par)

###########
# NESTED PARALLELISM WITH COPYDATFILE
# Usage in TestingFrameworkListOfTests.cmake: "baci_test_Nested_Par_CopyDat(<name_of_input_file> <num_proc> <num_groups> optional: <label>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <num_proc>: number of processors the test should use
# <num_groups>: the number of groups
# optional: <label>: add a label to the test
macro(baci_test_Nested_Par_CopyDat name_of_input_file num_proc num_groups)
  set(test_directory framework_test_output/${name_of_input_file}_p${num_proc})

  add_test(
    NAME ${name_of_input_file}-nestedPar_CopyDat-p${num_proc}
    COMMAND
      bash -c
      "mkdir -p ${PROJECT_BINARY_DIR}/${test_directory} &&  ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS} -np ${num_proc} $<TARGET_FILE:${baciname}> -ngroup=${num_groups} -nptype=copyDatFile ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat ${test_directory}/xxx"
    )

  require_fixture(${name_of_input_file}-nestedPar_CopyDat-p${num_proc} test_cleanup)
  set_processors(${name_of_input_file}-nestedPar_CopyDat-p${num_proc} ${num_proc})
  set_timeout(${name_of_input_file}-nestedPar_CopyDat-p${num_proc})
endmacro(baci_test_Nested_Par_CopyDat)

###########
# NESTED PARALLELISM WITH COPYDATFILE AND PRECURSOR SIMULATION
# Usage in TestingFrameworkListOfTests.cmake: "baci_test_Nested_Par_CopyDat_prepost(<name_of_input_file_precursor> <name_of_input_file> <name_of_input_file_post> <num_proc> <num_groups> <restart_step> optional: <label>)"
# <name_of_input_file_precursor>: is the inputfile of a precursor simulation (must equal the name of a .dat file in Input; without ".dat")
# <name_of_input_file>: is an inputfile relying on the output of arg1 (must equal the name of a .dat file in Input; without ".dat")
# <name_of_input_file_post>: is the inputfile of a "postprocessing simulation" restarted from <name_of_input_file> (optional, must equal the name of a .dat file in Input; without ".dat")
# <num_proc>: is the number of procs
# <num_groups>: the number of groups
# <restart_step>: number of restart step for <name_of_input_file_post>
# optional: <label>: add a label to the test
macro(
  baci_test_Nested_Par_CopyDat_prepost
  name_of_input_file_precursor
  name_of_input_file
  name_of_input_file_post
  num_proc
  num_groups
  restart_step
  )
  set(test_directory framework_test_output/${name_of_input_file}_p${num_proc})
  set(source_file ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat)

  # precursor simulation
  add_test(
    NAME ${name_of_input_file}_precursor-p${num_proc}
    COMMAND
      bash -c
      "mkdir -p ${PROJECT_BINARY_DIR}/${test_directory} &&  ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS} -np ${num_proc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file_precursor}.dat ${test_directory}/xxx"
    )

  require_fixture(${name_of_input_file}_precursor-p${num_proc} test_cleanup)
  set_processors(${name_of_input_file}_precursor-p${num_proc} ${num_proc})
  define_setup_fixture(${name_of_input_file}_precursor-p${num_proc} ${name_of_input_file}_precursor)
  set_timeout(${name_of_input_file}_precursor-p${num_proc})

  if(NOT "${ARGN}" STREQUAL "")
    set_label(${name_of_input_file}_precursor-p${num_proc} ${ARGN})
  endif()

  # restart from precursor output
  add_test(
    NAME ${name_of_input_file}-p${num_proc}
    COMMAND
      bash -c
      "${MPI_RUN} ${MPIEXEC_EXTRA_OPTS} -np ${num_proc} $<TARGET_FILE:${baciname}> -ngroup=${num_groups} -nptype=copyDatFile ${source_file} ${test_directory}/xxx"
    )

  require_fixture(${name_of_input_file}-p${num_proc} "${name_of_input_file}_precursor;test_cleanup")
  set_processors(${name_of_input_file}-p${num_proc} ${num_proc})
  define_setup_fixture(${name_of_input_file}-p${num_proc} ${name_of_input_file})
  set_timeout(${name_of_input_file}-p${num_proc})

  if(NOT "${ARGN}" STREQUAL "")
    set_label(${name_of_input_file}-p${num_proc} ${ARGN})
  endif()

  # add postprocessing simulation in case it is required
  if(NOT ${name_of_input_file_post} STREQUAL "")
    add_test(
      NAME ${name_of_input_file}_postprocess-p${num_proc}
      COMMAND
      COMMAND
        bash -c
        "${MPI_RUN} ${MPIEXEC_EXTRA_OPTS} -np ${num_proc} $<TARGET_FILE:${baciname}> -ngroup=${num_groups} -nptype=copyDatFile ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file_post}.dat ${test_directory}/xxx restart=${restart_step}"
      )

    require_fixture(
      ${name_of_input_file}_postprocess-p${num_proc} "${name_of_input_file};test_cleanup"
      )

    set_processors(${name_of_input_file}_postprocess-p${num_proc} ${num_proc})
    set_timeout(${name_of_input_file}_postprocess-p${num_proc})
  endif()
endmacro(baci_test_Nested_Par_CopyDat_prepost)

###########
# FRAMEWORK TESTS - testing the whole framework: Cubit, pre_exodus, BACI, and post-filter
# Usage in TestingFrameworkListOfTests.cmake: "baci_framework_test(<name_of_input_file> <num_proc> <xml_filename>)"
# <name_of_input_file>: must equal the name of a .jou/.bc/.head file in directory tests/framework-test
# <num_proc>: number of processors the test should use
# <xml_filename>: copy any xml-file to the build directory. May also be ""
macro(baci_framework_test name_of_input_file num_proc xml_filename)
  set(name_of_test ${name_of_input_file}-p${num_proc}-fw)
  set(test_directory framework_test_output/${name_of_input_file})

  set(RUNCUBIT
      ${CUBIT_DIR}/cubit\ -batch\ -nographics\ -nojournal\ ${PROJECT_SOURCE_DIR}/tests/framework-test/${name_of_input_file}.jou
      ) # cubit is run to generate an exo file
  set(RUNPREEXODUS
      ./pre_exodus\ --exo=${test_directory}/xxx_${name_of_input_file}.e\ --bc=${PROJECT_SOURCE_DIR}/tests/framework-test/${name_of_input_file}.bc\ --head=${PROJECT_SOURCE_DIR}/tests/framework-test/${name_of_input_file}.head\ --dat=${test_directory}/xxx.dat
      ) # pre_exodus is run to generate a Dat file

  if(NOT ${xml_filename} STREQUAL "")
    file(COPY ${PROJECT_SOURCE_DIR}/Input/${xml_filename} DESTINATION ./${test_directory}/)

    # if a XML file name is given, it is copied from the baci input directory to the build directory
  endif(NOT ${xml_filename} STREQUAL "")

  set(RUNBACI
      ${MPI_RUN}\ ${MPIEXEC_EXTRA_OPTS}\ -np\ ${num_proc}\ $<TARGET_FILE:${baciname}>\ ${test_directory}/xxx.dat\ ${test_directory}/xxx
      ) # baci is run using the generated dat file
  set(RUNPOSTFILTER
      ${MPI_RUN}\ ${MPIEXEC_EXTRA_OPTS}\ -np\ ${num_proc}\ ./post_drt_ensight\ --file=${test_directory}/xxx
      ) # post_drt_ensight is run for the resulting output

  add_test(
    NAME ${name_of_test}
    COMMAND
      bash -c
      "mkdir -p ${PROJECT_BINARY_DIR}/${test_directory} && ${RUNCUBIT} && ${RUNPREEXODUS} && ${RUNBACI} && ${RUNPOSTFILTER}"
    )

  require_fixture(${name_of_test} test_cleanup)
  set_environment(${name_of_test})
  set_fail_expression(${name_of_test})
  set_processors(${name_of_test} ${num_proc})
  set_timeout(${name_of_test})
endmacro(baci_framework_test)

###########
# CUT TESTS
# Usage in TestingFrameworkListOfTests.cmake: "cut_test(<num_proc>)"
# <num_proc>: number of processors the test should use
macro(cut_test num_proc)
  set(name_of_test test-p${num_proc}-cut)
  set(test_directory ${PROJECT_BINARY_DIR}/framework_test_output/cut_test_p${num_proc})

  set(RUNTESTS ${MPI_RUN}\ ${MPIEXEC_EXTRA_OPTS}\ -np\ ${num_proc}\ ${PROJECT_BINARY_DIR}/cut_test)

  add_test(
    NAME ${name_of_test}
    COMMAND bash -c "mkdir -p ${test_directory} && cd ${test_directory} && ${RUNTESTS}"
    )

  require_fixture(${name_of_test} test_cleanup)
  set_fail_expression(${name_of_test})
  set_processors(${name_of_test} ${num_proc})
  set_timeout(${name_of_test})
endmacro(cut_test)

###########
# PREPROCESSING TEST - generate default header file and test pre_exo with it
# Usage in TestingFrameworkListOfTests.cmake: "pre_processing(<name_of_input_file> <num_proc>)"
# <name_of_input_file>: must equal the name of .e/.bc/.head file in tests/pre_processing_test
# <num_proc>: number of processors the test should use
macro(pre_processing name_of_input_file num_proc)
  set(name_of_test ${name_of_input_file}-p${num_proc}-pre_processing)
  set(test_directory ${PROJECT_BINARY_DIR}/framework_test_output/${name_of_input_file}_p${num_proc})

  set(RUNPREEXODUS_NOHEAD
      ./pre_exodus\ --exo=${PROJECT_SOURCE_DIR}/tests/pre_processing_test/${name_of_input_file}.e
      ) # run pre_exodus to generate default head and bc file
  set(RUNPREEXODUS_DEFAULTHEAD
      ./pre_exodus\ --exo=${PROJECT_SOURCE_DIR}/tests/pre_processing_test/${name_of_input_file}.e\ --bc=${PROJECT_SOURCE_DIR}/tests/pre_processing_test/${name_of_input_file}.bc\ --head=default.head\ --dat=${test_directory}/xxx.dat
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
endmacro(pre_processing)

###########
# POSTPROCESSING TEST - run ensight postprocessor on previous test
# CAUTION: This tests bases on results of a previous simulation/test
# Usage in TestingFrameworkListOfTests.cmake: "post_processing(<name_of_input_file> <num_proc> <stresstype> <straintype> <startstep> <optional: identifier> <optional: field>"
# <name_of_input_file>: must equal the name of a .dat file from a previous tests
# <num_proc>: number of processors the test should use
# <num_proc_base_run>: number of processors of precursor base run
# <stresstype>: use post processor with this stresstype
# <straintype>: use post processot with this straintype
# <startstep>: start post processing at this step
# <optional: identifier>: additional identifier that can be added to the test name
# <optional: field>: additional field name that can be added to the test name
macro(
  post_processing
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

  set(name_of_test ${name_of_input_file}${IDENTIFIER}${FIELD}-p${num_proc}-pp)
  # define macros for serial and parallel runs
  set(RUNPOSTFILTER_SER
      ./post_drt_ensight\ --file=${test_directory}/xxx${IDENTIFIER}\ --output=${test_directory}/xxx${IDENTIFIER}_SER_${name_of_input_file}\ --stress=${stresstype}\ --strain=${straintype}\ --start=${startstep}
      )
  set(RUNPOSTFILTER_PAR
      ${MPI_RUN}\ ${MPIEXEC_EXTRA_OPTS}\ -np\ ${num_proc}\ ./post_drt_ensight\ --file=${test_directory}/xxx${IDENTIFIER}\ --output=${test_directory}/xxx${IDENTIFIER}_PAR_${name_of_input_file}\ --stress=${stresstype}\ --strain=${straintype}\ --start=${startstep}
      )

  # specify test case
  add_test(
    NAME ${name_of_input_file}${IDENTIFIER}${FIELD}-p${num_proc}-pp
    COMMAND
      sh -c
      " ${RUNPOSTFILTER_PAR} && ${RUNPOSTFILTER_SER} && ${PVPYTHON} ${PROJECT_SOURCE_DIR}/tests/post_processing_test/comparison.py ${test_directory}/xxx${IDENTIFIER}_PAR_${name_of_input_file}${FIELD}*.case ${test_directory}/xxx${IDENTIFIER}_SER_${name_of_input_file}${FIELD}*.case ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}${IDENTIFIER}${FIELD}.csv ${test_directory}"
    )

  require_fixture(
    ${name_of_input_file}${IDENTIFIER}${FIELD}-p${num_proc}-pp
    "${name_of_input_file}-p${num_proc_base_run};test_cleanup"
    )
  set_environment(${name_of_test})
  set_processors(${name_of_test} ${num_proc})
  set_timeout(${name_of_test})

  # Set "RUN_SERIAL TRUE" because result files can only be read by one process.
  set_run_serial(${name_of_test})
endmacro(post_processing)

###########
# COMPARE ABSOULTE - compare arbitrary result csv file to corresponding reference file with absolute difference
# Implementation can be found in 'utilities/diff_with_tolerance.py'
# CAUTION: This tests bases on results of a previous simulation/test
# Usage in TestingFrameworkListOfTests.cmake: "result_file_abs(<name_of_input_file> <num_proc> <filetag> <resultfilename> <referencefilename> <tolerance>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <num_proc>: number of processors the test should use
# <num_proc_base_run>: number of processors of precursor base run
# <filetag>: add tag to test name
# <resultfilename>: file that should be compared
# <referencefilename>: file to compare with
# <tolerance>: difference the values in <resultfilename> may have to the values in <referencefilename>
macro(
  result_file_abs
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
      ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3
      ${PROJECT_SOURCE_DIR}/utilities/diff_with_tolerance.py ${tolerance}
      ${test_directory}/${resultfilename} ${PROJECT_SOURCE_DIR}/Input/${referencefilename} abs_tol
      0.0
    )

  require_fixture(
    ${name_of_input_file}-p${num_proc}-${filetag}
    "${name_of_input_file}-p${num_proc_base_run};test_cleanup"
    )
  set_processors(${name_of_test} 1)
  set_timeout(${name_of_test})
endmacro(result_file_abs)

###########
# COMPARE RELATIVE - compare arbitrary result csv file to corresponding reference file with relative difference
# Implementation can be found in 'utilities/diff_with_tolerance.py'
# CAUTION: This tests bases on results of a previous simulation/test
# Usage in TestingFrameworkListOfTests.cmake: "result_file_rel(<name_of_input_file> <num_proc> <filetag> <resultfilename> <referencefilename> <tolerance>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <num_proc>: number of processors the test should use
# <num_proc_base_run>: number of processors of precursor base run
# <filetag>: add tag to test name
# <resultfilename>: file that should be compared
# <referencefilename>: file to compare with
# <tolerance>: difference the values in <resultfilename> may have to the values in <referencefilename>
# <min_val>: minimum value of denominator for calculation of relative difference
macro(
  result_file_rel
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
      ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3
      ${PROJECT_SOURCE_DIR}/utilities/diff_with_tolerance.py ${tolerance}
      ${test_directory}/${resultfilename} ${PROJECT_SOURCE_DIR}/Input/${referencefilename} rel_tol
      ${min_val}
    )

  require_fixture(${name_of_test} "${name_of_input_file}-p${num_proc_base_run};test_cleanup")
  set_processors(${name_of_test} 1)
  set_timeout(${name_of_test})
endmacro(result_file_rel)

###########
# COMPARE VTK - compare XML formatted .vtk result data set referenced by .pvd files to corresponding reference files
# CAUTION: This tests bases on results of a previous simulation/test
# Implementation can be found in '/tests/output_test/vtk_compare.py'
# Usage in TestingFrameworkListOfTests.cmake: "vtk_test(<name_of_input_file> <num_proc> <filetag> <pvd_referencefilename> <tolerance> <optional: time_steps>)"
# <name_of_test>: name of this test
# <name_of_input_file>: must equal the name of a .dat file from a previous test
# <num_proc_base_run>: number of processors of precursor base run
# <pvd_referencefilename>: file to compare with
# <tolerance>: difference the values may have
# <optional: time_steps>: time steps when to compare
macro(
  vtk_test
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
    NAME ${name_of_test}-p${num_proc_base_run}
    COMMAND
      ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3
      ${PROJECT_SOURCE_DIR}/tests/output_test/vtk_compare.py ${test_directory}
      ${PROJECT_SOURCE_DIR}/Input/${pvd_referencefilename} ${tolerance} ${num_extra_args}
      ${extra_macro_args}
    )

  require_fixture(
    ${name_of_test}-p${num_proc_base_run}
    "${name_of_input_file}-p${num_proc_base_run};test_cleanup"
    )
  set_processors(${name_of_test}-p${num_proc_base_run} 1)
  set_timeout(${name_of_test}-p${num_proc_base_run})
endmacro(vtk_test)

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
    "if [ -f *_CUTFAIL.pos ]; then mkdir -p ../cut-debug ; cp *_CUTFAIL.pos ../cut-debug/ ; fi ; rm -vfr xxx* framework_test_output* core.* amesos-failure.dat default.bc default.head"
  )
set_processors(test_cleanup 1)
set_tests_properties(test_cleanup PROPERTIES FIXTURES_CLEANUP test_cleanup)
