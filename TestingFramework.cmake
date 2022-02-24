###------------------------------------------------------------------ Test definitions

# Determine timeout for each test. Use default one, if it is not passed from the outside.
if (DEFINED ENV{GLOBAL_TEST_TIMEOUT})
  set(GLOBAL_TEST_TIMEOUT $ENV{GLOBAL_TEST_TIMEOUT})
  message("Global test timeout is $ENV{GLOBAL_TEST_TIMEOUT} s (before scaling).")
else ()
  # default test timeout, if not passed as an environment variable
  set(GLOBAL_TEST_TIMEOUT 260) # Default timeout

  message("Global test timeout is not passed as an environment variable. It is set to the default ${GLOBAL_TEST_TIMEOUT} s (before scaling).")
endif ()

# Determine timeout scale factor. Use default one, if it is not passed from the outside
if (DEFINED ENV{GLOBAL_TEST_TIMEOUT_SCALE})
  set(GLOBAL_TEST_TIMEOUT_SCALE $ENV{GLOBAL_TEST_TIMEOUT_SCALE})
  message("Global test timeout scale is $ENV{GLOBAL_TEST_TIMEOUT_SCALE}.")
else ()
  # default test timeout scale, if not passed as an environment variable
  if ("${CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")
    set(GLOBAL_TEST_TIMEOUT_SCALE 4) # Default timeout scale for debug configuration
  else ()
    set(GLOBAL_TEST_TIMEOUT_SCALE 1) # Default timeout scale
  endif ()


  message("Global test timeout scale is not passed as an environment variable. It is set to the default ${GLOBAL_TEST_TIMEOUT_SCALE} for this kind of build.")
endif ()

# Determine scaled global test timeout
math(EXPR GLOBAL_TEST_TIMEOUT_SCALED "${GLOBAL_TEST_TIMEOUT}*${GLOBAL_TEST_TIMEOUT_SCALE}")
message("The scaled global test timeout is ${GLOBAL_TEST_TIMEOUT_SCALED} s.")

####################################################################
################        Definition of macros       #################
####################################################################

# The macros defined in this section can be used in the file 'TestingFrameworkListOfTests.cmake' to define tests

# DEFAULT BACI TEST - run simulation with .dat file
# Usage in TestingFrameworkListOfTests.cmake: "baci_test(<name_of_input_file> <num_proc> <restart_step> <optional: minimal>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <num_proc>: number of processors the test should use
# <restart_step>: number of restart step; <""> indicates no restart
# <optional: minimal>: add "minimal" to add test to list of minimal tests
macro (baci_test name_of_input_file num_proc restart_step)
  add_test(NAME ${name_of_input_file}-p${num_proc}
    COMMAND ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS_LIST} -np ${num_proc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat xxx)
  if( "${ARGN}" STREQUAL "minimal")
    set_tests_properties(${name_of_input_file}-p${num_proc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED} LABELS minimal)
  else ()
    set_tests_properties(${name_of_input_file}-p${num_proc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
  endif ()

  if (${restart_step})
    add_test(NAME ${name_of_input_file}-p${num_proc}-restart
      COMMAND ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS_LIST} -np ${num_proc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat xxx restart=${restart_step})
    set_tests_properties(${name_of_input_file}-p${num_proc}-restart PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
  endif (${restart_step})
endmacro (baci_test)

# BACI TEST TIMEOUT - run simulation with .dat file and manually defined time for timeout
# Usage in TestingFrameworkListOfTests.cmake: "baci_test_extended_timeout(<name_of_input_file> <num_proc> <restart_step> <testtimeout>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <num_proc>: number of processors the test should use
# <restart_step>: number of restart step; <""> indicates no restart
# <testtimeout>: manually defined duration for test timeout
macro (baci_test_extended_timeout name_of_input_file num_proc restart_step testtimeout)
  # scale testtimeout with the global test timeout scale
  math(EXPR actualtesttimeout "${GLOBAL_TEST_TIMEOUT_SCALE} * ${testtimeout}")

  add_test(NAME ${name_of_input_file}-p${num_proc}
    COMMAND ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS_LIST} -np ${num_proc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat xxx)
  set_tests_properties(${name_of_input_file}-p${num_proc} PROPERTIES TIMEOUT ${actualtesttimeout})

  if (${restart_step})
    add_test(NAME ${name_of_input_file}-p${num_proc}-restart
      COMMAND ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS_LIST} -np ${num_proc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat xxx restart=${restart_step})
    set_tests_properties(${name_of_input_file}-p${num_proc}-restart PROPERTIES TIMEOUT ${actualtesttimeout})
  endif (${restart_step})
endmacro (baci_test_extended_timeout)

# DEFAULT BACI TEST + POST ENSIGHT - run BACI test and subsequent post ensight test in serial and parallel
# Usage in TestingFrameworkListOfTests.cmake: "baci_test_and_post_ensight_test(<name_of_input_file> <num_proc> <restart_step> <optional: minimal>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <num_proc>: number of processors the test should use
# <restart_step>: number of restart step; <""> indicates no restart
# <optional: minimal>: add "minimal" to add test to list of minimal tests
macro (baci_test_and_post_ensight_test name_of_input_file num_proc restart_step)
  # run normal testing
  if( "${ARGN}" STREQUAL "minimal")
    baci_test(${name_of_input_file} ${num_proc} "${restart_step}" minimal)
  else ()
    baci_test(${name_of_input_file} ${num_proc} "${restart_step}")
  endif ()

  # additionally run postprocessing in serial mode
  set(RUNPOSTFILTER_SER ./post_drt_ensight\ --file=xxx\ --output=xxx_SER\ --outputtype=bin\ --stress=ndxyz)

  add_test(NAME ${name_of_input_file}-p${num_proc}-post_ensight-ser
    COMMAND sh -c " ${RUNPOSTFILTER_SER} && ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3 ${PROJECT_SOURCE_DIR}/tests/post_processing_test/ensight_comparison.py ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat xxx_SER_structure.case")
  set_tests_properties(${name_of_input_file}-p${num_proc}-post_ensight-ser PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})

  # additionally run postprocessing in parallel mode
  set(RUNPOSTFILTER_PAR ${MPI_RUN}\ ${MPIEXEC_EXTRA_OPTS}\ -np\ ${num_proc}\ ./post_drt_ensight\ --file=xxx\ --output=xxx_PAR\ --outputtype=bin\ --stress=ndxyz)

  add_test(NAME ${name_of_input_file}-p${num_proc}-post_ensight-par
    COMMAND sh -c " ${RUNPOSTFILTER_PAR} && ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3 ${PROJECT_SOURCE_DIR}/tests/post_processing_test/ensight_comparison.py ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat xxx_PAR_structure.case")
  set_tests_properties(${name_of_input_file}-p${num_proc}-post_ensight-par PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
endmacro (baci_test_and_post_ensight_test)

###########
# RESTART SIMULATION
# CAUTION: This tests bases on results of a previous simulation/test
# Usage in TestingFrameworkListOfTests.cmake: "baci_test_restartonly(<name_of_input_file> <num_proc> <restart_step> <optional: identifier>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <num_proc>: number of processors the test should use
# <restart_step>: number of restart step; <""> indicates no restart
# <optional: identifier>: add an identifier to the file results are read from
macro (baci_test_restartonly name_of_input_file num_proc restart_step)
  # set additional output prefix identifier to empty string "" in default claase or to specific string if specified as optional input argument
  if(${ARGC} GREATER 3)
    set(IDENTIFIER ${ARGV3})
  else()
    set(IDENTIFIER "")
  endif()

  add_test(NAME ${name_of_input_file}-p${num_proc}-restart
    COMMAND ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS_LIST} -np ${num_proc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat xxx${IDENTIFIER} restart=${restart_step})

  set_tests_properties(${name_of_input_file}-p${num_proc}-restart PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
endmacro (baci_test_restartonly)

###########
# NESTED PARALLELISM
# Usage in TestingFrameworkListOfTests.cmake: "baci_test_Nested_Par(<name_of_input_file_1> <name_of_input_file_2> <restart_step>)"
# <name_of_input_file_1>: must equal the name of a .dat file in directory Input for the first test; without ".dat"
# <name_of_input_file_2>: must equal the name of a .dat file in directory Input for the second test; without ".dat"
# <restart_step>: number of restart step; <""> indicates no restart
macro (baci_test_Nested_Par name_of_input_file_1 name_of_input_file_2 restart_step)
  add_test(NAME ${name_of_input_file_1}-nestedPar
    COMMAND ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS_LIST} -np 3 $<TARGET_FILE:${baciname}> -ngroup=2 -glayout=1,2 -nptype=separateDatFiles ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file_1}.dat xxx ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file_2}.dat xxxAdditional)

  if (${restart_step})
    add_test(NAME ${name_of_input_file_1}-nestedPar-restart
    COMMAND ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS_LIST} -np 3 $<TARGET_FILE:${baciname}> -ngroup=2 -glayout=1,2 -nptype=separateDatFiles ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file_1}.dat xxx restart=${restart_step} ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file_2}.dat xxxAdditional restart=${restart_step})
  endif (${restart_step})

  set_tests_properties(${name_of_input_file_1}-nestedPar PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
endmacro (baci_test_Nested_Par)

###########
# NESTED PARALLELISM WITH COPYDATFILE
# Usage in TestingFrameworkListOfTests.cmake: "baci_test_Nested_Par_CopyDat(<name_of_input_file> <num_proc> <num_groups> <optional: minimal>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <num_proc>: number of processors the test should use
# <num_groups>: the number of groups
# <optional: minimal>: add "minimal" to add test to list of minimal tests
macro (baci_test_Nested_Par_CopyDat name_of_input_file num_proc num_groups)
  add_test(NAME ${name_of_input_file}-nestedPar_CopyDat-p${num_proc}
    COMMAND ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS_LIST} -np ${num_proc} $<TARGET_FILE:${baciname}> -ngroup=${num_groups} -nptype=copyDatFile ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat xxx )
    if( "${ARGN}" STREQUAL "minimal")
      set_tests_properties(${name_of_input_file}-nestedPar_CopyDat-p${num_proc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED} LABELS minimal)
    else ()
      set_tests_properties(${name_of_input_file}-nestedPar_CopyDat-p${num_proc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
    endif ()
endmacro (baci_test_Nested_Par_CopyDat)

###########
# NESTED PARALLELISM WITH COPYDATFILE AND PRECURSOR SIMULATION
# Usage in TestingFrameworkListOfTests.cmake: "baci_test_Nested_Par_CopyDat_prepost(<name_of_input_file_precursor> <name_of_input_file> <name_of_input_file_post> <num_proc> <num_groups> <restart_step> <optional: minimal>)"
# <name_of_input_file_precursor>: is the inputfile of a precursor simulation (must equal the name of a .dat file in Input; without ".dat")
# <name_of_input_file>: is an inputfile relying on the output of arg1 (must equal the name of a .dat file in Input; without ".dat")
# <name_of_input_file_post>: is the inputfile of a "postprocessing simulation" restarted from <name_of_input_file> (optional, must equal the name of a .dat file in Input; without ".dat")
# <num_proc>: is the number of procs
# <num_groups>: the number of groups
# <restart_step>: number of restart step for <name_of_input_file_post>
# <optional: minimal>: add "minimal" to add test to list of minimal tests
macro (baci_test_Nested_Par_CopyDat_prepost name_of_input_file_precursor name_of_input_file name_of_input_file_post num_proc num_groups restart_step)
  # precursor simulation
  add_test(NAME ${name_of_input_file}_precursor-p${num_proc}
    COMMAND ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS_LIST} -np ${num_proc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file_precursor}.dat xxx)
  # restart from precursor output
  add_test(NAME ${name_of_input_file}-p${num_proc}
    COMMAND ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS_LIST} -np ${num_proc} $<TARGET_FILE:${baciname}> -ngroup=${num_groups} -nptype=copyDatFile ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}.dat xxx)

  # add postprocessing simulation in case
  if (NOT ${name_of_input_file_post} STREQUAL "")
    add_test(NAME ${name_of_input_file}_postprocess-p${num_proc}
      COMMAND ${MPI_RUN} ${MPIEXEC_EXTRA_OPTS_LIST} -np ${num_proc} $<TARGET_FILE:${baciname}> -ngroup=${num_groups} -nptype=copyDatFile ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file_post}.dat xxx restart=${restart_step})
  endif(NOT ${name_of_input_file_post} STREQUAL "")
  if( "${ARGN}" STREQUAL "minimal")
    set_tests_properties(${name_of_input_file}_precursor-p${num_proc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED} LABELS minimal)
    set_tests_properties(${name_of_input_file}-p${num_proc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED} LABELS minimal)
    if (NOT ${name_of_input_file_post} STREQUAL "")
      set_tests_properties(${name_of_input_file}_postprocess-p${num_proc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED} LABELS minimal)
    endif(NOT ${name_of_input_file_post} STREQUAL "")
  else ()
    set_tests_properties(${name_of_input_file}_precursor-p${num_proc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
    set_tests_properties(${name_of_input_file}-p${num_proc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
    if (NOT ${name_of_input_file_post} STREQUAL "")
      set_tests_properties(${name_of_input_file}_postprocess-p${num_proc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
    endif(NOT ${name_of_input_file_post} STREQUAL "")
  endif ()
endmacro (baci_test_Nested_Par_CopyDat_prepost)

###########
# FRAMEWORK TESTS - testing the whole framework: Cubit, pre_exodus, BACI, and post-filter
# Usage in TestingFrameworkListOfTests.cmake: "baci_framework_test(<name_of_input_file> <num_proc> <xml_filename>)"
# <name_of_input_file>: must equal the name of a .jou/.bc/.head file in directory tests/framework-test
# <num_proc>: number of processors the test should use
# <xml_filename>: copy any xml-file to the build directory. May also be ""
macro (baci_framework_test name_of_input_file num_proc xml_filename)
  set (RUNCUBIT ${CUBIT_DIR}/cubit\ -batch\ -nographics\ -nojournal\ ${PROJECT_SOURCE_DIR}/tests/framework-test/${name_of_input_file}.jou) # cubit is run to generate an exo file
  set (RUNPREEXODUS ./pre_exodus\ --exo=xxx_${name_of_input_file}.e\ --bc=${PROJECT_SOURCE_DIR}/tests/framework-test/${name_of_input_file}.bc\ --head=${PROJECT_SOURCE_DIR}/tests/framework-test/${name_of_input_file}.head\ --dat=xxx.dat) # pre_exodus is run to generate a Dat file
  if (NOT ${xml_filename} STREQUAL "")
    file(COPY ${PROJECT_SOURCE_DIR}/Input/${xml_filename} DESTINATION ./) # if a XML file name is given, it is copied from the baci input directory to the build directory
  endif(NOT ${xml_filename} STREQUAL "")
  set (RUNBACI ${MPI_RUN}\ ${MPIEXEC_EXTRA_OPTS}\ -np\ ${num_proc}\ $<TARGET_FILE:${baciname}>\ xxx.dat\ xxx) # baci is run using the generated dat file
  set (RUNPOSTFILTER ${MPI_RUN}\ ${MPIEXEC_EXTRA_OPTS}\ -np\ ${num_proc}\ ./post_drt_ensight\ --file=xxx) # post_drt_ensight is run for the resulting output
  add_test(NAME ${name_of_input_file}-p${num_proc}-fw
  COMMAND sh -c "${RUNCUBIT} && ${RUNPREEXODUS} && ${RUNBACI} && ${RUNPOSTFILTER}")

# note: for the clean-up job in the end, every generated intermediate file has to start with "xxx"

  set_tests_properties(${name_of_input_file}-p${num_proc}-fw PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
  set_tests_properties ( ${name_of_input_file}-p${num_proc}-fw PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR:; ERROR ;Error " )
  set_tests_properties ( ${name_of_input_file}-p${num_proc}-fw PROPERTIES ENVIRONMENT "PATH=$ENV{PATH}" )
endmacro (baci_framework_test)

###########
# CUT TESTS
# Usage in TestingFrameworkListOfTests.cmake: "cut_test(<num_proc>)"
# <num_proc>: number of processors the test should use
macro (cut_test num_proc)
  set (RUNTESTS ${MPI_RUN}\ ${MPIEXEC_EXTRA_OPTS}\ -np\ ${num_proc}\ cut_test)
  add_test(NAME test-p${num_proc}-cut
  COMMAND sh -c "${RUNTESTS}")
  set_tests_properties(test-p${num_proc}-cut PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
  set_tests_properties ( test-p${num_proc}-cut PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR:; ERROR " )
endmacro (cut_test)

###########
# PREPROCESSING TEST - generate default header file and test pre_exo with it
# Usage in TestingFrameworkListOfTests.cmake: "pre_processing(<name_of_input_file> <num_proc>)"
# <name_of_input_file>: must equal the name of .e/.bc/.head file in tests/pre_processing_test
# <num_proc>: number of processors the test should use
macro (pre_processing name_of_input_file num_proc)
  set (RUNPREEXODUS_NOHEAD ./pre_exodus\ --exo=${PROJECT_SOURCE_DIR}/tests/pre_processing_test/${name_of_input_file}.e) # run pre_exodus to generate default head and bc file
  set (RUNPREEXODUS_DEFAULTHEAD ./pre_exodus\ --exo=${PROJECT_SOURCE_DIR}/tests/pre_processing_test/${name_of_input_file}.e\ --bc=${PROJECT_SOURCE_DIR}/tests/pre_processing_test/${name_of_input_file}.bc\ --head=default.head\ --dat=xxx.dat) # run pre_exodus to generate dat file using the default head file
  add_test(NAME ${name_of_input_file}-p${num_proc}-pre_processing
  COMMAND sh -c "${RUNPREEXODUS_NOHEAD} && ${RUNPREEXODUS_DEFAULTHEAD}")
  set_tests_properties(${name_of_input_file}-p${num_proc}-pre_processing PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
  set_tests_properties ( ${name_of_input_file}-p${num_proc}-pre_processing PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR:; ERROR ;Error " )
  set_tests_properties ( ${name_of_input_file}-p${num_proc}-pre_processing PROPERTIES ENVIRONMENT "PATH=$ENV{PATH}" )
endmacro (pre_processing)

###########
# POSTPROCESSING TEST - run ensight postprocessor on previous test
# CAUTION: This tests bases on results of a previous simulation/test
# Usage in TestingFrameworkListOfTests.cmake: "post_processing(<name_of_input_file> <num_proc> <stresstype> <straintype> <startstep> <optional: identifier> <optional: field>"
# <name_of_input_file>: must equal the name of a .dat file from a previous tests
# <num_proc>: number of processors the test should use
# <stresstype>: use post processor with this stresstype
# <straintype>: use post processot with this straintype
# <startstep>: start post processing at this step
# <optional: identifier>: additional identifier that can be added to the test name
# <optional: field>: additional field name that can be added to the test name
macro(post_processing name_of_input_file num_proc stresstype straintype startstep)
  # set additional output prefix identifier to empty string "" in default case or to specific string if specified as optional input argument
  if(${ARGC} GREATER 6)
    set(IDENTIFIER ${ARGV6})
  else()
    set(IDENTIFIER "")
  endif()

  # set field name to empty string "" in default case or to specific string if specified as optional input argument
  if(${ARGC} GREATER 5)
    set(FIELD ${ARGV5})
  else()
    set(FIELD "")
  endif()

  # define macros for serial and parallel runs
  set(RUNPOSTFILTER_SER ./post_drt_ensight\ --file=xxx${IDENTIFIER}\ --output=xxx${IDENTIFIER}_SER_${name_of_input_file}\ --stress=${stresstype}\ --strain=${straintype}\ --start=${startstep})
  set(RUNPOSTFILTER_PAR ${MPI_RUN}\ ${MPIEXEC_EXTRA_OPTS}\ -np\ ${num_proc}\ ./post_drt_ensight\ --file=xxx${IDENTIFIER}\ --output=xxx${IDENTIFIER}_PAR_${name_of_input_file}\ --stress=${stresstype}\ --strain=${straintype}\ --start=${startstep})

  # specify test case
  add_test(NAME ${name_of_input_file}${IDENTIFIER}${FIELD}-p${num_proc}-pp
    COMMAND sh -c " ${RUNPOSTFILTER_PAR} && ${RUNPOSTFILTER_SER} && ${PVPYTHON} ${PROJECT_SOURCE_DIR}/tests/post_processing_test/comparison.py xxx${IDENTIFIER}_PAR_${name_of_input_file}${FIELD}*.case xxx${IDENTIFIER}_SER_${name_of_input_file}${FIELD}*.case ${PROJECT_SOURCE_DIR}/Input/${name_of_input_file}${IDENTIFIER}${FIELD}.csv")
  set_tests_properties(${name_of_input_file}${IDENTIFIER}${FIELD}-p${num_proc}-pp PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
  set_tests_properties(${name_of_input_file}${IDENTIFIER}${FIELD}-p${num_proc}-pp PROPERTIES ENVIRONMENT "PATH=$ENV{PATH}")
endmacro(post_processing)

###########
# COMPARE ABSOULTE - compare arbitrary result csv file to corresponding reference file with absolute difference
# Implementation can be found in 'utilities/diff_with_tolerance.py'
# CAUTION: This tests bases on results of a previous simulation/test
# Usage in TestingFrameworkListOfTests.cmake: "result_file_abs(<name_of_input_file> <num_proc> <filetag> <resultfilename> <referencefilename> <tolerance>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <num_proc>: number of processors the test should use
# <filetag>: add tag to test name
# <resultfilename>: file that should be compared
# <referencefilename>: file to compare with
# <tolerance>: difference the values in <resultfilename> may have to the values in <referencefilename>
macro(result_file_abs name_of_input_file num_proc filetag resultfilename referencefilename tolerance)
  # add test to testing framework
  add_test(NAME ${name_of_input_file}-p${num_proc}-${filetag} COMMAND ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3 ${PROJECT_SOURCE_DIR}/utilities/diff_with_tolerance.py ${tolerance} ${PROJECT_BINARY_DIR}/${resultfilename} ${PROJECT_SOURCE_DIR}/Input/${referencefilename} abs_tol 0.0)

  # set maximum test runtime
  set_tests_properties(${name_of_input_file}-p${num_proc}-${filetag} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
endmacro(result_file_abs)

###########
# COMPARE RELATIVE - compare arbitrary result csv file to corresponding reference file with relative difference
# Implementation can be found in 'utilities/diff_with_tolerance.py'
# CAUTION: This tests bases on results of a previous simulation/test
# Usage in TestingFrameworkListOfTests.cmake: "result_file_rel(<name_of_input_file> <num_proc> <filetag> <resultfilename> <referencefilename> <tolerance>)"
# <name_of_input_file>: must equal the name of a .dat file in directory Input; without ".dat"
# <num_proc>: number of processors the test should use
# <filetag>: add tag to test name
# <resultfilename>: file that should be compared
# <referencefilename>: file to compare with
# <tolerance>: difference the values in <resultfilename> may have to the values in <referencefilename>
# <min_val>: minimum value of denominator for calculation of relative difference
macro(result_file_rel name_of_input_file num_proc filetag resultfilename referencefilename tolerance min_val)
  # add test to testing framework
  add_test(NAME ${name_of_input_file}-p${num_proc}-${filetag} COMMAND ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3 ${PROJECT_SOURCE_DIR}/utilities/diff_with_tolerance.py ${tolerance} ${PROJECT_BINARY_DIR}/${resultfilename} ${PROJECT_SOURCE_DIR}/Input/${referencefilename} rel_tol ${min_val})

  # set maximum test runtime
  set_tests_properties(${name_of_input_file}-p${num_proc}-${filetag} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
endmacro(result_file_rel)

###########
# COMPARE VTK - compare XML formatted .vtk result data set referenced by .pvd files to corresponding reference files
# CAUTION: This tests bases on results of a previous simulation/test
# Implementation can be found in '/tests/output_test/vtk_compare.py'
# Usage in TestingFrameworkListOfTests.cmake: "vtk_test(<name_of_input_file> <num_proc> <filetag> <pvd_referencefilename> <tolerance> <optional: time_steps>)"
# <name_of_input_file>: must equal the name of a .dat file from a previous test
# <num_proc>: number of processors the test should use
# <pvd_referencefilename>: file to compare with
# <tolerance>: difference the values may have
# <optional: time_steps>: time steps when to compare
macro(vtk_test name_of_input_file num_proc pvd_resultfilename pvd_referencefilename tolerance)
# this test takes a list of times as extra arguments to check results at those timesteps
# if no extra arguments are given test checks every timestep
  set (extra_macro_args ${ARGN})
    # Did we get any optional args?
    list(LENGTH extra_macro_args num_extra_args)
    if (${num_extra_args} GREATER 0)
        list(GET extra_macro_args 0 optional_arg)
    endif ()
  # add test to testing framework
  add_test(NAME ${name_of_input_file}-p${num_proc} COMMAND ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3 ${PROJECT_SOURCE_DIR}/tests/output_test/vtk_compare.py ${PROJECT_BINARY_DIR}/${pvd_resultfilename} ${PROJECT_SOURCE_DIR}/Input/${pvd_referencefilename} ${tolerance} ${num_extra_args} ${extra_macro_args})

  # set maximum test runtime
  set_tests_properties(${name_of_input_file}-p${num_proc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
endmacro(vtk_test)

###------------------------------------------------------------------ List of tests
include (TestingFrameworkListOfTests.cmake)


###------------------------------------------------------------------ Final cleanup
# remove any output files from our tests
# autogenerated core files (generated by kernel)
# special files generated by trilinos (such as amesos-failure.dat from ML)
add_test(NAME test_cleanup COMMAND sh -c "if [ -f *_CUTFAIL.pos ]; then mkdir -p ../cut-debug ; cp *_CUTFAIL.pos ../cut-debug/ ; fi ; rm -vfr xxx* core.* amesos-failure.dat default.bc default.head")
