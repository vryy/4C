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

# Definition of default baci tests with restart (and optional "minimal" flag)
macro (baci_test arg nproc restart)
  add_test(NAME ${arg}-p${nproc}
    COMMAND ${MPI_RUN} -np ${nproc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${arg}.dat xxx)
  if( "${ARGN}" STREQUAL "minimal")
    set_tests_properties(${arg}-p${nproc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED} LABELS minimal)
  else ()
    set_tests_properties(${arg}-p${nproc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
  endif ()

  if (${restart})
    add_test(NAME ${arg}-p${nproc}-restart
      COMMAND ${MPI_RUN} -np ${nproc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${arg}.dat xxx restart=${restart})
    set_tests_properties(${arg}-p${nproc}-restart PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
  endif (${restart})
endmacro (baci_test)

# Definition of baci tests with restart and extended runtime
macro (baci_test_extended_timeout arg nproc restart testtimeout)
  # scale testtimeout with the global test timeout scale
  math(EXPR actualtesttimeout "${GLOBAL_TEST_TIMEOUT_SCALE} * ${testtimeout}")

  add_test(NAME ${arg}-p${nproc}
    COMMAND ${MPI_RUN} -np ${nproc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${arg}.dat xxx)
  set_tests_properties(${arg}-p${nproc} PROPERTIES TIMEOUT ${actualtesttimeout})

  if (${restart})
    add_test(NAME ${arg}-p${nproc}-restart
      COMMAND ${MPI_RUN} -np ${nproc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${arg}.dat xxx restart=${restart})
    set_tests_properties(${arg}-p${nproc}-restart PROPERTIES TIMEOUT ${actualtesttimeout})
  endif (${restart})
endmacro (baci_test_extended_timeout)

# Definition of default baci tests with restart and parallel and serial ensight postprocessing
macro (baci_test_and_post_ensight_test arg nproc restart)

  # run normal testing
  if( "${ARGN}" STREQUAL "minimal")
    baci_test(${arg} ${nproc} "${restart}" minimal)
  else ()
    baci_test(${arg} ${nproc} "${restart}")
  endif ()


  # additionally run postprocessing in serial mode
  set(RUNPOSTFILTER_SER ./post_drt_ensight\ --file=xxx\ --output=xxx_SER\ --outputtype=bin\ --stress=ndxyz)

  add_test(NAME ${arg}-p${nproc}-post_ensight-ser
    COMMAND sh -c " ${RUNPOSTFILTER_SER} && ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3 ${PROJECT_SOURCE_DIR}/tests/post_processing_test/ensight_comparison.py ${PROJECT_SOURCE_DIR}/Input/${arg}.dat xxx_SER_structure.case")
  set_tests_properties(${arg}-p${nproc}-post_ensight-ser PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})

  # additionally run postprocessing in parallel mode
  set(RUNPOSTFILTER_PAR ${MPI_RUN}\ -np\ ${nproc}\ ./post_drt_ensight\ --file=xxx\ --output=xxx_PAR\ --outputtype=bin\ --stress=ndxyz)

  add_test(NAME ${arg}-p${nproc}-post_ensight-par
    COMMAND sh -c " ${RUNPOSTFILTER_PAR} && ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3 ${PROJECT_SOURCE_DIR}/tests/post_processing_test/ensight_comparison.py ${PROJECT_SOURCE_DIR}/Input/${arg}.dat xxx_PAR_structure.case")
  set_tests_properties(${arg}-p${nproc}-post_ensight-par PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
endmacro (baci_test_and_post_ensight_test)

# Restart test case from test case previously run
macro (baci_test_restartonly arg nproc restart)
  # set additional output prefix identifier to empty string "" in default claase or to specific string if specified as optional input argument
  if(${ARGC} GREATER 3)
    set(IDENTIFIER ${ARGV3})
  else()
    set(IDENTIFIER "")
  endif()

  add_test(NAME ${arg}-p${nproc}-restart
    COMMAND ${MPI_RUN} -np ${nproc} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${arg}.dat xxx${IDENTIFIER} restart=${restart})

  set_tests_properties(${arg}-p${nproc}-restart PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
endmacro (baci_test_restartonly)

# Run test case for nested parallelism
macro (baci_test_Nested_Par arg1 arg2 restart)
  add_test(NAME ${arg1}-nestedPar
    COMMAND ${MPI_RUN} -np 3 $<TARGET_FILE:${baciname}> -ngroup=2 -glayout=1,2 -nptype=separateDatFiles ${PROJECT_SOURCE_DIR}/Input/${arg1}.dat xxx ${PROJECT_SOURCE_DIR}/Input/${arg2}.dat xxxAdditional)

  if (${restart})
    add_test(NAME ${arg1}-nestedPar-restart
    COMMAND ${MPI_RUN} -np 3 $<TARGET_FILE:${baciname}> -ngroup=2 -glayout=1,2 -nptype=separateDatFiles ${PROJECT_SOURCE_DIR}/Input/${arg1}.dat xxx restart=${restart} ${PROJECT_SOURCE_DIR}/Input/${arg2}.dat xxxAdditional restart=${restart})
  endif (${restart})

  set_tests_properties(${arg1}-nestedPar PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
endmacro (baci_test_Nested_Par)

# Run test case for nested parallelism with copydatfile
# arg1 is the inputfile
# arg2 is the number of procs
# arg3 is the number of groups
macro (baci_test_Nested_Par_CopyDat arg1 arg2 arg3)
  add_test(NAME ${arg1}-nestedPar_CopyDat-p${arg2}
    COMMAND ${MPI_RUN} -np ${arg2} $<TARGET_FILE:${baciname}> -ngroup=${arg3} -nptype=copyDatFile ${PROJECT_SOURCE_DIR}/Input/${arg1}.dat xxx )
    if( "${ARGN}" STREQUAL "minimal")
      set_tests_properties(${arg1}-nestedPar_CopyDat-p${arg2} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED} LABELS minimal)
    else ()
      set_tests_properties(${arg1}-nestedPar_CopyDat-p${arg2} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
    endif ()
endmacro (baci_test_Nested_Par_CopyDat)

# Run test case for nested parallelism with copydatfile
# + precursor simulation + (optional) "postprocesssing" simulation)
# arg1 is the inputfile of a precursor simulation
# arg2 is an inputfile relying on the output of arg1
# arg3 is the inputfile of a "postprocessing simulation" restarted from arg2 (optional)
# arg4 is the number of procs
# arg5 is the number of groups
# restart is a restart for arg2 to read from arg1-output
macro (baci_test_Nested_Par_CopyDat_prepost arg1 arg2 arg3 arg4 arg5 restart)
  # precursor simulation
  add_test(NAME ${arg2}_precursor-p${arg4}
    COMMAND ${MPI_RUN} -np ${arg4} $<TARGET_FILE:${baciname}> ${PROJECT_SOURCE_DIR}/Input/${arg1}.dat xxx)
  # restart from precursor output
  add_test(NAME ${arg2}-p${arg4}
    COMMAND ${MPI_RUN} -np ${arg4} $<TARGET_FILE:${baciname}> -ngroup=${arg5} -nptype=copyDatFile ${PROJECT_SOURCE_DIR}/Input/${arg2}.dat xxx)

  # add postprocessing simulation in case
  if (NOT ${arg3} STREQUAL "")
    add_test(NAME ${arg2}_postprocess-p${arg4}
      COMMAND ${MPI_RUN} -np ${arg4} $<TARGET_FILE:${baciname}> -ngroup=${arg5} -nptype=copyDatFile ${PROJECT_SOURCE_DIR}/Input/${arg3}.dat xxx restart=${restart})
  endif(NOT ${arg3} STREQUAL "")
  if( "${ARGN}" STREQUAL "minimal")
    set_tests_properties(${arg2}_precursor-p${arg4} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED} LABELS minimal)
    set_tests_properties(${arg2}-p${arg4} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED} LABELS minimal)
    if (NOT ${arg3} STREQUAL "")
      set_tests_properties(${arg2}_postprocess-p${arg4} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED} LABELS minimal)
    endif(NOT ${arg3} STREQUAL "")
  else ()
    set_tests_properties(${arg2}_precursor-p${arg4} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
    set_tests_properties(${arg2}-p${arg4} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
    if (NOT ${arg3} STREQUAL "")
      set_tests_properties(${arg2}_postprocess-p${arg4} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
    endif(NOT ${arg3} STREQUAL "")
  endif ()
endmacro (baci_test_Nested_Par_CopyDat_prepost)


# FRAMEWORK TESTS - testing the whole framework: from cubit via pre_exodus and baci to the post-filter
macro (baci_framework_test testname nproc xmlfilename)
  set (RUNCUBIT ${CUBIT_DIR}/cubit\ -batch\ -nographics\ -nojournal\ ${PROJECT_SOURCE_DIR}/tests/framework-test/${testname}.jou) # cubit is run to generate an exo file
  set (RUNPREEXODUS ./pre_exodus\ --exo=xxx_${testname}.e\ --bc=${PROJECT_SOURCE_DIR}/tests/framework-test/${testname}.bc\ --head=${PROJECT_SOURCE_DIR}/tests/framework-test/${testname}.head\ --dat=xxx.dat) # pre_exodus is run to generate a Dat file
  if (NOT ${xmlfilename} STREQUAL "")
    file(COPY ${PROJECT_SOURCE_DIR}/Input/${xmlfilename} DESTINATION ./) # if a XML file name is given, it is copied from the baci input directory to the build directory
  endif(NOT ${xmlfilename} STREQUAL "")
  set (RUNBACI ${MPI_RUN}\ -np\ ${nproc}\ $<TARGET_FILE:${baciname}>\ xxx.dat\ xxx) # baci is run using the generated dat file
  set (RUNPOSTFILTER ${MPI_RUN}\ -np\ ${nproc}\ ./post_drt_ensight\ --file=xxx) # post_drt_ensight is run for the resulting output
  add_test(NAME ${testname}-p${nproc}-fw
  COMMAND sh -c "${RUNCUBIT} && ${RUNPREEXODUS} && ${RUNBACI} && ${RUNPOSTFILTER}")

# note: for the clean-up job in the end, every generated intermediate file has to start with "xxx"

  set_tests_properties(${testname}-p${nproc}-fw PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
  set_tests_properties ( ${testname}-p${nproc}-fw PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR:; ERROR ;Error " )
  set_tests_properties ( ${testname}-p${nproc}-fw PROPERTIES ENVIRONMENT "PATH=$ENV{PATH}" )

endmacro (baci_framework_test)

# CUT TESTS
macro (cut_test nproc)
  set (RUNTESTS ${MPI_RUN}\ -np\ ${nproc}\ cut_test)
  add_test(NAME test-p${nproc}-cut
  COMMAND sh -c "${RUNTESTS}")
  set_tests_properties(test-p${nproc}-cut PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
  set_tests_properties ( test-p${nproc}-cut PROPERTIES FAIL_REGULAR_EXPRESSION "ERROR:; ERROR " )
endmacro (cut_test)

# POSTPROCESSING TEST
macro(post_processing arg nproc stresstype straintype startstep)
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
  set(RUNPOSTFILTER_SER ./post_drt_ensight\ --file=xxx${IDENTIFIER}\ --output=xxx${IDENTIFIER}_SER_${arg}\ --stress=${stresstype}\ --strain=${straintype}\ --start=${startstep})
  set(RUNPOSTFILTER_PAR ${MPI_RUN}\ -np\ ${nproc}\ ./post_drt_ensight\ --file=xxx${IDENTIFIER}\ --output=xxx${IDENTIFIER}_PAR_${arg}\ --stress=${stresstype}\ --strain=${straintype}\ --start=${startstep})

  # specify test case
  add_test(NAME ${arg}${IDENTIFIER}${FIELD}-p${nproc}-pp
    COMMAND sh -c " ${RUNPOSTFILTER_PAR} && ${RUNPOSTFILTER_SER} && ${PVPYTHON} ${PROJECT_SOURCE_DIR}/tests/post_processing_test/comparison.py xxx${IDENTIFIER}_PAR_${arg}${FIELD}*.case xxx${IDENTIFIER}_SER_${arg}${FIELD}*.case ${PROJECT_SOURCE_DIR}/Input/${arg}${IDENTIFIER}${FIELD}.csv")
  set_tests_properties(${arg}${IDENTIFIER}${FIELD}-p${nproc}-pp PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})
  set_tests_properties(${arg}${IDENTIFIER}${FIELD}-p${nproc}-pp PROPERTIES ENVIRONMENT "PATH=$ENV{PATH}")
endmacro(post_processing)

# compare arbitrary result file to corresponding reference file
macro(result_file arg nproc filetag resultfilename referencefilename tolerance)

  # add test to testing framework
  add_test(NAME ${arg}-p${nproc}-${filetag} COMMAND ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3 ${PROJECT_SOURCE_DIR}/utilities/diff_with_tolerance.py ${tolerance} ${PROJECT_BINARY_DIR}/${resultfilename} ${PROJECT_SOURCE_DIR}/Input/${referencefilename})

  # set maximum test runtime
  set_tests_properties(${arg}-p${nproc}-${filetag} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})

endmacro(result_file)

# compare XML formatted .vtk result data set referenced by .pvd files to corresponding reference files
macro(vtk_test name nproc pvd_resultfilename pvd_referencefilename tolerance)
# this test takes a list of times as extra arguments to check results at those timesteps
# if no extra arguments are given test checks every timestep

  set (extra_macro_args ${ARGN})
    # Did we get any optional args?
    list(LENGTH extra_macro_args num_extra_args)
    if (${num_extra_args} GREATER 0)
        list(GET extra_macro_args 0 optional_arg)
    endif ()
  # add test to testing framework
  add_test(NAME ${name}-p${nproc} COMMAND ${PROJECT_SOURCE_DIR}/utilities/baci-python-venv/bin/python3 ${PROJECT_SOURCE_DIR}/tests/output_test/vtk_compare.py ${PROJECT_BINARY_DIR}/${pvd_resultfilename} ${PROJECT_SOURCE_DIR}/Input/${pvd_referencefilename} ${tolerance} ${num_extra_args} ${extra_macro_args})

  # set maximum test runtime
  set_tests_properties(${name}-p${nproc} PROPERTIES TIMEOUT ${GLOBAL_TEST_TIMEOUT_SCALED})

endmacro(vtk_test)

###------------------------------------------------------------------ List of tests
include (TestingFrameworkListOfTests.cmake)


###------------------------------------------------------------------ Final cleanup
# remove any output files from our tests
# autogenerated core files (generated by kernel)
# special files generated by trilinos (such as amesos-failure.dat from ML)
add_test(NAME test_cleanup COMMAND sh -c "if [ -f *_CUTFAIL.pos ]; then mkdir -p ../cut-debug ; cp *_CUTFAIL.pos ../cut-debug/ ; fi ; rm -vfr xxx* core.* amesos-failure.dat")
