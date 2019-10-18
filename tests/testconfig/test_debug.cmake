# parameters
set(CTEST_SOURCE_DIRECTORY "$ENV{CI_PROJECT_DIR}")
set(CTEST_BINARY_DIRECTORY "$ENV{CI_PROJECT_DIR}/../baci-build")
set($ENV{LC_MESSAGES}      "en_EN" )    # set output to english such that ctest can analyze it

set(CTEST_SITE "$ENV{HOSTNAME}")
set(CTEST_BUILD_NAME "$ENV{CTEST_BUILD_NAME_GITLAB}")

# prepare environment
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(WITH_MEMCHECK TRUE)
set(WITH_COVERAGE TRUE)

ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

# prepare commands
set(CTEST_CONFIGURE_COMMAND "$ENV{CTEST_CONFIGURE_PREFIX}${CTEST_SOURCE_DIRECTORY}/do-configure -f --debug-optimized --trilinos-debug --config=${CTEST_SOURCE_DIRECTORY}/buildconfig/$ENV{CTEST_BUILD_CONFIG_GITLAB} $ENV{CTEST_CONFIGURE_POSTFIX}")
set(CTEST_BUILD_COMMAND     "$ENV{CTEST_MAKE}")
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS "5000")
set(CTEST_COMMAND "ctest -D Nightly")

# do the testing
ctest_start("${HOSTNAME}")
#ctest_update()
ctest_configure()
ctest_build(RETURN_VALUE testBuild NUMBER_WARNINGS numWarnings)
if (NOT $ENV{TEST_TAG} STREQUAL "") # Only execute if tests are selected as otherwise scripte returns err
  ctest_test(RETURN_VALUE testRes)
endif (NOT $ENV{TEST_TAG} STREQUAL "")

if ($ENV{CTEST_DROP_SITE_CDASH_GITLAB} EQUAL 1) # Send results to the Dashboard (cdash) if enabled
    ctest_submit()
endif ($ENV{CTEST_DROP_SITE_CDASH_GITLAB} EQUAL 1)

if (NOT ${testBuild} EQUAL 0) # Send error for a failed build
  message( SEND_ERROR "Baci build failed!" )
endif (NOT ${testBuild} EQUAL 0)

if (NOT ${numWarnings} EQUAL 0 AND $ENV{CTEST_FAIL_ON_WARNING} EQUAL 1) # Send error for warnings if enabled
  message( SEND_ERROR "Baci build issued build warnings!" )
endif (NOT ${numWarnings} EQUAL 0 AND $ENV{CTEST_FAIL_ON_WARNING} EQUAL 1)

if (NOT ${testRes} EQUAL 0 AND NOT $ENV{TEST_TAG} STREQUAL "") # Send error if test fail (if test were performed)
  message( SEND_ERROR "Baci tests failed!" )
endif (NOT ${testRes} EQUAL 0 AND NOT $ENV{TEST_TAG} STREQUAL "")

