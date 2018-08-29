# parameters
set(CTEST_SOURCE_DIRECTORY "$ENV{CI_PROJECT_DIR}")
set(CTEST_BINARY_DIRECTORY "$ENV{CI_PROJECT_DIR}/../baci-build")
set(CTEST_TIMEOUT          "1000")       # ctest timeout 500s
set($ENV{LC_MESSAGES}      "en_EN" )    # set output to english such that ctest can analyze it

set(CTEST_SITE "$ENV{HOSTNAME}")
set(CTEST_BUILD_NAME "$ENV{CTEST_BUILD_NAME_GITLAB}")

# prepare environment
set(CTEST_CMAKE_GENERATOR "Unix Makefiles")
set(WITH_MEMCHECK FALSE)
set(WITH_COVERAGE FALSE)

# ctest_empty_binary_directory(${CTEST_BINARY_DIRECTORY})

# prepare commands
set(CTEST_CONFIGURE_COMMAND "$ENV{CTEST_CONFIGURE_PREFIX}${CTEST_SOURCE_DIRECTORY}/do-configure -f --config=${CTEST_SOURCE_DIRECTORY}/buildconfig/$ENV{CTEST_BUILD_CONFIG_GITLAB} $ENV{CTEST_CONFIGURE_POSTFIX}")
set(CTEST_BUILD_COMMAND     "$ENV{CTEST_MAKE}")
set(CTEST_CUSTOM_MAXIMUM_NUMBER_OF_WARNINGS "5000")

# do the testing
ctest_start("${HOSTNAME}")
#ctest_update()
ctest_configure()
ctest_build(RETURN_VALUE testBuild)
#ctest_test()

if ($ENV{CTEST_DROP_SITE_CDASH_GITLAB} EQUAL 1)
    ctest_submit()
endif ($ENV{CTEST_DROP_SITE_CDASH_GITLAB} EQUAL 1)

if (NOT ${testBuild} EQUAL 0)
  message( SEND_ERROR "Baci build failed!" )
endif (NOT ${testBuild} EQUAL 0)

