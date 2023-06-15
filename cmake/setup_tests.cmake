# Determine timeout for each test. Use default one, if it is not passed from the outside.
if(DEFINED ENV{GLOBAL_TEST_TIMEOUT})
  set(GLOBAL_TEST_TIMEOUT $ENV{GLOBAL_TEST_TIMEOUT})
  message(STATUS "Global test timeout is $ENV{GLOBAL_TEST_TIMEOUT} s (before scaling).")
else()
  # default test timeout, if not passed as an environment variable
  set(GLOBAL_TEST_TIMEOUT 260) # Default timeout

  message(
    STATUS
      "Global test timeout is not passed as an environment variable. It is set to the default ${GLOBAL_TEST_TIMEOUT} s (before scaling)."
    )
endif()

# Determine timeout scale factor. Use default one, if it is not passed from the outside
if(DEFINED ENV{GLOBAL_TEST_TIMEOUT_SCALE})
  set(GLOBAL_TEST_TIMEOUT_SCALE $ENV{GLOBAL_TEST_TIMEOUT_SCALE})
  message(STATUS "Global test timeout scale is $ENV{GLOBAL_TEST_TIMEOUT_SCALE}.")
else()
  # default test timeout scale, if not passed as an environment variable
  if("${CMAKE_BUILD_TYPE}" STREQUAL "DEBUG")
    set(GLOBAL_TEST_TIMEOUT_SCALE 4) # Default timeout scale for debug configuration
  else()
    set(GLOBAL_TEST_TIMEOUT_SCALE 1) # Default timeout scale
  endif()

  message(
    STATUS
      "Global test timeout scale is not passed as an environment variable. It is set to the default ${GLOBAL_TEST_TIMEOUT_SCALE} for this kind of build."
    )
endif()

# Determine scaled global test timeout
math(EXPR GLOBAL_TEST_TIMEOUT_SCALED "${GLOBAL_TEST_TIMEOUT}*${GLOBAL_TEST_TIMEOUT_SCALE}")
message(STATUS "The scaled global test timeout is ${GLOBAL_TEST_TIMEOUT_SCALED} s.")

# Fetch GoogleTest and setup the unit tests if option is enabled
option(BACI_WITH_GOOGLETEST "Use GoogleTest for unit testing" ON)
if(BACI_WITH_GOOGLETEST)
  message(STATUS "Unit tests with GoogleTest: enabled")

  # Define a convenience target for all unit tests and add it to 'full'
  # All unit test executables should add themselves as a dependency to 'unittests'
  add_custom_target(unittests)
  add_dependencies(full unittests)

  if(TARGET gtest)
    message(
      FATAL_ERROR "A target <gtest> has already been included by a TPL." "This is not supported."
      )
  endif()

  include(FetchContent)
  fetchcontent_declare(
    googletest
    GIT_REPOSITORY https://github.com/google/googletest.git
    GIT_TAG release-1.11.0
    )
  fetchcontent_makeavailable(googletest)

  if(DEFINED ENV{UNITTEST_TIMEOUT})
    set(UNITTEST_TIMEOUT $ENV{UNITTEST_TIMEOUT})
    message(STATUS "GoogleTest: The unit test timeout is $ENV{UNITTEST_TIMEOUT} s.")
  else()
    # default unit test timeout, if not passed as an environment variable
    set(UNITTEST_TIMEOUT 10) # Default unit test timeout

    # Determine scaled global test timeout
    math(EXPR UNITTEST_TIMEOUT "10*${GLOBAL_TEST_TIMEOUT_SCALE}")
    message(STATUS "The scaled global test timeout is ${UNITTEST_TIMEOUT} s.")

    message(
      STATUS
        "GoogleTest: The unit test timeout is not passed as an environment variable. "
        "It is set to the default value 10 s * ${GLOBAL_TEST_TIMEOUT_SCALE} = ${UNITTEST_TIMEOUT} s."
      )
  endif()

else()
  message(STATUS "Unit tests with GoogleTest: disabled")
endif()
