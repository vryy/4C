if(${FOUR_C_BUILD_TYPE_UPPER} EQUAL "DEBUG")
  set(FOUR_C_TEST_TIMEOUT_SCALE
      4
      CACHE STRING "Scale timeout of tests by this factor."
      )
else()
  set(FOUR_C_TEST_TIMEOUT_SCALE
      1
      CACHE STRING "Scale timeout of tests by this factor."
      )
endif()
message(STATUS "Global test timeout scale is ${FOUR_C_TEST_TIMEOUT_SCALE}.")

math(EXPR FOUR_C_TEST_GLOBAL_TIMEOUT "120*${FOUR_C_TEST_TIMEOUT_SCALE}")
message(STATUS "The scaled global test timeout is ${FOUR_C_TEST_GLOBAL_TIMEOUT} s.")

# Fetch GoogleTest and setup the unit tests if option is enabled
four_c_process_global_option(FOUR_C_WITH_GOOGLETEST "Use GoogleTest for unit testing" ON)
if(FOUR_C_WITH_GOOGLETEST)
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

  math(EXPR UNITTEST_TIMEOUT "10*${FOUR_C_TEST_TIMEOUT_SCALE}")
  message(STATUS "The scaled unit test timeout is ${UNITTEST_TIMEOUT} s.")

else()
  message(STATUS "Unit tests with GoogleTest: disabled")
endif()
