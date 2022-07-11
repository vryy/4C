# Link a unit test executable with necessary BACI libraries
#
# Usage: BACI_LINK_GOOGLE_TEST_NECESSARY_LIBRARIES(<name> [lib1] [lib2 ...])
#
# where lib1, lib2, ... are libraries of BACI that are necessary to compile the test code.
# All TPLs as well as GoogleTest are automatically linked and do not need to be specified here.
function(BACI_LINK_GOOGLE_TEST_NECESSARY_LIBRARIES TESTNAME)
  target_link_libraries(
    ${TESTNAME}
    gtest
    gmock
    ${ARGN}
    ${LIBRARIES}
    )

  # Link to a special version of dserror which throws a std::runtime_exception
  target_link_libraries(${TESTNAME} drt_error_for_testing)
endfunction()
