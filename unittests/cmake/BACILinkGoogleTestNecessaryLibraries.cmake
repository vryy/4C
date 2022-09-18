# Link a unit test executable with necessary BACI libraries
#
# Usage: baci_link_google_test_necessary_libraries(<name> [lib1] [lib2 ...])
#
# where lib1, lib2, ... are libraries of BACI that are necessary to compile the test code.
# All TPLs as well as GoogleTest are automatically linked and do not need to be specified here.
function(baci_link_google_test_necessary_libraries TESTNAME)
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
