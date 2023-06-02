## Link a unit test executable with necessary BACI libraries
#
# Usage: baci_link_google_test_necessary_libraries(<name> [lib1] [lib2 ...])
#
# where lib1, lib2, ... are libraries of BACI that are necessary to compile the test code.
# GoogleTest and GoogleMock are automatically linked and do not need to be specified here.
function(baci_link_google_test_necessary_libraries TESTNAME)
  # All libraries are linked as PRIVATE since a unit test executable cannot be used as a dependency itself.
  target_link_libraries(${TESTNAME} PRIVATE gtest gmock ${ARGN})

  # Link to common helpers for unit tests
  target_link_libraries(${TESTNAME} PRIVATE unittests_common)
endfunction()
