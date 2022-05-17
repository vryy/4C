# Add a unit test executable
#
# Usage: BACI_ADD_GOOGLE_TEST_EXECUTABLE(<name> [source1] [source2 ...])
#
function(BACI_ADD_GOOGLE_TEST_EXECUTABLE TESTNAME)
    add_executable(${TESTNAME} ${PROJECT_SOURCE_DIR}/unittests/gtest_main_mpi.cpp ${ARGN})
    target_include_directories(${TESTNAME} PRIVATE ${PROJECT_SOURCE_DIR})

    add_test(${TESTNAME} ${TESTNAME} --gtest_output=xml:unittest_reports/${TESTNAME}_report.xml)
    set_tests_properties(${TESTNAME} PROPERTIES TIMEOUT ${UNITTEST_TIMEOUT} LABELS minimal)

    add_dependencies(unittests ${TESTNAME})
endfunction()
