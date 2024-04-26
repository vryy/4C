set(TESTNAME unittests_functions)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_function_library_test.cpp
    4C_functionvariables_test.cpp
    4C_symbolic_expression_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
