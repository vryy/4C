set(TESTNAME unittests_functions)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_function_library_test.cpp
    baci_functionvariables_test.cpp
    baci_symbolic_expression_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_target_requires_cyclic_legacy_linking(${TESTNAME})
