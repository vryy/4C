set(TESTNAME unittests_shell_kl_nurbs)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_shell_kl_nurbs_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
