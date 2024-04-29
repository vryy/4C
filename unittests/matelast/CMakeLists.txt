set(TESTNAME unittests_matelast)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_elast_coupanisoexpo_elementfibers_test.cpp
    4C_elast_coupanisoexpo_gausspointfibers_test.cpp
    4C_elast_coupanisoexpobase_test.cpp
    4C_elast_coupanisoexposhear_elementfibers_test.cpp
    4C_elast_coupanisoexposhear_gausspointfibers_test.cpp
    4C_elast_isoogden_test.cpp
    )

four_c_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
