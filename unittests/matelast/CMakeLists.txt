set(TESTNAME unittests_matelast)

set(SOURCE_LIST
    # cmake-format: sortable
    baci_elast_coupanisoexpo_elementfibers_test.cpp
    baci_elast_coupanisoexpo_gausspointfibers_test.cpp
    baci_elast_coupanisoexpobase_test.cpp
    baci_elast_coupanisoexposhear_elementfibers_test.cpp
    baci_elast_coupanisoexposhear_gausspointfibers_test.cpp
    baci_elast_isoogden_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(${TESTNAME})
baci_target_requires_cyclic_legacy_linking(${TESTNAME})
