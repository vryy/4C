set(TESTNAME unittests_drt_so3)

set(SOURCE_LIST
    unit_so_hex8.cpp
    unit_so_tet4.cpp
    unit_so_element_service.cpp
    unit_so_hex8_determinant_analysis.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(
  ${TESTNAME} -Wl,--start-group ${BACI_LIBRARIES} -Wl,--end-group
  )
