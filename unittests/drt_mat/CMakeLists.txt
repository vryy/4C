set(TESTNAME unittests_drt_mat)

set(SOURCE_LIST
    unit_elasthyper_service.cpp
    unit_electrode.cpp
    unit_inelastic_defgrad_factors.cpp
    unit_material_service.cpp
    unit_multiplicative_split_defgrad_elasthyper.cpp
    unit_multiplicative_split_defgrad_elasthyper_service.cpp
    unit_stvenantkirchhoff.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(
  ${TESTNAME}
  -Wl,--start-group
  ${BACI_LIBRARIES}
  -Wl,--end-group
  ${LIBRARIES}
  )
