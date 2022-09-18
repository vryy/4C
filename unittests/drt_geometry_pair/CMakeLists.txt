set(TESTNAME unittests_drt_geometry_pair)

set(SOURCE_LIST
    unit_geometry_pair_line_to_surface.cpp
    unit_geometry_pair_line_to_surface_patch_averaged_normals.cpp
    unit_geometry_pair_line_to_volume_segmentation.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})

baci_link_google_test_necessary_libraries(
  ${TESTNAME}
  -Wl,--start-group
  ${BACI_LIBRARIES}
  -Wl,--end-group
  ${LIBRARIES}
  )
