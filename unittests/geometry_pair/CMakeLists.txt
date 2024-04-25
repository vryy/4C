set(TESTNAME unittests_geometry_pair)

set(SOURCE_LIST
    # cmake-format: sortable
    4C_geometry_pair_line_to_surface_patch_averaged_normals_test.cpp
    4C_geometry_pair_line_to_surface_test.cpp
    4C_geometry_pair_line_to_volume_segmentation_test.cpp
    )

baci_add_google_test_executable(${TESTNAME} SOURCE ${SOURCE_LIST})
