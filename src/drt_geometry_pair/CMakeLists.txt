set(OBJS_DRT_GEOMETRY_PAIR
    # cmake-format: sortable
    geometry_pair.cpp
    geometry_pair_element_faces.cpp
    geometry_pair_evaluation_data_base.cpp
    geometry_pair_factory.cpp
    geometry_pair_line_projection.cpp
    geometry_pair_line_to_3D_evaluation_data.cpp
    geometry_pair_line_to_surface.cpp
    geometry_pair_line_to_surface_evaluation_data.cpp
    geometry_pair_line_to_surface_gauss_point_projection.cpp
    geometry_pair_line_to_surface_segmentation.cpp
    geometry_pair_line_to_volume.cpp
    geometry_pair_line_to_volume_gauss_point_projection.cpp
    geometry_pair_line_to_volume_gauss_point_projection_cross_section.cpp
    geometry_pair_line_to_volume_segmentation.cpp
    geometry_pair_utility_functions.cpp
    )

baci_add_library(drt_geometry_pair ${OBJS_DRT_GEOMETRY_PAIR})

set(LIBS_DRT_GEOMETRY_PAIR
    # cmake-format: sortable
    drt_nurbs_discret
    linalg
    )

target_link_libraries(drt_geometry_pair PUBLIC ${LIBS_DRT_GEOMETRY_PAIR})
