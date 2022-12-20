set(OBJS_POST_VTK
    # cmake-format: sortable
    post_drt_vti_writer.cpp
    post_drt_vtk_writer.cpp
    post_drt_vtu_writer.cpp
    post_drt_vtu_writer_node_based.cpp
    )

baci_add_library(post_drt_vtk_lib ${OBJS_POST_VTK})
target_link_libraries(post_drt_vtk_lib PUBLIC post_drt_common)
