set(OBJS_POST_VTK
    # cmake-format: sortable
    post_vtk_vti_writer.cpp
    post_vtk_vtu_writer.cpp
    post_vtk_vtu_writer_node_based.cpp
    post_vtk_writer.cpp
    )

baci_add_library(post_vtk_lib ${OBJS_POST_VTK})
baci_add_dependency(post_vtk_lib post_common)
