set(OBJS_POST_VTK
    # cmake-format: sortable
    baci_post_vtk_vti_writer.cpp
    baci_post_vtk_vtu_writer.cpp
    baci_post_vtk_vtu_writer_node_based.cpp
    baci_post_vtk_writer.cpp
    )

set(_headers
    baci_post_vtk_vti_writer.H
    baci_post_vtk_vtu_writer.H
    baci_post_vtk_vtu_writer_node_based.H
    baci_post_vtk_writer.H
    )

baci_add_library(post_vtk_lib SOURCES ${OBJS_POST_VTK} HEADERS)
baci_add_dependency(post_vtk_lib post_common)
