set(OBJS_POST_DRT_ENSIGHT
    # cmake-format: sortable
    post_drt_ensight_writer.cpp
    post_drt_ensight_writer_nurbs.cpp
    )

baci_add_library(post_drt_ensight_lib ${OBJS_POST_DRT_ENSIGHT})
target_link_libraries(post_drt_ensight_lib PUBLIC post_drt_common)
