set(OBJS_POST_ENSIGHT
    # cmake-format: sortable
    post_ensight_writer.cpp
    post_ensight_writer_nurbs.cpp
    )

baci_add_library(post_ensight_lib ${OBJS_POST_ENSIGHT})
baci_add_dependency(post_ensight_lib post_common)
