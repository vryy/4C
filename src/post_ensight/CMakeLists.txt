set(OBJS_POST_ENSIGHT
    # cmake-format: sortable
    baci_post_ensight_writer.cpp
    baci_post_ensight_writer_nurbs.cpp
    )

baci_add_library(
  post_ensight_lib
  SOURCES
  ${OBJS_POST_ENSIGHT}
  HEADERS
  baci_post_ensight_writer.H
  )
baci_add_dependency(post_ensight_lib post_common)
