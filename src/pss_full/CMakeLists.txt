set(OBJS_PSS_FULL pss_am.c pss_table.c pss_table_iter.c)

baci_add_library(pss_full ${OBJS_PSS_FULL})
target_link_libraries(pss_full PUBLIC headers)
