set(OBJS_PSS_FULL
    # cmake-format: sortable
    pss_am.c
    pss_table.c
    pss_table_iter.c
    )

baci_add_library(pss_full ${OBJS_PSS_FULL})

set(LIBS_PSS_FULL
    # cmake-format: sortable
    headers
    )

target_link_libraries(pss_full PUBLIC ${LIBS_PSS_FULL})
