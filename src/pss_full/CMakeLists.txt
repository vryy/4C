set(OBJS_PSS_FULL
    # cmake-format: sortable
    pss_full_am.c
    pss_full_table.c
    pss_full_table_iter.c
    )

baci_add_library(pss_full ${OBJS_PSS_FULL})

set(LIBS_PSS_FULL
    # cmake-format: sortable
    headers
    )

baci_add_dependency(pss_full ${LIBS_PSS_FULL})
