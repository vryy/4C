set(OBJS
    # cmake-format: sortable
    io_legacy_table.c
    io_legacy_table_iter.c
    )

baci_add_library(io_legacy ${OBJS})
