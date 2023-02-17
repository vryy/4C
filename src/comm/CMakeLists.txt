set(OBJS_COMM
    # cmake-format: sortable
    comm_utils.cpp
    )

baci_add_library(comm ${OBJS_COMM})

set(LIBS_COMM
    # cmake-format: sortable
    lib
    nurbs_discret
    rebalance
    )

baci_add_dependency(comm ${LIBS_COMM})
