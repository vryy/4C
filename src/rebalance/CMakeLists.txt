set(OBJS_REBALANCE
    # cmake-format: sortable
    rebalance.cpp
    rebalance_utils.cpp
    )

baci_add_library(rebalance ${OBJS_REBALANCE})

set(LIBS_REBALANCE
    # cmake-format: sortable
    linalg
    )

baci_add_dependency(rebalance ${LIBS_REBALANCE})
