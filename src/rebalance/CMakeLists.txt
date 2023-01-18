set(OBJS_REBALANCE
    # cmake-format: sortable
    rebalance.cpp
    )

baci_add_library(rebalance ${OBJS_REBALANCE})

set(LIBS_REBALANCE
    # cmake-format: sortable
    linalg
    )

target_link_libraries(rebalance PUBLIC ${LIBS_REBALANCE})
