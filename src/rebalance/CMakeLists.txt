set(OBJS_REBALANCE
    # cmake-format: sortable
    drt_utils_rebalancing.cpp
    )

baci_add_library(rebalance ${OBJS_REBALANCE})

set(LIBS_REBALANCE
    # cmake-format: sortable
    linalg
    )

target_link_libraries(rebalance PUBLIC ${LIBS_REBALANCE})
