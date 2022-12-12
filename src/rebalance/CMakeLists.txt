set(OBJS_REBALANCE drt_utils_rebalancing.cpp)

baci_add_library(rebalance ${OBJS_REBALANCE})
target_link_libraries(rebalance PUBLIC linalg)
