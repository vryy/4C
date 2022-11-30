set(OBJS_DRT_WEAR wear_algorithm.cpp wear_dyn.cpp wear_partitioned.cpp)

baci_add_library(alg_wear ${OBJS_DRT_WEAR})
target_link_libraries(alg_wear PUBLIC drt_adapter)
