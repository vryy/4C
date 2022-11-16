set(OBJS_DRT_TSI
    tsi_algorithm.cpp
    tsi_monolithic.cpp
    tsi_partitioned.cpp
    tsi_dyn.cpp
    tsi_utils.cpp
    )

baci_add_library(alg_tsi ${OBJS_DRT_TSI})
target_link_libraries(alg_tsi PUBLIC drt_adapter drt_volmortar)
