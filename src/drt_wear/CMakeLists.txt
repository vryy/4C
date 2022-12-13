set(OBJS_DRT_WEAR
    # cmake-format: sortable
    wear_algorithm.cpp
    wear_dyn.cpp
    wear_partitioned.cpp
    )

baci_add_library(alg_wear ${OBJS_DRT_WEAR})

set(LIBS_ALG_WEAR
    # cmake-format: sortable
    alg_fs3i
    drt_adapter
    )

target_link_libraries(alg_wear PUBLIC ${LIBS_ALG_WEAR})
