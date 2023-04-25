set(OBJS_WEAR
    # cmake-format: sortable
    wear_algorithm.cpp
    wear_dyn.cpp
    wear_partitioned.cpp
    )

baci_add_library(alg_wear ${OBJS_WEAR})

set(LIBS_ALG_WEAR
    # cmake-format: sortable
    adapter
    alg_fs3i
    )

baci_add_dependency(alg_wear ${LIBS_ALG_WEAR})
