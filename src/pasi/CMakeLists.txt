set(OBJS_PASI
    # cmake-format: sortable
    pasi_dyn.cpp
    pasi_partitioned.cpp
    pasi_partitioned_onewaycoup.cpp
    pasi_partitioned_twowaycoup.cpp
    pasi_str_model_evaluator_partitioned.cpp
    pasi_utils.cpp
    )

baci_add_library(alg_pasi ${OBJS_PASI})

set(LIBS_ALG_PASI
    # cmake-format: sortable
    adapter
    )

baci_add_dependency(alg_pasi ${LIBS_ALG_PASI})
