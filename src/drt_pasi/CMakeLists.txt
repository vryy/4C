set(OBJS_DRT_PASI
    # cmake-format: sortable
    pasi_dyn.cpp
    pasi_partitioned.cpp
    pasi_partitioned_onewaycoup.cpp
    pasi_partitioned_twowaycoup.cpp
    pasi_str_model_evaluator_partitioned.cpp
    pasi_utils.cpp
    )

baci_add_library(alg_pasi ${OBJS_DRT_PASI})

set(LIBS_ALG_PASI
    # cmake-format: sortable
    drt_adapter
    )

target_link_libraries(alg_pasi PUBLIC ${LIBS_ALG_PASI})
