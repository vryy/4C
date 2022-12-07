set(OBJS_DRT_PASI
    pasi_dyn.cpp
    pasi_partitioned_onewaycoup.cpp
    pasi_partitioned_twowaycoup.cpp
    pasi_partitioned.cpp
    pasi_str_model_evaluator_partitioned.cpp
    pasi_utils.cpp
    )

baci_add_library(alg_pasi ${OBJS_DRT_PASI})
target_link_libraries(alg_pasi PUBLIC drt_adapter)
