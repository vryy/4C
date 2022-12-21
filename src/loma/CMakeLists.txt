set(OBJS_LOMA
    # cmake-format: sortable
    loma_algorithm.cpp
    loma_dyn.cpp
    )

baci_add_library(alg_loma ${OBJS_LOMA})

set(LIBS_ALG_LOMA
    # cmake-format: sortable
    linalg
    )

target_link_libraries(alg_loma PUBLIC ${LIBS_ALG_LOMA})
