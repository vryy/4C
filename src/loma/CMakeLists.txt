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

baci_add_dependency(alg_loma ${LIBS_ALG_LOMA})
