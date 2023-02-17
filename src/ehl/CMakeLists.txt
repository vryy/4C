set(OBJS_EHL
    # cmake-format: sortable
    ehl_base.cpp
    ehl_dyn.cpp
    ehl_monolithic.cpp
    ehl_partitioned.cpp
    ehl_utils.cpp
    )

baci_add_library(alg_ehl ${OBJS_EHL})

set(LIBS_ALG_EHL
    # cmake-format: sortable
    linalg
    )

baci_add_dependency(alg_ehl ${LIBS_ALG_EHL})
