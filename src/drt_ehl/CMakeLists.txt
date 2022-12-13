set(OBJS_DRT_EHL
    # cmake-format: sortable
    ehl_base.cpp
    ehl_dyn.cpp
    ehl_monolithic.cpp
    ehl_partitioned.cpp
    ehl_utils.cpp
    )

baci_add_library(alg_ehl ${OBJS_DRT_EHL})

set(LIBS_ALG_EHL
    # cmake-format: sortable
    linalg
    )

target_link_libraries(alg_ehl PUBLIC ${LIBS_ALG_EHL})
