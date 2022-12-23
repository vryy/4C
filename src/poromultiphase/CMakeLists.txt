set(OBJS_POROMULTIPHASE
    # cmake-format: sortable
    poromultiphase_base.cpp
    poromultiphase_dyn.cpp
    poromultiphase_monolithic.cpp
    poromultiphase_monolithic_twoway.cpp
    poromultiphase_partitioned.cpp
    poromultiphase_partitioned_twoway.cpp
    poromultiphase_utils.cpp
    poromultiphase_utils_clonestrategy.cpp
    )

baci_add_library(alg_poromultiphase ${OBJS_POROMULTIPHASE})

set(LIBS_ALG_PROMULTIPHASE
    # cmake-format: sortable
    inpar
    )

target_link_libraries(alg_poromultiphase PUBLIC ${LIBS_ALG_PROMULTIPHASE})
