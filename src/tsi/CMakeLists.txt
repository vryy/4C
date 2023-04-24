set(OBJS_TSI
    # cmake-format: sortable
    tsi_algorithm.cpp
    tsi_dyn.cpp
    tsi_monolithic.cpp
    tsi_partitioned.cpp
    tsi_utils.cpp
    )

baci_add_library(alg_tsi ${OBJS_TSI})

set(LIBS_ALG_TSI
    # cmake-format: sortable
    adapter
    ele_thermo
    volmortar
    )

baci_add_dependency(alg_tsi ${LIBS_ALG_TSI})
