set(OBJS_DRT_POROELAST
    # cmake-format: sortable
    poro_base.cpp
    poro_dyn.cpp
    poro_monolithicfluidsplit.cpp
    poro_monolithicmeshtying.cpp
    poro_monolithicsplit.cpp
    poro_monolithicsplit_nopenetration.cpp
    poro_monolithicstructuresplit.cpp
    poro_partitioned.cpp
    poro_scatra_base.cpp
    poro_scatra_monolithic.cpp
    poro_scatra_part.cpp
    poro_scatra_part_1wc.cpp
    poro_scatra_part_2wc.cpp
    poro_utils_clonestrategy.cpp
    poroelast_monolithic.cpp
    poroelast_utils.cpp
    )
baci_add_library(alg_poroelast ${OBJS_DRT_POROELAST})

set(LIBS_ALG_POROELAST
    # cmake-format: sortable
    drt_inpar
    ele_so
    ele_w1
    ele_w1
    )

target_link_libraries(alg_poroelast PUBLIC ${LIBS_ALG_POROELAST})
