set(OBJS_DRT_POROELAST
    poro_dyn.cpp
    poro_base.cpp
    poroelast_monolithic.cpp
    poroelast_utils.cpp
    poro_utils_clonestrategy.cpp
    poro_partitioned.cpp
    poro_monolithicsplit.cpp
    poro_monolithicstructuresplit.cpp
    poro_monolithicfluidsplit.cpp
    poro_monolithicsplit_nopenetration.cpp
    poro_monolithicmeshtying.cpp
    poro_scatra_base.cpp
    poro_scatra_part.cpp
    poro_scatra_part_1wc.cpp
    poro_scatra_part_2wc.cpp
    poro_scatra_monolithic.cpp
    )
baci_add_library(alg_poroelast ${OBJS_DRT_POROELAST})
target_link_libraries(alg_poroelast PUBLIC drt_inpar)
