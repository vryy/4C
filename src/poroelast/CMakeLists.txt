set(OBJS_POROELAST
    # cmake-format: sortable
    poroelast_base.cpp
    poroelast_dyn.cpp
    poroelast_monolithic.cpp
    poroelast_monolithicfluidsplit.cpp
    poroelast_monolithicmeshtying.cpp
    poroelast_monolithicsplit.cpp
    poroelast_monolithicsplit_nopenetration.cpp
    poroelast_monolithicstructuresplit.cpp
    poroelast_partitioned.cpp
    poroelast_utils.cpp
    poroelast_utils_clonestrategy.cpp
    )
baci_add_library(alg_poroelast ${OBJS_POROELAST})

set(LIBS_ALG_POROELAST
    # cmake-format: sortable
    ele_so
    ele_w1
    ele_w1
    inpar
    )

baci_add_dependency(alg_poroelast ${LIBS_ALG_POROELAST})
