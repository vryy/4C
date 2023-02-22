set(OBJS_POROELAST_SCATRA
    # cmake-format: sortable
    poroelast_scatra_base.cpp
    poroelast_scatra_dyn.cpp
    poroelast_scatra_monolithic.cpp
    poroelast_scatra_part.cpp
    poroelast_scatra_part_1wc.cpp
    poroelast_scatra_part_2wc.cpp
    poroelast_scatra_utils.cpp
    poroelast_scatra_utils_clonestrategy.cpp
    )
baci_add_library(alg_poroelast_scatra ${OBJS_POROELAST_SCATRA})

set(LIBS_ALG_POROELAST_SCATRA
    # cmake-format: sortable
    alg_poroelast
    lib
    )

target_link_libraries(alg_poroelast_scatra PUBLIC ${LIBS_ALG_POROELAST_SCATRA})
