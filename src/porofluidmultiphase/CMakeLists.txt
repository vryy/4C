set(OBJS_POROFLUIDMULTIPHASE
    # cmake-format: sortable
    porofluidmultiphase_dyn.cpp
    porofluidmultiphase_meshtying_strategy_artery.cpp
    porofluidmultiphase_meshtying_strategy_std.cpp
    porofluidmultiphase_resulttest.cpp
    porofluidmultiphase_timint_implicit.cpp
    porofluidmultiphase_timint_ost.cpp
    porofluidmultiphase_utils.cpp
    )
baci_add_library(alg_porofluidmultiphase ${OBJS_POROFLUIDMULTIPHASE})

set(LIBS_ALG_POROFLUIDMULTIPHASE
    # cmake-format: sortable
    inpar
    )

target_link_libraries(alg_porofluidmultiphase PUBLIC ${LIBS_ALG_POROFLUIDMULTIPHASE})
