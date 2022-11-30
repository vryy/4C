set(OBJS_DRT_POROFLUIDMULTIPHASE
    porofluidmultiphase_dyn.cpp
    porofluidmultiphase_timint_implicit.cpp
    porofluidmultiphase_timint_ost.cpp
    porofluidmultiphase_meshtying_strategy_std.cpp
    porofluidmultiphase_meshtying_strategy_artery.cpp
    porofluidmultiphase_resulttest.cpp
    porofluidmultiphase_utils.cpp
    )
baci_add_library(alg_porofluidmultiphase ${OBJS_DRT_POROFLUIDMULTIPHASE})

target_link_libraries(alg_porofluidmultiphase PUBLIC drt_lib_inpar)
