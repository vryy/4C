set(OBJS_FLUID_TURBULENCE
    # cmake-format: sortable
    fluid_turbulence_boxfilter.cpp
    fluid_turbulence_dyn_smag.cpp
    fluid_turbulence_dyn_vreman.cpp
    fluid_turbulence_hit_forcing.cpp
    fluid_turbulence_hit_initial_field.cpp
    fluid_turbulence_statistic_manager.cpp
    fluid_turbulence_statistics_bfda.cpp
    fluid_turbulence_statistics_ccy.cpp
    fluid_turbulence_statistics_cha.cpp
    fluid_turbulence_statistics_hit.cpp
    fluid_turbulence_statistics_ldc.cpp
    fluid_turbulence_statistics_mean_general.cpp
    fluid_turbulence_statistics_ph.cpp
    fluid_turbulence_statistics_sqc.cpp
    fluid_turbulence_statistics_tgv.cpp
    fluid_turbulence_transfer_turb_inflow.cpp
    fluid_turbulence_tstatistics_bfs.cpp
    fluid_turbulence_turbulent_flow_algorithm.cpp
    )
baci_add_library(alg_fluid_turbulence ${OBJS_FLUID_TURBULENCE})

set(LIBS_ALG_FLUID_TURBULENCE
    # cmake-format: sortable
    lib
    )

baci_add_dependency(alg_fluid_turbulence ${LIBS_ALG_FLUID_TURBULENCE})
