set(OBJS_DRT_FLUID_TURBULENCE
    dyn_smag.cpp
    dyn_vreman.cpp
    boxfilter.cpp
    drt_transfer_turb_inflow.cpp
    turbulence_hit_initial_field.cpp
    turbulence_hit_forcing.cpp
    turbulent_flow_algorithm.cpp
    turbulence_statistic_manager.cpp
    turbulence_statistics_bfs.cpp
    turbulence_statistics_ph.cpp
    turbulence_statistics_bfda.cpp
    turbulence_statistics_ccy.cpp
    turbulence_statistics_cha.cpp
    turbulence_statistics_hit.cpp
    turbulence_statistics_tgv.cpp
    turbulence_statistics_ldc.cpp
    turbulence_statistics_mean_general.cpp
    turbulence_statistics_sqc.cpp
    )
baci_add_library(alg_fluid_turbulence ${OBJS_DRT_FLUID_TURBULENCE})
target_link_libraries(alg_fluid_turbulence PUBLIC drt_lib)
