set(OBJS_DRT_FLUID_TURBULENCE
    # cmake-format: sortable
    boxfilter.cpp
    drt_transfer_turb_inflow.cpp
    dyn_smag.cpp
    dyn_vreman.cpp
    turbulence_hit_forcing.cpp
    turbulence_hit_initial_field.cpp
    turbulence_statistic_manager.cpp
    turbulence_statistics_bfda.cpp
    turbulence_statistics_bfs.cpp
    turbulence_statistics_ccy.cpp
    turbulence_statistics_cha.cpp
    turbulence_statistics_hit.cpp
    turbulence_statistics_ldc.cpp
    turbulence_statistics_mean_general.cpp
    turbulence_statistics_ph.cpp
    turbulence_statistics_sqc.cpp
    turbulence_statistics_tgv.cpp
    turbulent_flow_algorithm.cpp
    )
baci_add_library(alg_fluid_turbulence ${OBJS_DRT_FLUID_TURBULENCE})

set(LIBS_ALG_FLUID_TURBULENCE
    # cmake-format: sortable
    drt_lib
    )

target_link_libraries(alg_fluid_turbulence PUBLIC ${LIBS_ALG_FLUID_TURBULENCE})
