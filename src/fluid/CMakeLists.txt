set(OBJS_FLUID
    # cmake-format: sortable
    fluid_coupling_red_models.cpp
    fluid_discret_extractor.cpp
    fluid_discretization_runtime_vtu_output_params.cpp
    fluid_dyn_nln_drt.cpp
    fluid_functions.cpp
    fluid_impedancecondition.cpp
    fluid_implicit_integration.cpp
    fluid_meshtying.cpp
    fluid_MHD_evaluate.cpp
    fluid_result_test.cpp
    fluid_rotsym_periodicbc_utils.cpp
    fluid_timint.cpp
    fluid_timint_ac.cpp
    fluid_timint_ac_ost.cpp
    fluid_timint_bdf2.cpp
    fluid_timint_genalpha.cpp
    fluid_timint_hdg.cpp
    fluid_timint_hdg_weak_comp.cpp
    fluid_timint_loma.cpp
    fluid_timint_loma_bdf2.cpp
    fluid_timint_loma_genalpha.cpp
    fluid_timint_loma_ost.cpp
    fluid_timint_ost.cpp
    fluid_timint_poro.cpp
    fluid_timint_poro_genalpha.cpp
    fluid_timint_poro_ost.cpp
    fluid_timint_poro_stat.cpp
    fluid_timint_red.cpp
    fluid_timint_red_bdf2.cpp
    fluid_timint_red_genalpha.cpp
    fluid_timint_red_ost.cpp
    fluid_timint_red_stat.cpp
    fluid_timint_stat.cpp
    fluid_timint_stat_hdg.cpp
    fluid_timint_two_phase.cpp
    fluid_timint_two_phase_genalpha.cpp
    fluid_timint_two_phase_ost.cpp
    fluid_timint_two_phase_stat.cpp
    fluid_utils.cpp
    fluid_utils_infnormscaling.cpp
    fluid_utils_mapextractor.cpp
    fluid_volumetric_surfaceFlow_condition.cpp
    fluid_xwall.cpp
    )

baci_add_library(alg_fluid ${OBJS_FLUID})

set(LIBS_ALG_FLUID
    # cmake-format: sortable
    lib
    lib_artnet
    rebalance
    red_airways
    )

baci_add_dependency(alg_fluid ${LIBS_ALG_FLUID})
