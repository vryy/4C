set(OBJS_STATINV_ANALYSIS
    # cmake-format: sortable
    inv_analysis_auglagr.cpp
    inv_analysis_base.cpp
    inv_analysis_cal_drt.cpp
    inv_analysis_chol_factor.cpp
    inv_analysis_chol_factory.cpp
    inv_analysis_control.cpp
    inv_analysis_DcsMatrix.cpp
    inv_analysis_factory.cpp
    inv_analysis_initial_guess.cpp
    inv_analysis_likelihood_evaluation.cpp
    inv_analysis_matpar_manager.cpp
    inv_analysis_matpar_manager_elementwise.cpp
    inv_analysis_matpar_manager_patchwise.cpp
    inv_analysis_matpar_manager_tvsvd.cpp
    inv_analysis_matpar_manager_uniform.cpp
    inv_analysis_metropolis_kernel.cpp
    inv_analysis_objective_funct_disp.cpp
    inv_analysis_objective_funct_surfcurr.cpp
    inv_analysis_optimizer_base.cpp
    inv_analysis_optimizer_bruteforce.cpp
    inv_analysis_optimizer_factory.cpp
    inv_analysis_optimizer_lbfgs.cpp
    inv_analysis_optimizer_mh.cpp
    inv_analysis_optimizer_smc.cpp
    inv_analysis_optimizer_smc_predict.cpp
    inv_analysis_particle_comm.cpp
    inv_analysis_particle_data.cpp
    inv_analysis_particle_group.cpp
    inv_analysis_regularization_base.cpp
    inv_analysis_regularization_tikhonov.cpp
    inv_analysis_regularization_totalvariation.cpp
    inv_analysis_resulttest.cpp
    inv_analysis_timint_adjoint.cpp
    inv_analysis_timint_adjoint_prestress.cpp
    inv_analysis_utils.cpp
    inv_analysis_writer.cpp
    )

baci_add_library(alg_statinvanalysis ${OBJS_STATINV_ANALYSIS})

set(LIBS_ALG_STATINVANALYSIS
    # cmake-format: sortable
    inpar
    io
    lib
    )

baci_add_dependency(alg_statinvanalysis ${LIBS_ALG_STATINVANALYSIS})

set(OBJS_STATINV_ANALYSIS_LIB
    # cmake-format: sortable
    inv_analysis_particle_data.cpp
    )

baci_add_library(lib_statinv_ana ${OBJS_STATINV_ANALYSIS_LIB})

set(LIBS_LIB_STATINV_ANA
    # cmake-format: sortable
    lib
    )

baci_add_dependency(lib_statinv_ana ${LIBS_LIB_STATINV_ANA})
