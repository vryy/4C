set(OBJS_STATINV_ANALYSIS
    # cmake-format: sortable
    chol_factor.cpp
    chol_factory.cpp
    DcsMatrix.cpp
    initial_guess.cpp
    invana_auglagr.cpp
    invana_base.cpp
    invana_cal_drt.cpp
    invana_control.cpp
    invana_factory.cpp
    invana_resulttest.cpp
    invana_utils.cpp
    invana_writer.cpp
    likelihood_evaluation.cpp
    matpar_manager.cpp
    matpar_manager_elementwise.cpp
    matpar_manager_patchwise.cpp
    matpar_manager_tvsvd.cpp
    matpar_manager_uniform.cpp
    metropolis_kernel.cpp
    objective_funct_disp.cpp
    objective_funct_surfcurr.cpp
    optimizer_base.cpp
    optimizer_bruteforce.cpp
    optimizer_factory.cpp
    optimizer_lbfgs.cpp
    optimizer_mh.cpp
    optimizer_smc.cpp
    optimizer_smc_predict.cpp
    particle_comm.cpp
    particle_data.cpp
    particle_group.cpp
    regularization_base.cpp
    regularization_tikhonov.cpp
    regularization_totalvariation.cpp
    timint_adjoint.cpp
    timint_adjoint_prestress.cpp
    )

baci_add_library(alg_statinvanalysis ${OBJS_STATINV_ANALYSIS})

set(LIBS_ALG_STATINVANALYSIS
    # cmake-format: sortable
    drt_inpar
    drt_io
    drt_lib
    )

target_link_libraries(alg_statinvanalysis PUBLIC ${LIBS_ALG_STATINVANALYSIS})

set(OBJS_STATINV_ANALYSIS_LIB
    # cmake-format: sortable
    particle_data.cpp
    )

baci_add_library(drt_lib_statinv_ana ${OBJS_STATINV_ANALYSIS_LIB})

set(LIBS_DRT_LIB_STATINV_ANA
    # cmake-format: sortable
    drt_lib
    )

target_link_libraries(drt_lib_statinv_ana PUBLIC ${LIBS_DRT_LIB_STATINV_ANA})
