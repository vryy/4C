set(OBJS_STATINV_ANALYSIS
    invana_cal_drt.cpp
    invana_control.cpp
    invana_writer.cpp
    timint_adjoint.cpp
    timint_adjoint_prestress.cpp
    objective_funct_disp.cpp
    objective_funct_surfcurr.cpp
    matpar_manager.cpp
    matpar_manager_uniform.cpp
    matpar_manager_elementwise.cpp
    matpar_manager_patchwise.cpp
    matpar_manager_tvsvd.cpp
    regularization_base.cpp
    regularization_totalvariation.cpp
    regularization_tikhonov.cpp
    invana_resulttest.cpp
    invana_utils.cpp
    invana_factory.cpp
    invana_base.cpp
    invana_auglagr.cpp
    initial_guess.cpp
    optimizer_factory.cpp
    optimizer_base.cpp
    optimizer_lbfgs.cpp
    optimizer_smc.cpp
    optimizer_smc_predict.cpp
    optimizer_mh.cpp
    optimizer_bruteforce.cpp
    particle_group.cpp
    particle_data.cpp
    particle_comm.cpp
    likelihood_evaluation.cpp
    DcsMatrix.cpp
    chol_factory.cpp
    chol_factor.cpp
    metropolis_kernel.cpp
    )

baci_add_library(alg_statinvanalysis ${OBJS_STATINV_ANALYSIS})
target_link_libraries(alg_statinvanalysis PUBLIC drt_lib drt_inpar drt_io)

set(OBJS_STATINV_ANALYSIS_LIB particle_data.cpp)
baci_add_library(drt_lib_statinv_ana ${OBJS_STATINV_ANALYSIS_LIB})
target_link_libraries(drt_lib_statinv_ana PUBLIC drt_lib)
