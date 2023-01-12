set(OBJS_FSI
    # cmake-format: sortable
    fs_monolithic.cpp
    fsi_algorithm.cpp
    fsi_constr_overlapprec.cpp
    fsi_constrmonolithic.cpp
    fsi_constrmonolithic_fluidsplit.cpp
    fsi_constrmonolithic_structuresplit.cpp
    fsi_debugwriter.cpp
    fsi_dirichletneumann.cpp
    fsi_dirichletneumann_disp.cpp
    fsi_dirichletneumann_factory.cpp
    fsi_dirichletneumann_vel.cpp
    fsi_dirichletneumann_volcoupl.cpp
    fsi_dirichletneumannslideale.cpp
    fsi_dyn.cpp
    fsi_fluid_ale.cpp
    fsi_fluidfluidmonolithic_fluidsplit.cpp
    fsi_fluidfluidmonolithic_fluidsplit_nonox.cpp
    fsi_fluidfluidmonolithic_structuresplit.cpp
    fsi_fluidfluidmonolithic_structuresplit_nonox.cpp
    fsi_lung_overlapprec.cpp
    fsi_lungmonolithic.cpp
    fsi_lungmonolithic_fluidsplit.cpp
    fsi_lungmonolithic_structuresplit.cpp
    fsi_monolithic.cpp
    fsi_monolithic_linearsystem.cpp
    fsi_monolithic_nonox.cpp
    fsi_monolithic_redistribution.cpp
    fsi_monolithic_timintada.cpp
    fsi_monolithicfluidsplit.cpp
    fsi_monolithicstructuresplit.cpp
    fsi_mortarmonolithic_fluidsplit.cpp
    fsi_mortarmonolithic_structuresplit.cpp
    fsi_nox_aitken.cpp
    fsi_nox_epsilon.cpp
    fsi_nox_fixpoint.cpp
    fsi_nox_group.cpp
    fsi_nox_jacobian.cpp
    fsi_nox_linearsystem_gcr.cpp
    fsi_nox_mpe.cpp
    fsi_nox_newton.cpp
    fsi_nox_sd.cpp
    fsi_noxlinsys.cpp
    fsi_overlapprec.cpp
    fsi_overlapprec_amgnxn.cpp
    fsi_overlapprec_fsiamg.cpp
    fsi_overlapprec_fsiamg_analysis.cpp
    fsi_overlapprec_fsiamg_smoothers.cpp
    fsi_overlapprec_fsiamg_vcycles.cpp
    fsi_overlapprec_hybrid.cpp
    fsi_partitioned.cpp
    fsi_resulttest.cpp
    fsi_slidingmonolithic_fluidsplit.cpp
    fsi_slidingmonolithic_structuresplit.cpp
    fsi_statustest.cpp
    fsi_str_model_evaluator_partitioned.cpp
    fsi_structureale.cpp
    fsi_utils.cpp
    )

baci_add_library(alg_fsi ${OBJS_FSI})

set(LIBS_ALG_FSI
    # cmake-format: sortable
    adapter
    alg_fbi
    alg_fsi_xfem
    linalg
    solver
    )

target_link_libraries(alg_fsi PUBLIC ${LIBS_ALG_FSI})
