set(OBJS_STRUCTURE_NEW
    # cmake-format: sortable
    structure_new_dbc.cpp
    structure_new_discretization_runtime_vtu_output_params.cpp
    structure_new_factory.cpp
    structure_new_functions.cpp
    structure_new_gauss_point_data_output_manager.cpp
    structure_new_impl_gemm.cpp
    structure_new_impl_genalpha.cpp
    structure_new_impl_genalpha_liegroup.cpp
    structure_new_impl_generic.cpp
    structure_new_impl_ost.cpp
    structure_new_impl_prestress.cpp
    structure_new_impl_statics.cpp
    structure_new_integrator.cpp
    structure_new_model_evaluator.cpp
    structure_new_model_evaluator_contact.cpp
    structure_new_model_evaluator_contactdata.cpp
    structure_new_model_evaluator_data.cpp
    structure_new_model_evaluator_factory.cpp
    structure_new_model_evaluator_generic.cpp
    structure_new_model_evaluator_lagpenconstraint.cpp
    structure_new_model_evaluator_meshtying.cpp
    structure_new_model_evaluator_multiphysics.cpp
    structure_new_model_evaluator_springdashpot.cpp
    structure_new_model_evaluator_structure.cpp
    structure_new_monitor_dbc.cpp
    structure_new_nln_linearsystem_scaling.cpp
    structure_new_nln_solver_factory.cpp
    structure_new_nln_solver_fullnewton.cpp
    structure_new_nln_solver_generic.cpp
    structure_new_nln_solver_nox.cpp
    structure_new_nln_solver_ptc.cpp
    structure_new_nln_solver_utils.cpp
    structure_new_nox_nln_str_linearsystem.cpp
    structure_new_predict_constdisvelaccpress.cpp
    structure_new_predict_factory.cpp
    structure_new_predict_generic.cpp
    structure_new_predict_tangdis.cpp
    structure_new_resulttest.cpp
    structure_new_solver_factory.cpp
    structure_new_timint_base.cpp
    structure_new_timint_basedataglobalstate.cpp
    structure_new_timint_basedataio.cpp
    structure_new_timint_basedataio_monitor_dbc.cpp
    structure_new_timint_basedataio_runtime_vtk_output.cpp
    structure_new_timint_basedataio_runtime_vtp_output.cpp
    structure_new_timint_basedatasdyn.cpp
    structure_new_timint_explicit.cpp
    structure_new_timint_factory.cpp
    structure_new_timint_implicit.cpp
    structure_new_timint_implicitbase.cpp
    structure_new_timint_loca_continuation.cpp
    structure_new_timint_locainterface.cpp
    structure_new_timint_noxinterface.cpp
    structure_new_utils.cpp
    )

baci_add_library(alg_structure_new ${OBJS_STRUCTURE_NEW})

set(LIBS_ALG_STRUCTURE_NEW
    # cmake-format: sortable
    adapter
    inpar
    lib
    lib_browniandyn
    locacontinuation
    nlnsolver_nox
    timestepping
    )

baci_add_dependency(alg_structure_new ${LIBS_ALG_STRUCTURE_NEW})
