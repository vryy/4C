set(OBJS_DRT_STRUCTURE_NEW
    # cmake-format: sortable
    gauss_point_data_output_manager.cpp
    nox_nln_str_linearsystem.cpp
    str_dbc.cpp
    str_discretization_runtime_vtu_output_params.cpp
    str_factory.cpp
    str_functions.cpp
    str_impl_gemm.cpp
    str_impl_genalpha.cpp
    str_impl_genalpha_liegroup.cpp
    str_impl_generic.cpp
    str_impl_ost.cpp
    str_impl_prestress.cpp
    str_impl_statics.cpp
    str_integrator.cpp
    str_model_evaluator.cpp
    str_model_evaluator_cardiovascular0d.cpp
    str_model_evaluator_contact.cpp
    str_model_evaluator_contactdata.cpp
    str_model_evaluator_data.cpp
    str_model_evaluator_factory.cpp
    str_model_evaluator_generic.cpp
    str_model_evaluator_lagpenconstraint.cpp
    str_model_evaluator_meshtying.cpp
    str_model_evaluator_multiphysics.cpp
    str_model_evaluator_springdashpot.cpp
    str_model_evaluator_structure.cpp
    str_monitor_dbc.cpp
    str_nln_linearsystem_scaling.cpp
    str_nln_solver_factory.cpp
    str_nln_solver_fullnewton.cpp
    str_nln_solver_generic.cpp
    str_nln_solver_nox.cpp
    str_nln_solver_ptc.cpp
    str_nln_solver_utils.cpp
    str_predict_constdisvelaccpress.cpp
    str_predict_factory.cpp
    str_predict_generic.cpp
    str_predict_tangdis.cpp
    str_resulttest.cpp
    str_solver_factory.cpp
    str_timint_base.cpp
    str_timint_basedataglobalstate.cpp
    str_timint_basedataio.cpp
    str_timint_basedataio_monitor_dbc.cpp
    str_timint_basedataio_runtime_vtk_output.cpp
    str_timint_basedataio_runtime_vtp_output.cpp
    str_timint_basedatasdyn.cpp
    str_timint_explicit.cpp
    str_timint_factory.cpp
    str_timint_implicit.cpp
    str_timint_implicitbase.cpp
    str_timint_loca_continuation.cpp
    str_timint_locainterface.cpp
    str_timint_noxinterface.cpp
    str_utils.cpp
    )

baci_add_library(alg_structure_new ${OBJS_DRT_STRUCTURE_NEW})

set(LIBS_ALG_STRUCTURE_NEW
    # cmake-format: sortable
    drt_adapter
    drt_inpar
    drt_lib
    drt_lib_browniandyn
    drt_timestepping
    locacontinuation
    nlnsolver_nox
    )

target_link_libraries(alg_structure_new PUBLIC ${LIBS_ALG_STRUCTURE_NEW})
