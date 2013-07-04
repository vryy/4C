/*!---------------------------------------------------------------------------
\file fluid_ele_action.H

\brief provides enum of actions.

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*---------------------------------------------------------------------------*/

#ifndef FLUID_ELE_ACTION_H
#define FLUID_ELE_ACTION_H

namespace FLD
{

/*--------------------------------------------------------------------------
 | enum that provides all possible fluid actions
 *--------------------------------------------------------------------------*/
enum Action {
    none,
    calc_fluid_systemmat_and_residual,
    calc_porousflow_fluid_coupling,
    calc_loma_mono_odblock,
    calc_fluid_genalpha_sysmat_and_residual,
    calc_fluid_genalpha_update_for_subscales,
    calc_div_u,
    calc_fluid_error,
    calc_mean_Cai,
    set_mean_Cai,
    calc_turbulence_statistics,
    calc_loma_statistics,
    calc_turbscatra_statistics,
    calc_dissipation,
    calc_model_params_mfsubgr_scales,
    calc_fluid_box_filter,
    calc_smagorinsky_const,
    integrate_shape,
    calc_divop,
    calc_mat_derivative_u,
    calc_node_normal,
    set_general_fluid_parameter,
    set_time_parameter,
    set_turbulence_parameter,
    set_loma_parameter,
    set_poro_parameter,
    calc_adjoint_systemmat_and_residual,
//    calc_adjoint_neumann,
    set_general_adjoint_parameter,
    set_adjoint_time_parameter
}; // enum Action

/*--------------------------------------------------------------------------
 | enum that provides all possible fluid actions on a boundary
 *--------------------------------------------------------------------------*/
enum BoundaryAction {
    ba_none,
    integrate_Shapefunction,
    calc_area,
    calc_flowrate,
    flowratederiv,
    Outletimpedance,
    ba_calc_node_normal,
    calc_node_curvature,
    calc_surface_tension,
    enforce_weak_dbc,
    evaluate_nitsche_par,
    mixed_hybrid_dbc,
    flow_dep_pressure_bc,
    conservative_outflow_bc,
    calc_Neumann_inflow,
    calc_pressure_bou_int,
    center_of_mass_calc,
    traction_velocity_component,
    traction_Uv_integral_component,
    no_penetration,
    no_penetrationIDs,
    poro_boundary,
    poro_prescoupl,
    ba_calc_adjoint_neumann
}; // enum BoundaryAction

/*--------------------------------------------------------------------------
 | enum that provides all possible fluid actions on a element interfaces
 *--------------------------------------------------------------------------*/
enum IntFaceAction {
    ifa_none,
    EOS_and_GhostPenalty_stabilization
}; // enum IntFaceAction

}

#endif
