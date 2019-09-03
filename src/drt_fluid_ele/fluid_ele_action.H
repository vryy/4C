/*----------------------------------------------------------------------*/
/*!

\brief Enumeration of actions provided by the fluid element

\maintainer Martin Kronbichler

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FLUID_ELE_ACTION_H
#define FLUID_ELE_ACTION_H

namespace FLD
{
  /*--------------------------------------------------------------------------
   | enum that provides all possible fluid actions
   *--------------------------------------------------------------------------*/
  enum Action
  {
    none,
    calc_fluid_systemmat_and_residual,
    calc_porousflow_fluid_coupling,
    calc_poroscatra_mono_odblock,
    calc_loma_mono_odblock,
    calc_fluid_genalpha_sysmat_and_residual,
    calc_fluid_genalpha_update_for_subscales,
    calc_div_u,
    calc_mass_matrix,
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
    calc_vreman_const,
    integrate_shape,
    calc_divop,
    calc_mat_deriv_u_and_rot_u,
    calc_press_grad_and_div_eps,
    calc_volume_gaussint,
    calc_fluidfrac_projection,
    calc_node_normal,
    set_general_fluid_parameter,
    set_time_parameter,
    set_turbulence_parameter,
    set_loma_parameter,
    set_poro_parameter,
    set_general_face_fluid_parameter,
    set_general_face_xfem_parameter,
    set_general_fluid_xfem_parameter,
    set_topopt_parameter,
    calc_adjoint_systemmat_and_residual,
    //    calc_adjoint_neumann,
    set_general_adjoint_parameter,
    set_adjoint_time_parameter,
    calc_volume,
    interpolate_velocity_to_given_point_immersed,
    search_immersed_boundary_elements,
    least_squares_matrix_rhs_immersed_boundary,
    update_immersed_information,
    interpolate_velgrad_to_given_point,
    calc_artificial_velocity_divergence,
    interpolate_pressure_to_given_point,
    correct_immersed_fluid_bound_vel,
    interpolate_velocity_to_given_point,
    xwall_l2_projection,
    xwall_calc_mk,
    tauw_via_gradient,
    set_two_phase_parameter,
    project_fluid_field,
    interpolate_hdg_to_node,
    interpolate_hdg_for_hit,
    project_hdg_force_on_dof_vec_for_hit,
    project_hdg_initial_field_for_hit,
    velgradient_projection,
    presgradient_projection,
    calc_dt_via_cfl,
    calc_mass_flow_periodic_hill,
    reset_immersed_ele,
    calc_velgrad_ele_center,
    calc_pressure_average,
    update_local_solution
  };  // enum Action

  /*--------------------------------------------------------------------------
   | enum that provides all possible fluid actions on a boundary
   *--------------------------------------------------------------------------*/
  enum BoundaryAction
  {
    ba_none,
    integrate_Shapefunction,
    calc_area,
    calc_flowrate,
    flowratederiv,
    Outletimpedance,
    dQdu,
    ba_calc_node_normal,
    calc_node_curvature,
    calc_surface_tension,
    enforce_weak_dbc,
    estimate_Nitsche_trace_maxeigenvalue_,
    mixed_hybrid_dbc,
    flow_dep_pressure_bc,
    slip_supp_bc,
    navier_slip_bc,
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
    poro_splitnopenetration,
    poro_splitnopenetration_OD,
    poro_splitnopenetration_ODpres,
    poro_splitnopenetration_ODdisp,
    ba_calc_adjoint_neumann,
    fpsi_coupling
  };  // enum BoundaryAction

  /*--------------------------------------------------------------------------
   | enum that provides all possible fluid actions on a element interfaces
   *--------------------------------------------------------------------------*/
  enum IntFaceAction
  {
    ifa_none,
    EOS_and_GhostPenalty_stabilization
  };  // enum IntFaceAction

}  // namespace FLD

#endif
