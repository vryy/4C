/*----------------------------------------------------------------------*/
/*! \file

\brief Enumeration of actions provided by the fluid element


\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FLUID_ELE_ACTION_HPP
#define FOUR_C_FLUID_ELE_ACTION_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace FLD
{
  /*--------------------------------------------------------------------------
   | enum that provides all possible fluid actions
   *--------------------------------------------------------------------------*/
  enum Action
  {
    calc_dissipation,
    calc_div_u,
    calc_divop,
    calc_dt_via_cfl,
    calc_fluid_box_filter,
    calc_fluid_error,
    calc_fluid_genalpha_sysmat_and_residual,
    calc_fluid_genalpha_update_for_subscales,
    calc_fluid_systemmat_and_residual,
    calc_loma_mono_odblock,
    calc_loma_statistics,
    calc_mass_flow_periodic_hill,
    calc_mass_matrix,
    calc_mean_Cai,
    calc_model_params_mfsubgr_scales,
    calc_node_normal,
    calc_poroscatra_mono_odblock,
    calc_porousflow_fluid_coupling,
    calc_pressure_average,
    calc_smagorinsky_const,
    calc_turbscatra_statistics,
    calc_turbulence_statistics,
    calc_velgrad_ele_center,
    calc_volume,
    calc_vreman_const,
    correct_immersed_fluid_bound_vel,
    integrate_shape,
    interpolate_hdg_for_hit,
    interpolate_hdg_to_node,
    interpolate_pressure_to_given_point,
    interpolate_velgrad_to_given_point,
    interpolate_velocity_to_given_point,
    interpolate_velocity_to_given_point_immersed,
    none,
    presgradient_projection,
    project_fluid_field,
    project_hdg_force_on_dof_vec_for_hit,
    project_hdg_initial_field_for_hit,
    reset_immersed_ele,
    set_general_face_fluid_parameter,
    set_general_face_xfem_parameter,
    set_general_fluid_parameter,
    set_general_fluid_xfem_parameter,
    set_loma_parameter,
    set_mean_Cai,
    set_poro_parameter,
    set_time_parameter,
    set_turbulence_parameter,
    tauw_via_gradient,
    update_immersed_information,
    update_local_solution,
    velgradient_projection,
    xwall_calc_mk,
    xwall_l2_projection,
  };  // enum Action

  /*--------------------------------------------------------------------------
   | enum that provides all possible fluid actions on a boundary
   *--------------------------------------------------------------------------*/
  enum BoundaryAction
  {
    Outletimpedance,
    ba_calc_node_normal,
    ba_none,
    calc_Neumann_inflow,
    calc_area,
    calc_flowrate,
    calc_node_curvature,
    calc_pressure_bou_int,
    calc_surface_tension,
    center_of_mass_calc,
    dQdu,
    enforce_weak_dbc,
    estimate_Nitsche_trace_maxeigenvalue_,
    flow_dep_pressure_bc,
    flowratederiv,
    fpsi_coupling,
    integrate_Shapefunction,
    mixed_hybrid_dbc,
    navier_slip_bc,
    no_penetration,
    no_penetrationIDs,
    poro_boundary,
    poro_prescoupl,
    poro_splitnopenetration,
    poro_splitnopenetration_OD,
    poro_splitnopenetration_ODdisp,
    poro_splitnopenetration_ODpres,
    slip_supp_bc,
    traction_Uv_integral_component,
    traction_velocity_component,
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

FOUR_C_NAMESPACE_CLOSE

#endif
