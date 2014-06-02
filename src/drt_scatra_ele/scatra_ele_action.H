/*!---------------------------------------------------------------------------
\file scatra_ele_action.H

\brief provides enum of actions. Motivation was to have a quick enum
comparison instead of a cumbersome std::string compare.

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*---------------------------------------------------------------------------*/

#ifndef SCATRA_ELE_ACTION_H
#define SCATRA_ELE_ACTION_H

namespace SCATRA
{

/*--------------------------------------------------------------------------*/
/*!
 * \brief enum that provides all possible scatra actions
 *///                                                             nis Mar12
/*--------------------------------------------------------------------------*/
enum Action {
  // domain action
  set_general_scatra_parameter,     // set general parameters for element evaluation
  set_turbulence_scatra_parameter,  // set turbulence parameters for element evaluation
  set_time_parameter,               // set time-integration parameters for element evaluation
  set_elch_scatra_parameter,        // set elch parameters for element evaluation
  set_lsreinit_scatra_parameter,    // set reinitialization parameters for element evaluation
  calc_mat_and_rhs,                // calc_condif_systemmat_and_residual,
  calc_subgrid_diffusivity_matrix, // calc_subgrid_diffusivity_matrix,
  calc_domain_and_bodyforce,       // calc_domain_and_bodyforce,
  get_material_parameters,         // get_material_parameters,
  get_material_internal_state,     // get material internal state (for postprocessing for example),
  set_material_internal_state,     // set material internal state (for restart for example),
  get_material_ionic_currents,     // get material ionic currents (for postprocessing for example),
  integrate_shape_functions,       // integrate_shape_functions,
  time_update_material,            // time update for materials
  check_scatra_element_parameter,  // set parameter defined in every single element
  calc_initial_time_deriv,         // calc_initial_time_deriv,
  calc_flux_domain,                // calc_condif_flux,
  calc_mean_scalars,               // calc_mean_scalars,
  calc_error,                      // calc_error,
  calc_elch_conductivity,          // calc_elch_conductivity,
  calc_elch_initial_potential,     // calc_initial_potential_field,
  calc_scatra_box_filter,          // calc_scatra_box_filter,
  calc_turbulent_prandtl_number,   // calc_turbulent_prandtl_number
  calc_vreman_scatra,              // calc_vreman_scatra,
  calc_TG_mat_and_rhs,             // levelset_TaylorGalerkin,
  reinitialize_levelset,           // reinitialize_levelset,
  calc_error_reinit,               // calc_error_reinit,
  calc_mean_Cai,                   // calc_mean_Cai (multifractal subgrid-scales)
  calc_dissipation,                // calculate dissipation introduced by stabilization and turbulence models
  calc_integr_reaction,            // calculate integral reaction
  calc_elch_electrode_kinetics,    // calc_elch_electrode_kinetics,
  set_mean_Cai,                    // special action for mfs
  calc_mat_and_rhs_lsreinit_correction_step, // level-set action: setup of system for correction step
  calc_node_based_reinit_velocity, // level-set action: compute node-based velocity field for reinitialization equation
  calc_integr_grad_reac,           // calculate regularization part in gradient for adjoint photoacoustics
  calc_integr_pat_rhsvec           // integrate rhs vec for adjoint pat
}; // enum Action

/*--------------------------------------------------------------------------*/
/*!
 * \brief enum that provides all possible scatra actions on a boundary
 *///                                                             nis Mar12
/*--------------------------------------------------------------------------*/
enum BoundaryAction {
// new action                       //PRE_BACI_CLEANING NAME:
  bd_integrate_shape_functions,     //integrate_shape_functions,
  bd_calc_normal_vectors,           //calc_normal_vectors,
  bd_calc_weak_Dirichlet,           //WeakDirichlet,
  bd_calc_Neumann_inflow,           //calc_Neumann_inflow,
  bd_add_convective_mass_flux,      //add_convective_mass_flux
  bd_calc_convective_heat_transfer, //calc_convective_heat_transfer,
  bd_calc_elch_electrode_kinetics,  //calc_elch_electrode_kinetics,
  bd_calc_elch_linearize_nernst,    //calc_elch_electrode_kinetics,
  bd_calc_loma_therm_press,         //calc_therm_press,
  bd_calc_surface_permeability,     //calc_surface_permeability,
  bd_calc_fpssi_conditions

}; // enum SCATRA::BoundaryAction
} // namespace SCATRA

#endif // SCATRA_ELE_ACTION_H
