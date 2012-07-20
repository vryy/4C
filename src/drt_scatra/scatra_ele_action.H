/*!---------------------------------------------------------------------------
\file scatra_ele_action.H

\brief provides enum of actions. Motivation was to have a quick enum
comparison instead of a cumberso me string compare.

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
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
  calc_mat_and_rhs,                // calc_condif_systemmat_and_residual,
  calc_subgrid_diffusivity_matrix, // calc_subgrid_diffusivity_matrix,
  calc_domain_and_bodyforce,       // calc_domain_and_bodyforce,
  get_material_parameters,         // get_material_parameters,
  integrate_shape_functions,       // integrate_shape_functions,
  time_update_material,            // time update for materials
  calc_initial_time_deriv,         // calc_initial_time_deriv,
  calc_flux_domain,                // calc_condif_flux,
  calc_mean_scalars,               // calc_mean_scalars,
  calc_error,                      // calc_error,
  calc_elch_conductivity,          // calc_elch_conductivity,
  calc_elch_initial_potential,     // calc_initial_potential_field,
  calc_TG_mat_and_rhs,             // levelset_TaylorGalerkin,
  reinitialize_levelset,           // reinitialize_levelset,
  calc_error_reinit                // calc_error_reinit,
}; // enum Action

/*--------------------------------------------------------------------------*/
/*!
 * \brief enum that provides all possible scatra actions on a boundary
 *///                                                             nis Mar12
/*--------------------------------------------------------------------------*/
enum BoundaryAction {
  bd_calc_normal_vectors,                 //calc_normal_vectors,
  bd_calc_elch_electrode_kinetics,        //calc_elch_electrode_kinetics,
  bd_calc_loma_therm_press,               //calc_therm_press,
  bd_integrate_shape_functions,           //integrate_shape_functions_boundary,
  bd_calc_Neumann_inflow,                 //calc_Neumann_inflow,
  bd_calc_convective_heat_transfer,       //calc_convective_heat_transfer,
  bd_calc_weak_Dirichlet,                 //WeakDirichlet,
  bd_calc_surface_permeability,           //calc_surface_permeability,
  bd_calc_TG_outflow,                     //levelset_TaylorGalerkin_boundary,
  bd_reinitialize_levelset,               //reinitialize_levelset_boundary,
  bd_add_convective_mass_flux             //add_convective_mass_flux
}; // enum BoundaryAction


} // namespace SCATRA

#endif // SCATRA_ELE_ACTION_H
