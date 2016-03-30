/*!---------------------------------------------------------------------------
\file scatra_ele_action.H

\brief provides enum of actions. Motivation was to have a quick enum
comparison instead of a cumbersome std::string compare.

<pre>
\maintainer Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
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
    calc_cell_mechanotransduction,               // calculate response of cell internal transport to mechanical stimulus
    calc_dissipation,                            // calculate dissipation introduced by stabilization and turbulence models
    calc_domain_and_bodyforce,                   // calc_domain_and_bodyforce,
    calc_domain_integral,                        // n/a
    calc_elch_boundary_kinetics_point,           // evaluate electrode kinetics point conditions
    calc_elch_conductivity,                      // calc_elch_conductivity,
    calc_elch_domain_kinetics,                   // calculate electrode kinetics in domain
    calc_elch_electrode_soc,                     // calculate electrode state of charge
    calc_error,                                  // calc_error,
    calc_error_reinit,                           // calc_error_reinit,
    calc_flux_domain,                            // calc_condif_flux,
    calc_immersed_element_source,                // calculate sink/source term due to segregation from an immersed discretizaton
    calc_immersed_phi_at_given_point,            // calculate concentration at a given point during immersed cell migration
    calc_initial_time_deriv,                     // calc_initial_time_deriv,
    calc_integr_grad_diff,                       // calculate diffusion coefficient gradient for photoacoustics
    calc_integr_grad_reac,                       // calculate reaction coefficient gradient for photoacoustics
    calc_integr_pat_rhsvec,                      // integrate rhs vec for adjoint pat
    calc_integr_reaction,                        // calculate integral reaction
    calc_mass_center_smoothingfunct,             // calculate mass center of minus domain, using the smoothing function. (Interface thickness required)
    calc_mat_and_rhs,                            // calc_condif_systemmat_and_residual,
    calc_mat_and_rhs_lsreinit_correction_step,   // level-set action: setup of system for correction step
    calc_mean_Cai,                               // calc_mean_Cai (multifractal subgrid-scales)
    calc_total_and_mean_scalars,                 // calculate total and mean values of transported scalars
    calc_node_based_reinit_velocity,             // level-set action: compute node-based velocity field for reinitialization equation
    calc_scatra_box_filter,                      // calc_scatra_box_filter,
    calc_scatra_mono_odblock_fluid,              // calculate off diagonal matrix block for coupling with fluid
    calc_scatra_mono_odblock_mesh,               // calculate off diagonal matrix block for coupling with mesh (shape derivatives)
    calc_scatra_mono_odblock_scatrathermo,       // calculate off-diagonal matrix block (derivatives of scatra residuals w.r.t. thermo dofs) for scatra-thermo interaction
    calc_scatra_mono_odblock_thermoscatra,       // calculate off-diagonal matrix block (derivatives of thermo residuals w.r.t. scatra dofs) for scatra-thermo interaction
    calc_subgrid_diffusivity_matrix,             // calc_subgrid_diffusivity_matrix,
    calc_TG_mat_and_rhs,                         // levelset_TaylorGalerkin,
    calc_turbulent_prandtl_number,               // calc_turbulent_prandtl_number
    calc_vreman_scatra,                          // calc_vreman_scatra,
    check_scatra_element_parameter,              // set parameter defined in every single element
    get_material_internal_state,                 // get material internal state (for postprocessing for example),
    get_material_ionic_currents,                 // get material ionic currents (for postprocessing for example),
    get_material_parameters,                     // get_material_parameters,
    integrate_shape_functions,                   // integrate_shape_functions,
    recon_curvature_at_nodes,                    // reconstruct curvature at node points through L2-projection
    recon_gradients_at_nodes,                    // reconstruct gradients at node points through L2-projection
    reinitialize_levelset,                       // reinitialize_levelset,
    set_diffcond_scatra_parameter,               // set diffusion-conduction parameters for element evaluation
    set_elch_scatra_parameter,                   // set elch parameters for element evaluation
    set_general_scatra_parameter,                // set general parameters for element evaluation
    set_lsreinit_scatra_parameter,               // set reinitialization parameters for element evaluation
    set_material_internal_state,                 // set material internal state (for restart for example),
    set_mean_Cai,                                // special action for mfs
    set_time_parameter,                          // set time-integration parameters for element evaluation
    set_turbulence_scatra_parameter,             // set turbulence parameters for element evaluation
    time_update_material                         // time update for materials
  }; // enum Action


  /*--------------------------------------------------------------------------*/
  /*!
   * \brief enum that provides all possible scatra actions on a boundary
   *///                                                             nis Mar12
  /*--------------------------------------------------------------------------*/
  enum BoundaryAction {
    bd_add_convective_mass_flux,
    bd_calc_boundary_integral,
    bd_calc_convective_heat_transfer,
    bd_calc_elch_boundary_kinetics,
    bd_calc_elch_cell_voltage,
    bd_calc_elch_linearize_nernst,
    bd_calc_fps3i_surface_permeability,   // Kedem-Katchalsky equations --> FPS3I problem
    bd_calc_fs3i_surface_permeability,
    bd_calc_loma_therm_press,
    bd_calc_Neumann,
    bd_calc_Neumann_inflow,
    bd_calc_mechanotransduction,
    bd_calc_normal_vectors,
    bd_calc_Robin,                        // Robin boundary condition
    bd_calc_s2icoupling,                  // evaluate scatra-scatra interface coupling conditions
    bd_calc_s2icoupling_od,               // evaluate derivatives of scatra-scatra interface coupling conditions w.r.t. secondary dofs
    bd_calc_weak_Dirichlet,
    bd_integrate_shape_functions
  }; // enum SCATRA::BoundaryAction
} // namespace SCATRA
#endif // SCATRA_ELE_ACTION_H
