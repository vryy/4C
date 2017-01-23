/*!---------------------------------------------------------------------
\file drt_globalproblem_enums.H

\brief Some important enums used in conjunction with DRT::Problem

\level 1

\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235

---------------------------------------------------------------------*/

#ifndef DRT_GLOBALPROBLEM_ENUMS_H
#define DRT_GLOBALPROBLEM_ENUMS_H

#include <string>
#include <map>

/*----------------------------------------------------------------------*
 | PROBLEM TYPES                                          m.gee 7/01    |
 *----------------------------------------------------------------------*/
typedef enum _PROBLEM_TYP
{
  // problem types in alphabetical order
  prb_none,                 // not a problem at all
  prb_ac_fsi,               // atherosclerosis growth/multiscale problem
  prb_acou,                 // acoustic problem
  prb_ale,                  // pure ale problem
  prb_art_net,              // arterial network problem _1D_ARTERY_
  prb_biofilm_fsi,          // biofilm growth problem
  prb_cardiac_monodomain,   // Cardiac electrophsiology problem
  prb_cavitation,           // cavitation problem
  prb_combust,              // combustion problem
  prb_ehl,                  // elastohydrodynamic lubrication problem (or lubrication structure interaction)
  prb_elch,                 // electrochemical problem
  prb_fluid,                // fluid problem
  prb_fluid_ale,            // fluid on an ale mesh (no structure)
  prb_fluid_redmodels,      // fluid_redairways problem
  prb_fluid_topopt,         // fluid topology optimization
  prb_fluid_xfem,           // fluid problem including XFEM interfaces
  prb_fluid_xfem_ls,        // xfluid calculations using levelset for cut.
  prb_fps3i,                // fluid porous structure scatra scatra interaction
  prb_fpsi,                 // fluid porous structure interaction problem
  prb_fpsi_xfem,            // fluid poro structure interaction problem including XFEM interfaces (atm just for FSI Interface!)
  prb_freesurf,             // free surface fluid
  prb_fsi,                  // fluid structure interaction problem
  prb_fsi_lung,             // airway fsi problem with attached parenchyma balloon
  prb_fsi_redmodels,        // fluid structure interaction problem
  prb_fsi_xfem,             // fluid structure interaction problem including XFEM interfaces
  prb_gas_fsi,              // fsi with gas transport
  prb_immersed_ale_fsi,     // immersed fsi coupled with ale fsi
  prb_immersed_cell,        // immersed cell migration process (poroelast-structure interaction)
  prb_immersed_fsi,         // immersed fsi
  prb_immersed_membrane_fsi,// immersed fsi with membrane finite elements
  prb_invana,               // inverse analysis of structures(so far)
  prb_level_set,            // level-set problem
  prb_loma,                 // low-Mach-number flow problem
  prb_lubrication,          // lubrication problem (reduced fluid model for elastohydrodynamic lubrication)
  prb_meshfree,             // particleMeshFree simulation
  prb_np_support,           // supporting procs for nested parallelism
  prb_opt,                  // structural optimization  problem
  prb_particle,             // particle simulation
  prb_polymernetwork,       // polymer network
  prb_poroelast,            // poroelasticity
  prb_poroscatra,           // passive scalar transport in porous media
  prb_porofluidmultiphase,  // multiphase flow in porous media
  prb_poromultiphase,       // multiphase flow in elastic porous media
  prb_poromultiphasescatra, // multiphase flow in elastic porous media with transport of species
  prb_red_airways,          // reduced dimensional airways
  prb_redairways_tissue,    // coupling of reduced-dimensional airways with 3D tissue model
  prb_scatra,               // scalar transport problem (e.g. convection-diffusion)
  prb_scatra_endoexocytosis,// biochemical transport problem (e.g. convection-diffusion) for endo-/and exocytosis
  prb_ssi,                  // scalar structure interaction
  prb_sti,                  // scalar-thermo interaction
  prb_struct_ale,           // structural problem, ale formulation
  prb_structure,            // structural problem
  prb_tfsi_aero,            // thermo fluid structure interaction problem including the aero code
  prb_thermo,               // thermal problem
  prb_thermo_fsi,           // thermo-fluid-structure-interaction problem
  prb_tsi,                  // thermal structure interaction
  prb_tutorial,             // student's c++/baci tutorial
  prb_two_phase_flow,       // two phase flow problem using the fluid and levelset modules
  prb_uq                    // uncertainty quantification
} PROBLEM_TYP;


/*----------------------------------------------------------------------*/
/* The known nptypes for nested parallelism                 ghamm 02/12 */
/*----------------------------------------------------------------------*/
typedef enum _NP_TYPE
{
  copy_dat_file,
  every_group_read_dat_file,
  separate_dat_files,
  no_nested_parallelism
} NP_TYPE;

#endif
