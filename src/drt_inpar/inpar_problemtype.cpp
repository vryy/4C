/*----------------------------------------------------------------------*/
/*!
 *
\file inpar_problemtype.cpp
\brief convert problem type string to enum
\level 1
\maintainer Martin Kronbichler
*/

/*----------------------------------------------------------------------*/

#include "inpar_problemtype.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string,PROBLEM_TYP> DRT::StringToProblemTypeMap()
{
  static std::map<std::string,PROBLEM_TYP> string2prbtype;

  if(string2prbtype.size() == 0)
  {
    // problem types in alphabetical order
    string2prbtype["Acoustics"]                                        = prb_acou;
    string2prbtype["Ale"]                                              = prb_ale;
    string2prbtype["ArterialNetwork"]                                  = prb_art_net;
    string2prbtype["Atherosclerosis_Fluid_Structure_Interaction"]      = prb_ac_fsi;
    string2prbtype["Biofilm_Fluid_Structure_Interaction"]              = prb_biofilm_fsi;
    string2prbtype["Cardiac_Monodomain"]                               = prb_cardiac_monodomain;
    string2prbtype["Cavitation"]                                       = prb_cavitation;
    string2prbtype["Elastohydrodynamic_Lubrication"]                   = prb_ehl;
    string2prbtype["Electrochemistry"]                                 = prb_elch;
    string2prbtype["Electromagnetics"]                                 = prb_elemag;
    string2prbtype["Fluid"]                                            = prb_fluid;
    string2prbtype["Fluid_Ale"]                                        = prb_fluid_ale;
    string2prbtype["Fluid_Freesurface"]                                = prb_freesurf;
    string2prbtype["Fluid_Poro_Structure_Interaction_XFEM"]            = prb_fpsi_xfem;
    string2prbtype["Fluid_Porous_Structure_Interaction"]               = prb_fpsi;
    string2prbtype["Fluid_Porous_Structure_Scalar_Scalar_Interaction"] = prb_fps3i;
    string2prbtype["Fluid_RedModels"]                                  = prb_fluid_redmodels;
    string2prbtype["Fluid_Structure_Interaction"]                      = prb_fsi;
    string2prbtype["Fluid_Structure_Interaction_Lung"]                 = prb_fsi_lung;
    string2prbtype["Fluid_Structure_Interaction_RedModels"]            = prb_fsi_redmodels;
    string2prbtype["Fluid_Structure_Interaction_XFEM"]                 = prb_fsi_xfem;
    string2prbtype["Fluid_Top_Opt"]                                    = prb_fluid_topopt;
    string2prbtype["Fluid_XFEM"]                                       = prb_fluid_xfem;
    string2prbtype["Fluid_XFEM_LevelSet"]                              = prb_fluid_xfem_ls;
    string2prbtype["Gas_Fluid_Structure_Interaction"]                  = prb_gas_fsi;
    string2prbtype["Immersed_ALE_FSI"]                                 = prb_immersed_ale_fsi;
    string2prbtype["Immersed_CellMigration"]                           = prb_immersed_cell;
    string2prbtype["Immersed_FSI"]                                     = prb_immersed_fsi;
    string2prbtype["Immersed_Membrane_FSI"]                            = prb_immersed_membrane_fsi;
    string2prbtype["Inverse_Analysis"]                                 = prb_invana;
    string2prbtype["Level_Set"]                                        = prb_level_set;
    string2prbtype["Low_Mach_Number_Flow"]                             = prb_loma;
    string2prbtype["Lubrication"]                                      = prb_lubrication;
    string2prbtype["NP_Supporting_Procs"]                              = prb_np_support;
    string2prbtype["Particle"]                                         = prb_particle;
    string2prbtype["Particle_Structure_Interaction"]                   = prb_pasi;
    string2prbtype["Polymer_Network"]                                  = prb_polymernetwork;
    string2prbtype["Poroelastic_scalar_transport"]                     = prb_poroscatra;
    string2prbtype["Poroelasticity"]                                   = prb_poroelast;
    string2prbtype["Multiphase_Poroelasticity"]                        = prb_poromultiphase;
    string2prbtype["Multiphase_Poroelasticity_ScaTra"]                 = prb_poromultiphasescatra;
    string2prbtype["Multiphase_Porous_Flow"]                           = prb_porofluidmultiphase;
    string2prbtype["RedAirways_Tissue"]                                = prb_redairways_tissue;
    string2prbtype["ReducedDimensionalAirWays"]                        = prb_red_airways;
    string2prbtype["Scalar_Thermo_Interaction"]                        = prb_sti;
    string2prbtype["Scalar_Transport"]                                 = prb_scatra;
    string2prbtype["Scalar_Transport_EndoExocytosis"]                  = prb_scatra_endoexocytosis;
    string2prbtype["Structure"]                                        = prb_structure;
    string2prbtype["Structure_Ale"]                                    = prb_struct_ale;
    string2prbtype["Structure_Scalar_Interaction"]                     = prb_ssi;
    string2prbtype["Thermo"]                                           = prb_thermo;
    string2prbtype["Thermo_Structure_Interaction"]                     = prb_tsi;
    string2prbtype["Thermo_Fluid_Structure_Interaction"]               = prb_thermo_fsi;
    string2prbtype["Tutorial"]                                         = prb_tutorial;
    string2prbtype["Two_Phase_Flow"]                                   = prb_two_phase_flow;
    string2prbtype["UQ"]                                               = prb_uq;
    string2prbtype["Variational_Chemical_Diffusion"]                   = prb_var_chemdiff;
    string2prbtype["XContact"]                                         = prb_xcontact;
  }

  return string2prbtype;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PROBLEM_TYP DRT::StringToProblemType(std::string name)
{
  std::map<std::string,PROBLEM_TYP> map = StringToProblemTypeMap();
  std::map<std::string,PROBLEM_TYP>::const_iterator i = map.find(name);
  if (i!=map.end())
    return i->second;
  dserror("unsupported problem name '%s'",name.c_str());

  return prb_none;
}
