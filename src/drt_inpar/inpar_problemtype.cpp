/*----------------------------------------------------------------------*/
/*!
\file inpar_problemtype.H

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/

#include "inpar_problemtype.H"
#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string,PROBLEM_TYP> DRT::StringToProblemTypeMap()
{
  static std::map<std::string,PROBLEM_TYP> string2prbtype;
  if (string2prbtype.size()==0)
  {
    string2prbtype["Structure"] =                                 prb_structure;
    string2prbtype["Structure_Ale"] =                             prb_struct_ale;
    string2prbtype["Fluid"] =                                     prb_fluid;
    string2prbtype["Fluid_XFEM"] =                                prb_fluid_xfem;
    string2prbtype["Fluid_Fluid_Ale"] =                           prb_fluid_fluid_ale;
    string2prbtype["Fluid_Fluid"] =                               prb_fluid_fluid;
    string2prbtype["Fluid_Fluid_FSI"] =                           prb_fluid_fluid_fsi;
    string2prbtype["Fluid_Ale"] =                                 prb_fluid_ale;
    string2prbtype["Fluid_RedModels"] =                          prb_fluid_redmodels;
    string2prbtype["Fluid_Freesurface"] =                         prb_freesurf;
    string2prbtype["Scalar_Transport"] =                          prb_scatra;
    string2prbtype["Fluid_Structure_Interaction"] =               prb_fsi;
    string2prbtype["Fluid_Structure_Interaction_RedModels"] =    prb_fsi_redmodels;
    string2prbtype["Fluid_Structure_Interaction_XFEM"] =          prb_fsi_xfem;
    string2prbtype["Ale"] =                                       prb_ale;
    string2prbtype["Thermo_Structure_Interaction"] =              prb_tsi;
    string2prbtype["Thermo"] =                                    prb_thermo;
    string2prbtype["Low_Mach_Number_Flow"] =                      prb_loma;
    string2prbtype["Electrochemistry"] =                          prb_elch;
    string2prbtype["Combustion"] =                                prb_combust;
    string2prbtype["ArterialNetwork"] =                           prb_art_net;
    string2prbtype["Fluid_Structure_Interaction_Lung"] =          prb_fsi_lung;
    string2prbtype["AeroCode_Thermo_Fluid_Structure_Interaction"] = prb_tfsi_aero;
    string2prbtype["ReducedDimensionalAirWays"] =                 prb_red_airways;
    string2prbtype["Gas_Fluid_Structure_Interaction"] =           prb_gas_fsi;
    string2prbtype["Biofilm_Fluid_Structure_Interaction"] =       prb_biofilm_fsi;
    string2prbtype["Thermo_Fluid_Structure_Interaction"] =        prb_thermo_fsi;
    string2prbtype["Fluid_Top_Opt"] =                             prb_fluid_topopt;
    string2prbtype["Poroelasticity"] =                            prb_poroelast;
    string2prbtype["Poroelastic_scalar_transport"] =              prb_poroscatra;
    string2prbtype["Fluid_Porous_Structure_Interaction"] =        prb_fpsi;
    string2prbtype["Structure_Scalar_Interaction"] =              prb_ssi;
    string2prbtype["NP_Supporting_Procs"] =                       prb_np_support;
    string2prbtype["RedAirways_Tissue"] =                         prb_redairways_tissue;
    string2prbtype["Particle"] =                                  prb_particle;
    string2prbtype["Cavitation"] =                                prb_cavitation;
    string2prbtype["Level_Set"] =                                 prb_level_set;
    string2prbtype["Fluid_Structure_Crack_Interaction"] =         prb_fsi_crack;
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
