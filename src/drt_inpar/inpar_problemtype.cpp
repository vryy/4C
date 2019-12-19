/*----------------------------------------------------------------------*/
/*! \file
\brief convert problem type string to enum
\level 1
\maintainer Martin Kronbichler
*/

/*----------------------------------------------------------------------*/

#include "inpar_problemtype.H"

#include "drt_validparameters.H"

#include "../drt_lib/drt_dserror.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void INPAR::PROBLEMTYPE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& type = list->sublist("PROBLEM TYP", false, "");

  {
    Teuchos::Array<std::string> name;
    Teuchos::Array<int> label;
    Teuchos::Array<std::string> shape_name;
    Teuchos::Array<int> shape_label;

    // fill the arrays
    {
      std::map<std::string, ProblemType> map = StringToProblemTypeMap();
      std::map<std::string, ProblemType>::const_iterator i;
      for (i = map.begin(); i != map.end(); ++i)
      {
        name.push_back(i->first);
        label.push_back(i->second);
      }
      std::map<std::string, ShapeFunctionType> shape_map = StringToShapeFunctionTypeMap();
      std::map<std::string, ShapeFunctionType>::const_iterator j;
      for (j = shape_map.begin(); j != shape_map.end(); ++j)
      {
        shape_name.push_back(j->first);
        shape_label.push_back((int)j->second);
      }
    }

    setStringToIntegralParameter<int>(
        "PROBLEMTYP", "Fluid_Structure_Interaction", "", name, label, &type);

    setStringToIntegralParameter<int>("SHAPEFCT", "Polynomial",
        "Defines the function spaces for the spatial approximation", shape_name, shape_label,
        &type);
  }

  IntParameter("RESTART", 0, "", &type);
  DoubleParameter("RESTARTTIME", -1.0, "Used defined restart time", &type);
  IntParameter("RANDSEED", -1, "Set the random seed. If < 0 use current time.", &type);

#if 0  // currently not in use
//  BoolParameter("BANDWITHOPT","No","Do bandwith optimization of dof numbering",&type);
  setStringToIntegralParameter<int>("BANDWIDTHOPT","No",
                                    "Do bandwith optimization of dof numbering",
                                    yesnotuple,yesnovalue,&type);
#endif
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, ProblemType> INPAR::PROBLEMTYPE::StringToProblemTypeMap()
{
  static std::map<std::string, ProblemType> string2prbtype;

  if (string2prbtype.size() == 0)
  {
    // problem types in alphabetical order
    string2prbtype["Acoustics"] = prb_acou;
    string2prbtype["Ale"] = prb_ale;
    string2prbtype["ArterialNetwork"] = prb_art_net;
    string2prbtype["Atherosclerosis_Fluid_Structure_Interaction"] = prb_ac_fsi;
    string2prbtype["Biofilm_Fluid_Structure_Interaction"] = prb_biofilm_fsi;
    string2prbtype["Cardiac_Monodomain"] = prb_cardiac_monodomain;
    string2prbtype["Elastohydrodynamic_Lubrication"] = prb_ehl;
    string2prbtype["Electrochemistry"] = prb_elch;
    string2prbtype["Electromagnetics"] = prb_elemag;
    string2prbtype["Fluid"] = prb_fluid;
    string2prbtype["Fluid_Ale"] = prb_fluid_ale;
    string2prbtype["Fluid_Beam_Interaction"] = prb_fbi;
    string2prbtype["Fluid_Freesurface"] = prb_freesurf;
    string2prbtype["Fluid_Poro_Structure_Interaction_XFEM"] = prb_fpsi_xfem;
    string2prbtype["Fluid_Porous_Structure_Interaction"] = prb_fpsi;
    string2prbtype["Fluid_Porous_Structure_Scalar_Scalar_Interaction"] = prb_fps3i;
    string2prbtype["Fluid_RedModels"] = prb_fluid_redmodels;
    string2prbtype["Fluid_Structure_Interaction"] = prb_fsi;
    string2prbtype["Fluid_Structure_Interaction_Lung"] = prb_fsi_lung;
    string2prbtype["Fluid_Structure_Interaction_RedModels"] = prb_fsi_redmodels;
    string2prbtype["Fluid_Structure_Interaction_XFEM"] = prb_fsi_xfem;
    string2prbtype["Fluid_Top_Opt"] = prb_fluid_topopt;
    string2prbtype["Fluid_XFEM"] = prb_fluid_xfem;
    string2prbtype["Fluid_XFEM_LevelSet"] = prb_fluid_xfem_ls;
    string2prbtype["Gas_Fluid_Structure_Interaction"] = prb_gas_fsi;
    string2prbtype["Immersed_ALE_FSI"] = prb_immersed_ale_fsi;
    string2prbtype["Immersed_CellMigration"] = prb_immersed_cell;
    string2prbtype["Immersed_FSI"] = prb_immersed_fsi;
    string2prbtype["Immersed_Membrane_FSI"] = prb_immersed_membrane_fsi;
    string2prbtype["Inverse_Analysis"] = prb_invana;
    string2prbtype["Level_Set"] = prb_level_set;
    string2prbtype["Low_Mach_Number_Flow"] = prb_loma;
    string2prbtype["Lubrication"] = prb_lubrication;
    string2prbtype["NP_Supporting_Procs"] = prb_np_support;
    string2prbtype["Particle"] = prb_particle;
    string2prbtype["Particle_Structure_Interaction"] = prb_pasi;
    string2prbtype["Polymer_Network"] = prb_polymernetwork;
    string2prbtype["Poroelastic_scalar_transport"] = prb_poroscatra;
    string2prbtype["Poroelasticity"] = prb_poroelast;
    string2prbtype["Multiphase_Poroelasticity"] = prb_poromultiphase;
    string2prbtype["Multiphase_Poroelasticity_ScaTra"] = prb_poromultiphasescatra;
    string2prbtype["Multiphase_Porous_Flow"] = prb_porofluidmultiphase;
    string2prbtype["RedAirways_Tissue"] = prb_redairways_tissue;
    string2prbtype["ReducedDimensionalAirWays"] = prb_red_airways;
    string2prbtype["Scalar_Thermo_Interaction"] = prb_sti;
    string2prbtype["Scalar_Transport"] = prb_scatra;
    string2prbtype["Scalar_Transport_EndoExocytosis"] = prb_scatra_endoexocytosis;
    string2prbtype["Structure"] = prb_structure;
    string2prbtype["Structure_Ale"] = prb_struct_ale;
    string2prbtype["Structure_Scalar_Interaction"] = prb_ssi;
    string2prbtype["Thermo"] = prb_thermo;
    string2prbtype["Thermo_Structure_Interaction"] = prb_tsi;
    string2prbtype["Thermo_Fluid_Structure_Interaction"] = prb_thermo_fsi;
    string2prbtype["Tutorial"] = prb_tutorial;
    string2prbtype["Two_Phase_Flow"] = prb_two_phase_flow;
    string2prbtype["UQ"] = prb_uq;
    string2prbtype["Variational_Chemical_Diffusion"] = prb_var_chemdiff;
    string2prbtype["XContact"] = prb_xcontact;
  }

  return string2prbtype;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ProblemType INPAR::PROBLEMTYPE::StringToProblemType(std::string name)
{
  std::map<std::string, ProblemType> map = StringToProblemTypeMap();
  std::map<std::string, ProblemType>::const_iterator i = map.find(name);
  if (i != map.end()) return i->second;
  dserror("unsupported problem name '%s'", name.c_str());

  return prb_none;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, ShapeFunctionType> INPAR::PROBLEMTYPE::StringToShapeFunctionTypeMap()
{
  static std::map<std::string, ShapeFunctionType> string2shapefuntype;

  if (string2shapefuntype.size() == 0)
  {
    // problem types in alphabetical order
    string2shapefuntype["Polynomial"] = ShapeFunctionType::shapefunction_polynomial;
    string2shapefuntype["Nurbs"] = ShapeFunctionType::shapefunction_nurbs;
    string2shapefuntype["Meshfree"] = ShapeFunctionType::shapefunction_meshfree;
    string2shapefuntype["HDG"] = ShapeFunctionType::shapefunction_hdg;
  }

  return string2shapefuntype;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ShapeFunctionType INPAR::PROBLEMTYPE::StringToShapeFunctionType(std::string name)
{
  std::map<std::string, ShapeFunctionType> map = StringToShapeFunctionTypeMap();
  std::map<std::string, ShapeFunctionType>::const_iterator i = map.find(name);
  if (i != map.end()) return i->second;
  dserror(
      "'%s' does not name a shape function type. Check for typos or consider adding the shape "
      "function type to the map.",
      name.c_str());

  return ShapeFunctionType::shapefunction_undefined;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::string INPAR::PROBLEMTYPE::ShapeFunctionTypeToString(ShapeFunctionType shapefunctiontype)
{
  std::map<std::string, ShapeFunctionType> map = StringToShapeFunctionTypeMap();

  std::map<ShapeFunctionType, std::string> reversed;
  for (std::map<std::string, ShapeFunctionType>::iterator i = map.begin(); i != map.end(); ++i)
    reversed[i->second] = i->first;

  std::map<ShapeFunctionType, std::string>::const_iterator i = reversed.find(shapefunctiontype);
  if (i != reversed.end()) return i->second;
  dserror(
      "Could not find the name of the given shape function type or the shapefunction is "
      "undefined.");

  return "";
}
