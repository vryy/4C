/*----------------------------------------------------------------------*/
/*! \file
\brief convert problem type string to enum
\level 1
*/

/*----------------------------------------------------------------------*/

#include "baci_inpar_problemtype.H"

#include "baci_discretization_fem_general_shape_function_type.H"
#include "baci_inpar_validparameters.H"
#include "baci_utils_exceptions.H"

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
    using IntegerType = std::underlying_type_t<ProblemType>;
    Teuchos::Array<std::string> name;
    Teuchos::Array<IntegerType> label;

    for (const auto& [prb_name, prb_enum] : StringToProblemTypeMap())
    {
      name.push_back(prb_name);
      label.push_back(static_cast<IntegerType>(prb_enum));
    }

    setStringToIntegralParameter<IntegerType>(
        "PROBLEMTYP", "Fluid_Structure_Interaction", "", name, label, &type);
  }

  {
    using IntegerType = std::underlying_type_t<CORE::FE::ShapeFunctionType>;
    Teuchos::Array<std::string> name;
    Teuchos::Array<IntegerType> label;

    for (const auto& [prb_name, prb_enum] : CORE::FE::StringToShapeFunctionTypeMap())
    {
      name.push_back(prb_name);
      label.push_back(static_cast<IntegerType>(prb_enum));
    }

    setStringToIntegralParameter<IntegerType>("SHAPEFCT", "Polynomial",
        "Defines the function spaces for the spatial approximation", name, label, &type);
  }

  IntParameter("RESTART", 0, "", &type);
  DoubleParameter("RESTARTTIME", -1.0, "Used defined restart time", &type);
  IntParameter("RANDSEED", -1, "Set the random seed. If < 0 use current time.", &type);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, ProblemType> INPAR::PROBLEMTYPE::StringToProblemTypeMap()
{
  static std::map<std::string, ProblemType> string2prbtype;

  if (string2prbtype.size() == 0)
  {
    // problem types in alphabetical order
    string2prbtype["Ale"] = ProblemType::ale;
    string2prbtype["ArterialNetwork"] = ProblemType::art_net;
    string2prbtype["Atherosclerosis_Fluid_Structure_Interaction"] = ProblemType::ac_fsi;
    string2prbtype["Biofilm_Fluid_Structure_Interaction"] = ProblemType::biofilm_fsi;
    string2prbtype["Cardiac_Monodomain"] = ProblemType::cardiac_monodomain;
    string2prbtype["Elastohydrodynamic_Lubrication"] = ProblemType::ehl;
    string2prbtype["Electrochemistry"] = ProblemType::elch;
    string2prbtype["Electromagnetics"] = ProblemType::elemag;
    string2prbtype["Fluid"] = ProblemType::fluid;
    string2prbtype["Fluid_Ale"] = ProblemType::fluid_ale;
    string2prbtype["Fluid_Beam_Interaction"] = ProblemType::fbi;
    string2prbtype["Fluid_Freesurface"] = ProblemType::freesurf;
    string2prbtype["Fluid_Poro_Structure_Interaction_XFEM"] = ProblemType::fpsi_xfem;
    string2prbtype["Fluid_Porous_Structure_Interaction"] = ProblemType::fpsi;
    string2prbtype["Fluid_Porous_Structure_Scalar_Scalar_Interaction"] = ProblemType::fps3i;
    string2prbtype["Fluid_RedModels"] = ProblemType::fluid_redmodels;
    string2prbtype["Fluid_Structure_Interaction"] = ProblemType::fsi;
    string2prbtype["Fluid_Structure_Interaction_Lung"] = ProblemType::fsi_lung;
    string2prbtype["Fluid_Structure_Interaction_RedModels"] = ProblemType::fsi_redmodels;
    string2prbtype["Fluid_Structure_Interaction_XFEM"] = ProblemType::fsi_xfem;
    string2prbtype["Fluid_XFEM"] = ProblemType::fluid_xfem;
    string2prbtype["Fluid_XFEM_LevelSet"] = ProblemType::fluid_xfem_ls;
    string2prbtype["Gas_Fluid_Structure_Interaction"] = ProblemType::gas_fsi;
    string2prbtype["Immersed_FSI"] = ProblemType::immersed_fsi;
    string2prbtype["Level_Set"] = ProblemType::level_set;
    string2prbtype["Low_Mach_Number_Flow"] = ProblemType::loma;
    string2prbtype["Lubrication"] = ProblemType::lubrication;
    string2prbtype["NP_Supporting_Procs"] = ProblemType::np_support;
    string2prbtype["Particle"] = ProblemType::particle;
    string2prbtype["Particle_Structure_Interaction"] = ProblemType::pasi;
    string2prbtype["Polymer_Network"] = ProblemType::polymernetwork;
    string2prbtype["Poroelastic_scalar_transport"] = ProblemType::poroscatra;
    string2prbtype["Poroelasticity"] = ProblemType::poroelast;
    string2prbtype["Multiphase_Poroelasticity"] = ProblemType::poromultiphase;
    string2prbtype["Multiphase_Poroelasticity_ScaTra"] = ProblemType::poromultiphasescatra;
    string2prbtype["Multiphase_Porous_Flow"] = ProblemType::porofluidmultiphase;
    string2prbtype["RedAirways_Tissue"] = ProblemType::redairways_tissue;
    string2prbtype["ReducedDimensionalAirWays"] = ProblemType::red_airways;
    string2prbtype["Scalar_Thermo_Interaction"] = ProblemType::sti;
    string2prbtype["Scalar_Transport"] = ProblemType::scatra;
    string2prbtype["Structure"] = ProblemType::structure;
    string2prbtype["Structure_Ale"] = ProblemType::struct_ale;
    string2prbtype["Structure_Scalar_Interaction"] = ProblemType::ssi;
    string2prbtype["Structure_Scalar_Thermo_Interaction"] = ProblemType::ssti;
    string2prbtype["Thermo"] = ProblemType::thermo;
    string2prbtype["Thermo_Structure_Interaction"] = ProblemType::tsi;
    string2prbtype["Thermo_Fluid_Structure_Interaction"] = ProblemType::thermo_fsi;
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

  return ProblemType::none;
}
