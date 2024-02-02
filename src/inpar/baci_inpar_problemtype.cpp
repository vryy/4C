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

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void INPAR::PROBLEMTYPE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace INPUT;
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& type = list->sublist("PROBLEM TYP", false, "");

  {
    using IntegerType = std::underlying_type_t<GLOBAL::ProblemType>;
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
  IntParameter("RANDSEED", -1, "Set the random seed. If < 0 use current time.", &type);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, GLOBAL::ProblemType> INPAR::PROBLEMTYPE::StringToProblemTypeMap()
{
  static std::map<std::string, GLOBAL::ProblemType> string2prbtype;

  if (string2prbtype.size() == 0)
  {
    // problem types in alphabetical order
    string2prbtype["Ale"] = GLOBAL::ProblemType::ale;
    string2prbtype["ArterialNetwork"] = GLOBAL::ProblemType::art_net;
    string2prbtype["Atherosclerosis_Fluid_Structure_Interaction"] = GLOBAL::ProblemType::ac_fsi;
    string2prbtype["Biofilm_Fluid_Structure_Interaction"] = GLOBAL::ProblemType::biofilm_fsi;
    string2prbtype["Cardiac_Monodomain"] = GLOBAL::ProblemType::cardiac_monodomain;
    string2prbtype["Elastohydrodynamic_Lubrication"] = GLOBAL::ProblemType::ehl;
    string2prbtype["Electrochemistry"] = GLOBAL::ProblemType::elch;
    string2prbtype["Electromagnetics"] = GLOBAL::ProblemType::elemag;
    string2prbtype["Fluid"] = GLOBAL::ProblemType::fluid;
    string2prbtype["Fluid_Ale"] = GLOBAL::ProblemType::fluid_ale;
    string2prbtype["Fluid_Beam_Interaction"] = GLOBAL::ProblemType::fbi;
    string2prbtype["Fluid_Freesurface"] = GLOBAL::ProblemType::freesurf;
    string2prbtype["Fluid_Poro_Structure_Interaction_XFEM"] = GLOBAL::ProblemType::fpsi_xfem;
    string2prbtype["Fluid_Porous_Structure_Interaction"] = GLOBAL::ProblemType::fpsi;
    string2prbtype["Fluid_Porous_Structure_Scalar_Scalar_Interaction"] = GLOBAL::ProblemType::fps3i;
    string2prbtype["Fluid_RedModels"] = GLOBAL::ProblemType::fluid_redmodels;
    string2prbtype["Fluid_Structure_Interaction"] = GLOBAL::ProblemType::fsi;
    string2prbtype["Fluid_Structure_Interaction_Lung"] = GLOBAL::ProblemType::fsi_lung;
    string2prbtype["Fluid_Structure_Interaction_RedModels"] = GLOBAL::ProblemType::fsi_redmodels;
    string2prbtype["Fluid_Structure_Interaction_XFEM"] = GLOBAL::ProblemType::fsi_xfem;
    string2prbtype["Fluid_XFEM"] = GLOBAL::ProblemType::fluid_xfem;
    string2prbtype["Fluid_XFEM_LevelSet"] = GLOBAL::ProblemType::fluid_xfem_ls;
    string2prbtype["Gas_Fluid_Structure_Interaction"] = GLOBAL::ProblemType::gas_fsi;
    string2prbtype["Immersed_FSI"] = GLOBAL::ProblemType::immersed_fsi;
    string2prbtype["Level_Set"] = GLOBAL::ProblemType::level_set;
    string2prbtype["Low_Mach_Number_Flow"] = GLOBAL::ProblemType::loma;
    string2prbtype["Lubrication"] = GLOBAL::ProblemType::lubrication;
    string2prbtype["NP_Supporting_Procs"] = GLOBAL::ProblemType::np_support;
    string2prbtype["Particle"] = GLOBAL::ProblemType::particle;
    string2prbtype["Particle_Structure_Interaction"] = GLOBAL::ProblemType::pasi;
    string2prbtype["Polymer_Network"] = GLOBAL::ProblemType::polymernetwork;
    string2prbtype["Poroelastic_scalar_transport"] = GLOBAL::ProblemType::poroscatra;
    string2prbtype["Poroelasticity"] = GLOBAL::ProblemType::poroelast;
    string2prbtype["Multiphase_Poroelasticity"] = GLOBAL::ProblemType::poromultiphase;
    string2prbtype["Multiphase_Poroelasticity_ScaTra"] = GLOBAL::ProblemType::poromultiphasescatra;
    string2prbtype["Multiphase_Porous_Flow"] = GLOBAL::ProblemType::porofluidmultiphase;
    string2prbtype["RedAirways_Tissue"] = GLOBAL::ProblemType::redairways_tissue;
    string2prbtype["ReducedDimensionalAirWays"] = GLOBAL::ProblemType::red_airways;
    string2prbtype["Scalar_Thermo_Interaction"] = GLOBAL::ProblemType::sti;
    string2prbtype["Scalar_Transport"] = GLOBAL::ProblemType::scatra;
    string2prbtype["Structure"] = GLOBAL::ProblemType::structure;
    string2prbtype["Structure_Ale"] = GLOBAL::ProblemType::struct_ale;
    string2prbtype["Structure_Scalar_Interaction"] = GLOBAL::ProblemType::ssi;
    string2prbtype["Structure_Scalar_Thermo_Interaction"] = GLOBAL::ProblemType::ssti;
    string2prbtype["Thermo"] = GLOBAL::ProblemType::thermo;
    string2prbtype["Thermo_Structure_Interaction"] = GLOBAL::ProblemType::tsi;
    string2prbtype["Thermo_Fluid_Structure_Interaction"] = GLOBAL::ProblemType::thermo_fsi;
  }

  return string2prbtype;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
GLOBAL::ProblemType INPAR::PROBLEMTYPE::StringToProblemType(std::string name)
{
  std::map<std::string, GLOBAL::ProblemType> map = StringToProblemTypeMap();
  std::map<std::string, GLOBAL::ProblemType>::const_iterator i = map.find(name);
  if (i != map.end()) return i->second;
  dserror("unsupported problem name '%s'", name.c_str());

  return GLOBAL::ProblemType::none;
}

BACI_NAMESPACE_CLOSE
