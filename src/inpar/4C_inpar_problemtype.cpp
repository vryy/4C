/*----------------------------------------------------------------------*/
/*! \file
\brief convert problem type string to enum
\level 1
*/

/*----------------------------------------------------------------------*/

#include "4C_inpar_problemtype.hpp"

#include "4C_discretization_fem_general_shape_function_type.hpp"
#include "4C_utils_exceptions.hpp"
#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void INPAR::PROBLEMTYPE::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using Teuchos::setStringToIntegralParameter;
  using Teuchos::tuple;

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& type = list->sublist("PROBLEM TYP", false, "");

  {
    using IntegerType = std::underlying_type_t<CORE::ProblemType>;
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

  CORE::UTILS::IntParameter("RESTART", 0, "", &type);
  CORE::UTILS::IntParameter("RANDSEED", -1, "Set the random seed. If < 0 use current time.", &type);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, CORE::ProblemType> INPAR::PROBLEMTYPE::StringToProblemTypeMap()
{
  static std::map<std::string, CORE::ProblemType> string2prbtype;

  if (string2prbtype.size() == 0)
  {
    // problem types in alphabetical order
    string2prbtype["Ale"] = CORE::ProblemType::ale;
    string2prbtype["ArterialNetwork"] = CORE::ProblemType::art_net;
    string2prbtype["Atherosclerosis_Fluid_Structure_Interaction"] = CORE::ProblemType::ac_fsi;
    string2prbtype["Biofilm_Fluid_Structure_Interaction"] = CORE::ProblemType::biofilm_fsi;
    string2prbtype["Cardiac_Monodomain"] = CORE::ProblemType::cardiac_monodomain;
    string2prbtype["Elastohydrodynamic_Lubrication"] = CORE::ProblemType::ehl;
    string2prbtype["Electrochemistry"] = CORE::ProblemType::elch;
    string2prbtype["Electromagnetics"] = CORE::ProblemType::elemag;
    string2prbtype["Fluid"] = CORE::ProblemType::fluid;
    string2prbtype["Fluid_Ale"] = CORE::ProblemType::fluid_ale;
    string2prbtype["Fluid_Beam_Interaction"] = CORE::ProblemType::fbi;
    string2prbtype["Fluid_Freesurface"] = CORE::ProblemType::freesurf;
    string2prbtype["Fluid_Poro_Structure_Interaction_XFEM"] = CORE::ProblemType::fpsi_xfem;
    string2prbtype["Fluid_Porous_Structure_Interaction"] = CORE::ProblemType::fpsi;
    string2prbtype["Fluid_Porous_Structure_Scalar_Scalar_Interaction"] = CORE::ProblemType::fps3i;
    string2prbtype["Fluid_RedModels"] = CORE::ProblemType::fluid_redmodels;
    string2prbtype["Fluid_Structure_Interaction"] = CORE::ProblemType::fsi;
    string2prbtype["Fluid_Structure_Interaction_Lung"] = CORE::ProblemType::fsi_lung;
    string2prbtype["Fluid_Structure_Interaction_RedModels"] = CORE::ProblemType::fsi_redmodels;
    string2prbtype["Fluid_Structure_Interaction_XFEM"] = CORE::ProblemType::fsi_xfem;
    string2prbtype["Fluid_XFEM"] = CORE::ProblemType::fluid_xfem;
    string2prbtype["Fluid_XFEM_LevelSet"] = CORE::ProblemType::fluid_xfem_ls;
    string2prbtype["Gas_Fluid_Structure_Interaction"] = CORE::ProblemType::gas_fsi;
    string2prbtype["Immersed_FSI"] = CORE::ProblemType::immersed_fsi;
    string2prbtype["Level_Set"] = CORE::ProblemType::level_set;
    string2prbtype["Low_Mach_Number_Flow"] = CORE::ProblemType::loma;
    string2prbtype["Lubrication"] = CORE::ProblemType::lubrication;
    string2prbtype["NP_Supporting_Procs"] = CORE::ProblemType::np_support;
    string2prbtype["Particle"] = CORE::ProblemType::particle;
    string2prbtype["Particle_Structure_Interaction"] = CORE::ProblemType::pasi;
    string2prbtype["Polymer_Network"] = CORE::ProblemType::polymernetwork;
    string2prbtype["Poroelastic_scalar_transport"] = CORE::ProblemType::poroscatra;
    string2prbtype["Poroelasticity"] = CORE::ProblemType::poroelast;
    string2prbtype["Multiphase_Poroelasticity"] = CORE::ProblemType::poromultiphase;
    string2prbtype["Multiphase_Poroelasticity_ScaTra"] = CORE::ProblemType::poromultiphasescatra;
    string2prbtype["Multiphase_Porous_Flow"] = CORE::ProblemType::porofluidmultiphase;
    string2prbtype["RedAirways_Tissue"] = CORE::ProblemType::redairways_tissue;
    string2prbtype["ReducedDimensionalAirWays"] = CORE::ProblemType::red_airways;
    string2prbtype["Scalar_Thermo_Interaction"] = CORE::ProblemType::sti;
    string2prbtype["Scalar_Transport"] = CORE::ProblemType::scatra;
    string2prbtype["Structure"] = CORE::ProblemType::structure;
    string2prbtype["Structure_Ale"] = CORE::ProblemType::struct_ale;
    string2prbtype["Structure_Scalar_Interaction"] = CORE::ProblemType::ssi;
    string2prbtype["Structure_Scalar_Thermo_Interaction"] = CORE::ProblemType::ssti;
    string2prbtype["Thermo"] = CORE::ProblemType::thermo;
    string2prbtype["Thermo_Structure_Interaction"] = CORE::ProblemType::tsi;
    string2prbtype["Thermo_Fluid_Structure_Interaction"] = CORE::ProblemType::thermo_fsi;
  }

  return string2prbtype;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CORE::ProblemType INPAR::PROBLEMTYPE::StringToProblemType(std::string name)
{
  std::map<std::string, CORE::ProblemType> map = StringToProblemTypeMap();
  std::map<std::string, CORE::ProblemType>::const_iterator i = map.find(name);
  if (i != map.end()) return i->second;
  FOUR_C_THROW("unsupported problem name '%s'", name.c_str());

  return CORE::ProblemType::none;
}

FOUR_C_NAMESPACE_CLOSE
