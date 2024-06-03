/*----------------------------------------------------------------------*/
/*! \file

\brief mesh clone strategy for scalar transport problems

\level 1

*/
/*----------------------------------------------------------------------*/


#include "4C_scatra_utils_clonestrategy.hpp"

#include "4C_discretization_fem_general_element.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_scatra_ele.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, std::string> SCATRA::ScatraFluidCloneStrategy::conditions_to_copy() const
{
  return {{"TransportDirichlet", "Dirichlet"}, {"TransportPointNeumann", "PointNeumann"},
      {"TransportLineNeumann", "LineNeumann"}, {"TransportSurfaceNeumann", "SurfaceNeumann"},
      {"TransportVolumeNeumann", "VolumeNeumann"},
      {"TransportNeumannInflow", "TransportNeumannInflow"},
      {"TaylorGalerkinOutflow", "TaylorGalerkinOutflow"},
      {"TaylorGalerkinNeumannInflow", "TaylorGalerkinNeumannInflow"},
      {"ReinitializationTaylorGalerkin", "ReinitializationTaylorGalerkin"},
      {"LinePeriodic", "LinePeriodic"}, {"SurfacePeriodic", "SurfacePeriodic"},
      {"TurbulentInflowSection", "TurbulentInflowSection"}, {"LineNeumann", "FluidLineNeumann"},
      {"SurfaceNeumann", "FluidSurfaceNeumann"}, {"VolumeNeumann", "FluidVolumeNeumann"},
      {"KrylovSpaceProjection", "KrylovSpaceProjection"},
      {"ElchBoundaryKinetics", "ElchBoundaryKinetics"}, {"ScaTraFluxCalc", "ScaTraFluxCalc"},
      {"Initfield", "Initfield"}, {"TransportRobin", "TransportRobin"},
      {"FSICoupling", "FSICoupling"}, {"Mortar", "Mortar"}, {"ScaTraCoupling", "ScaTraCoupling"},
      {"LsContact", "LsContact"}, {"SPRboundary", "SPRboundary"},
      {"XFEMLevelsetTwophase", "XFEMLevelsetTwophase"},
      {"XFEMLevelsetCombustion", "XFEMLevelsetCombustion"}};
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScatraFluidCloneStrategy::check_material_type(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  CORE::Materials::MaterialType mtype =
      GLOBAL::Problem::Instance()->Materials()->ParameterById(matid)->Type();
  if ((mtype != CORE::Materials::m_scatra) && (mtype != CORE::Materials::m_sutherland) &&
      (mtype != CORE::Materials::m_ion) && (mtype != CORE::Materials::m_th_fourier_iso) &&
      (mtype != CORE::Materials::m_thermostvenant) && (mtype != CORE::Materials::m_matlist) &&
      (mtype != CORE::Materials::m_matlist_reactions) && (mtype != CORE::Materials::m_myocard) &&
      (mtype != CORE::Materials::m_scatra_multiporo_fluid) &&
      (mtype != CORE::Materials::m_scatra_multiporo_volfrac) &&
      (mtype != CORE::Materials::m_scatra_multiporo_solid) &&
      (mtype != CORE::Materials::m_scatra_multiporo_temperature))
    FOUR_C_THROW("Material with ID %d is not admissible for scalar transport elements", matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScatraFluidCloneStrategy::set_element_data(
    Teuchos::RCP<CORE::Elements::Element> newele, CORE::Elements::Element* oldele, const int matid,
    const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the transport element!
  DRT::ELEMENTS::Transport* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (trans != nullptr)
  {
    trans->SetMaterial(matid, oldele);
    trans->SetDisType(oldele->Shape());  // set distype as well!
  }
  else
  {
    FOUR_C_THROW("unsupported element type '%s'", typeid(*newele).name());
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SCATRA::ScatraFluidCloneStrategy::determine_ele_type(
    CORE::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // note: ismyele, actele remain unused here! Used only for ALE creation

  // we only support transport elements here
  eletype.push_back("TRANSP");

  return true;  // yes, we copy EVERY element (no submeshes)
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, std::string> SCATRA::ScatraReactionCloneStrategy::conditions_to_copy() const
{
  return {{"TransportDirichlet", "Dirichlet"}, {"TransportPointNeumann", "PointNeumann"},
      {"TransportLineNeumann", "LineNeumann"}, {"TransportSurfaceNeumann", "SurfaceNeumann"},
      {"TransportVolumeNeumann", "VolumeNeumann"},
      {"KrylovSpaceProjection", "KrylovSpaceProjection"}, {"ScaTraFluxCalc", "ScaTraFluxCalc"},
      {"Initfield", "Initfield"}, {"ScaTraCoupling", "ScaTraCoupling"},
      {"SSICoupling", "SSICoupling"}, {"SSICouplingScatraToSolid", "SSICouplingScatraToSolid"},
      {"SSICouplingSolidToScatra", "SSICouplingSolidToScatra"},
      {"ScatraHeteroReactionMaster", "ScatraHeteroReactionMaster"},
      {"ScatraHeteroReactionSlave", "ScatraHeteroReactionSlave"}};
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScatraReactionCloneStrategy::check_material_type(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  CORE::Materials::MaterialType mtype =
      GLOBAL::Problem::Instance()->Materials()->ParameterById(matid)->Type();
  if ((mtype != CORE::Materials::m_scatra) && (mtype != CORE::Materials::m_matlist) &&
      (mtype != CORE::Materials::m_matlist_reactions))
    FOUR_C_THROW("Material with ID %d is not admissible for scalar transport elements", matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScatraReactionCloneStrategy::set_element_data(
    Teuchos::RCP<CORE::Elements::Element> newele, CORE::Elements::Element* oldele, const int matid,
    const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the transport element!
  DRT::ELEMENTS::Transport* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (trans != nullptr)
  {
    trans->SetMaterial(matid, oldele);
    trans->SetDisType(oldele->Shape());  // set distype as well!
  }
  else
  {
    FOUR_C_THROW("unsupported element type '%s'", typeid(*newele).name());
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SCATRA::ScatraReactionCloneStrategy::determine_ele_type(
    CORE::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // note: ismyele, actele remain unused here! Used only for ALE creation

  // we only support transport elements here
  eletype.push_back("TRANSP");

  return true;  // yes, we copy EVERY element (no submeshes)
}

FOUR_C_NAMESPACE_CLOSE
