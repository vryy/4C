/*----------------------------------------------------------------------*/
/*! \file

\brief strategy for cloning thermo discretization from scatra discretization

\level 2

*/
/*----------------------------------------------------------------------*/
#include "4C_sti_clonestrategy.hpp"

#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_scatra_ele.hpp"

FOUR_C_NAMESPACE_OPEN

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
void STI::ScatraThermoCloneStrategy::CheckMaterialType(const int matid)
{
  // check whether material with specified ID is compatible with cloned element or not
  switch (GLOBAL::Problem::Instance()->Materials()->ById(matid)->Type())
  {
    case CORE::Materials::m_soret:
    case CORE::Materials::m_th_fourier_iso:
      // do nothing in case of compatible material
      break;

    default:
    {
      // throw error in case of incompatible material
      FOUR_C_THROW("Material with ID %d is not compatible with cloned transport element!", matid);
      break;
    }
  }
}

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
std::map<std::string, std::string> STI::ScatraThermoCloneStrategy::ConditionsToCopy() const
{
  return {{"PointThermoCoupling", "PointCoupling"}, {"S2IKinetics", "S2IKinetics"},
      {"S2IMeshtying", "S2IMeshtying"}, {"ScaTraFluxCalc", "ScaTraFluxCalc"},
      {"ThermoDirichlet", "Dirichlet"}, {"ThermoPointNeumann", "PointNeumann"},
      {"ThermoLineNeumann", "LineNeumann"}, {"ThermoSurfaceNeumann", "SurfaceNeumann"},
      {"ThermoVolumeNeumann", "VolumeNeumann"}, {"ThermoInitfield", "Initfield"},
      {"ThermoRobin", "TransportRobin"}, {"ScatraPartitioning", "ScatraPartitioning"}};
}  // STI::ScatraThermoCloneStrategy::ConditionsToCopy()

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
bool STI::ScatraThermoCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // set type of cloned element to transport type
  eletype.emplace_back("TRANSP");

  // element should always be cloned
  return true;
}  // STI::ScatraThermoCloneStrategy::DetermineEleType

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
void STI::ScatraThermoCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbs)
{
  // cast pointers to current element on source discretization and to current cloned element on
  // target discretization
  auto* oldele_transport = dynamic_cast<DRT::ELEMENTS::Transport*>(oldele);
  Teuchos::RCP<DRT::ELEMENTS::Transport> newele_transport =
      Teuchos::rcp_dynamic_cast<DRT::ELEMENTS::Transport>(newele);

  // safety check
  if (oldele_transport == nullptr or newele_transport == Teuchos::null)
    FOUR_C_THROW(
        "Expected transport element, but received element of type '%s'!", typeid(*newele).name());

  // provide cloned element with material
  newele_transport->SetMaterial(matid, oldele);

  // provide cloned element with discretization type
  newele_transport->SetDisType(oldele->Shape());

  // provide cloned element with physical implementation type
  switch (oldele_transport->ImplType())
  {
    case INPAR::SCATRA::impltype_elch_diffcond_thermo:
    case INPAR::SCATRA::impltype_elch_diffcond:
    {
      newele_transport->SetImplType(INPAR::SCATRA::impltype_thermo_elch_diffcond);
      break;
    }
    case INPAR::SCATRA::impltype_elch_electrode_thermo:
    case INPAR::SCATRA::impltype_elch_electrode:
    {
      newele_transport->SetImplType(INPAR::SCATRA::impltype_thermo_elch_electrode);
      break;
    }
    default:
    {
      FOUR_C_THROW(
          "Scatra-thermo interaction not yet implemented for given element implementation type!");
      break;
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
