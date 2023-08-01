/*----------------------------------------------------------------------*/
/*! \file

\brief mesh clone strategy for scalar transport problems

\level 1

*/
/*----------------------------------------------------------------------*/


#include "baci_scatra_utils_clonestrategy.H"

#include "baci_lib_element.H"
#include "baci_lib_globalproblem.H"
#include "baci_mat_par_bundle.H"
#include "baci_mat_par_material.H"
#include "baci_scatra_ele.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, std::string> SCATRA::ScatraFluidCloneStrategy::ConditionsToCopy() const
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
void SCATRA::ScatraFluidCloneStrategy::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if ((mtype != INPAR::MAT::m_scatra) && (mtype != INPAR::MAT::m_mixfrac) &&
      (mtype != INPAR::MAT::m_sutherland) && (mtype != INPAR::MAT::m_tempdepwater) &&
      (mtype != INPAR::MAT::m_arrhenius_pv) && (mtype != INPAR::MAT::m_ferech_pv) &&
      (mtype != INPAR::MAT::m_ion) && (mtype != INPAR::MAT::m_th_fourier_iso) &&
      (mtype != INPAR::MAT::m_thermostvenant) && (mtype != INPAR::MAT::m_yoghurt) &&
      (mtype != INPAR::MAT::m_matlist) && (mtype != INPAR::MAT::m_matlist_reactions) &&
      (mtype != INPAR::MAT::m_myocard) && (mtype != INPAR::MAT::m_scatra_multiporo_fluid) &&
      (mtype != INPAR::MAT::m_scatra_multiporo_volfrac) &&
      (mtype != INPAR::MAT::m_scatra_multiporo_solid) &&
      (mtype != INPAR::MAT::m_scatra_multiporo_temperature))
    dserror("Material with ID %d is not admissible for scalar transport elements", matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScatraFluidCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the transport element!
  DRT::ELEMENTS::Transport* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (trans != NULL)
  {
    trans->SetMaterial(matid, oldele);
    trans->SetDisType(oldele->Shape());  // set distype as well!
  }
  else
  {
    dserror("unsupported element type '%s'", typeid(*newele).name());
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SCATRA::ScatraFluidCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // note: ismyele, actele remain unused here! Used only for ALE creation

  // we only support transport elements here
  eletype.push_back("TRANSP");

  return true;  // yes, we copy EVERY element (no submeshes)
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::map<std::string, std::string> SCATRA::ScatraReactionCloneStrategy::ConditionsToCopy() const
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
void SCATRA::ScatraReactionCloneStrategy::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if ((mtype != INPAR::MAT::m_scatra) && (mtype != INPAR::MAT::m_matlist) &&
      (mtype != INPAR::MAT::m_matlist_reactions))
    dserror("Material with ID %d is not admissible for scalar transport elements", matid);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SCATRA::ScatraReactionCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the transport element!
  DRT::ELEMENTS::Transport* trans = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (trans != NULL)
  {
    trans->SetMaterial(matid, oldele);
    trans->SetDisType(oldele->Shape());  // set distype as well!
  }
  else
  {
    dserror("unsupported element type '%s'", typeid(*newele).name());
  }
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool SCATRA::ScatraReactionCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // note: ismyele, actele remain unused here! Used only for ALE creation

  // we only support transport elements here
  eletype.push_back("TRANSP");

  return true;  // yes, we copy EVERY element (no submeshes)
}
