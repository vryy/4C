/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_utils_clonestrategy.cpp

 \brief utils methods for cloning the porofluid discretization

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/


#include "poromultiphase_utils_clonestrategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_porofluidmultiphase_ele/porofluidmultiphase_ele.H"
#include "../drt_lib/drt_element.H"


/*----------------------------------------------------------------------*
 | define conditions to copy to the cloned discretization    vuong 08/16 |
 *----------------------------------------------------------------------*/
std::map<std::string, std::string>
POROMULTIPHASE::UTILS::PoroFluidMultiPhaseCloneStrategy::ConditionsToCopy()
{
  std::map<std::string, std::string> conditions_to_copy;

  conditions_to_copy.insert(std::pair<std::string, std::string>("PoroDirichlet", "Dirichlet"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("PoroPointNeumann", "PointNeumann"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("PoroLineNeumann", "LineNeumann"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("PoroSurfaceNeumann", "SurfaceNeumann"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("PoroVolumeNeumann", "VolumeNeumann"));
  conditions_to_copy.insert(std::pair<std::string, std::string>("Initfield", "Initfield"));
  conditions_to_copy.insert(
      std::pair<std::string, std::string>("ArtPorofluidCouplCon", "ArtPorofluidCouplCon"));

  return conditions_to_copy;
}


/*----------------------------------------------------------------------*
 | check for correct material                               vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::UTILS::PoroFluidMultiPhaseCloneStrategy::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if ((mtype != INPAR::MAT::m_fluidporo_multiphase) and
      (mtype != INPAR::MAT::m_fluidporo_multiphase_reactions))
    dserror("Material with ID %d is not admissible for porofluid multiphase elements", matid);
}


/*----------------------------------------------------------------------*
 | set element-specific data (material etc.)                 vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::UTILS::PoroFluidMultiPhaseCloneStrategy::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: SetMaterial() was reimplemented by the PoroFluidMultiPhase element!
  DRT::ELEMENTS::PoroFluidMultiPhase* porofluidele =
      dynamic_cast<DRT::ELEMENTS::PoroFluidMultiPhase*>(newele.get());
  if (porofluidele != NULL)
  {
    porofluidele->SetMaterial(matid);
    porofluidele->SetDisType(oldele->Shape());  // set distype as well!
  }
  else
  {
    dserror("unsupported element type '%s'", typeid(*newele).name());
  }
  return;
}


/*----------------------------------------------------------------------*
 | determine whether element is copied or not               vuong 08/16 |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASE::UTILS::PoroFluidMultiPhaseCloneStrategy::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // note: ismyele, actele remain unused here! Used only for ALE creation

  // we only support POROFLUIDMULTIPHASE elements here
  eletype.push_back("POROFLUIDMULTIPHASE");

  // all elements are copied
  return true;
}
