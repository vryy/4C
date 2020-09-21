/*----------------------------------------------------------------------------*/
/*! \file

\brief Strategy to clone scatra discretization from elemag discretization

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "elemag_utils_clonestrategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_mat/matpar_material.H"
#include "../drt_mat/matpar_bundle.H"

// we need to know all element types for the scatra mesh creation
#include "../drt_elemag/elemag_ele.H"
#include "../drt_scatra_ele/scatra_ele.H"
#include "../drt_scatra_ele/scatra_ele_hdg.H"

#include "../drt_lib/drt_globalproblem_enums.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <ShapeFunctionType sft>
std::map<std::string, std::string> ELEMAG::UTILS::ScatraCloneStrategy<sft>::ConditionsToCopy()
{
  std::map<std::string, std::string> conditions_to_copy;

  conditions_to_copy.insert(std::pair<std::string, std::string>("Dirichlet", "Dirichlet"));

  return conditions_to_copy;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <ShapeFunctionType sft>
void ELEMAG::UTILS::ScatraCloneStrategy<sft>::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = DRT::Problem::Instance()->Materials()->ById(matid)->Type();
  if (mtype != INPAR::MAT::m_scatra)
    dserror("Material with ID %d is not admissible for TRANSP elements", matid);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <ShapeFunctionType sft>
void ELEMAG::UTILS::ScatraCloneStrategy<sft>::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool nurbsdis)
{
  auto Transport = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (Transport != NULL)
  {
    Transport->SetDisType(oldele->Shape());
    Transport->SetMaterial(matid);
    if (sft == ShapeFunctionType::shapefunction_hdg)
    {
      auto scatraele = dynamic_cast<DRT::ELEMENTS::ScaTraHDG*>(Transport);
      scatraele->SetImplType(INPAR::SCATRA::impltype_std_hdg);
      scatraele->SetDegree(oldele->Degree());
      scatraele->SetCompletePolynomialSpace(false);
    }
    else
    {
      Transport->SetImplType(INPAR::SCATRA::impltype_std);
    }
  }
  else
    dserror("unsupported ale element type '%s'", typeid(*newele).name());

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <ShapeFunctionType sft>
bool ELEMAG::UTILS::ScatraCloneStrategy<sft>::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // Clone it now.
  if (ismyele)
  {
    sft == ShapeFunctionType::shapefunction_hdg ? eletype.push_back("TRANSPHDG")
                                                : eletype.push_back("TRANSP");
  }

  return true;
}

// template classes
template class ELEMAG::UTILS::ScatraCloneStrategy<ShapeFunctionType::shapefunction_polynomial>;
template class ELEMAG::UTILS::ScatraCloneStrategy<ShapeFunctionType::shapefunction_hdg>;