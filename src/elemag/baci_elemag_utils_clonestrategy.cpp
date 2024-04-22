/*----------------------------------------------------------------------------*/
/*! \file

\brief Strategy to clone scatra discretization from elemag discretization

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "baci_elemag_utils_clonestrategy.hpp"

#include "baci_elemag_ele.hpp"
#include "baci_global_data.hpp"
#include "baci_global_data_enums.hpp"
#include "baci_mat_par_bundle.hpp"
#include "baci_mat_par_material.hpp"
#include "baci_scatra_ele.hpp"
#include "baci_scatra_ele_hdg.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <CORE::FE::ShapeFunctionType sft>
std::map<std::string, std::string> ELEMAG::UTILS::ScatraCloneStrategy<sft>::ConditionsToCopy() const
{
  return {{"Dirichlet", "Dirichlet"}};
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <CORE::FE::ShapeFunctionType sft>
void ELEMAG::UTILS::ScatraCloneStrategy<sft>::CheckMaterialType(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  INPAR::MAT::MaterialType mtype = GLOBAL::Problem::Instance()->Materials()->ById(matid)->Type();
  if (mtype != INPAR::MAT::m_scatra)
    FOUR_C_THROW("Material with ID %d is not admissible for TRANSP elements", matid);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <CORE::FE::ShapeFunctionType sft>
void ELEMAG::UTILS::ScatraCloneStrategy<sft>::SetElementData(
    Teuchos::RCP<DRT::Element> newele, DRT::Element* oldele, const int matid, const bool nurbsdis)
{
  auto Transport = dynamic_cast<DRT::ELEMENTS::Transport*>(newele.get());
  if (Transport != nullptr)
  {
    Transport->SetDisType(oldele->Shape());
    Transport->SetMaterial(matid);
    if (sft == CORE::FE::ShapeFunctionType::hdg)
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
    FOUR_C_THROW("unsupported ale element type '%s'", typeid(*newele).name());

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <CORE::FE::ShapeFunctionType sft>
bool ELEMAG::UTILS::ScatraCloneStrategy<sft>::DetermineEleType(
    DRT::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // Clone it now.
  if (ismyele)
  {
    sft == CORE::FE::ShapeFunctionType::hdg ? eletype.push_back("TRANSPHDG")
                                            : eletype.push_back("TRANSP");
  }

  return true;
}

// template classes
template class ELEMAG::UTILS::ScatraCloneStrategy<CORE::FE::ShapeFunctionType::polynomial>;
template class ELEMAG::UTILS::ScatraCloneStrategy<CORE::FE::ShapeFunctionType::hdg>;
FOUR_C_NAMESPACE_CLOSE
