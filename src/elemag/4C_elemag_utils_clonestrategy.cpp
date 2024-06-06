/*----------------------------------------------------------------------------*/
/*! \file

\brief Strategy to clone scatra discretization from elemag discretization

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_elemag_utils_clonestrategy.hpp"

#include "4C_elemag_ele.hpp"
#include "4C_global_data.hpp"
#include "4C_legacy_enum_definitions_problem_type.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_hdg.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <Core::FE::ShapeFunctionType sft>
std::map<std::string, std::string> EleMag::UTILS::ScatraCloneStrategy<sft>::conditions_to_copy()
    const
{
  return {{"Dirichlet", "Dirichlet"}};
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <Core::FE::ShapeFunctionType sft>
void EleMag::UTILS::ScatraCloneStrategy<sft>::check_material_type(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  Core::Materials::MaterialType mtype =
      Global::Problem::Instance()->Materials()->ParameterById(matid)->Type();
  if (mtype != Core::Materials::m_scatra)
    FOUR_C_THROW("Material with ID %d is not admissible for TRANSP elements", matid);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <Core::FE::ShapeFunctionType sft>
void EleMag::UTILS::ScatraCloneStrategy<sft>::set_element_data(
    Teuchos::RCP<Core::Elements::Element> newele, Core::Elements::Element* oldele, const int matid,
    const bool nurbsdis)
{
  auto Transport = dynamic_cast<Discret::ELEMENTS::Transport*>(newele.get());
  if (Transport != nullptr)
  {
    Transport->SetDisType(oldele->Shape());
    Transport->SetMaterial(0, Mat::Factory(matid));
    if (sft == Core::FE::ShapeFunctionType::hdg)
    {
      auto scatraele = dynamic_cast<Discret::ELEMENTS::ScaTraHDG*>(Transport);
      scatraele->SetImplType(Inpar::ScaTra::impltype_std_hdg);
      scatraele->SetDegree(oldele->Degree());
      scatraele->set_complete_polynomial_space(false);
    }
    else
    {
      Transport->SetImplType(Inpar::ScaTra::impltype_std);
    }
  }
  else
    FOUR_C_THROW("unsupported ale element type '%s'", typeid(*newele).name());

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <Core::FE::ShapeFunctionType sft>
bool EleMag::UTILS::ScatraCloneStrategy<sft>::determine_ele_type(
    Core::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // Clone it now.
  if (ismyele)
  {
    sft == Core::FE::ShapeFunctionType::hdg ? eletype.push_back("TRANSPHDG")
                                            : eletype.push_back("TRANSP");
  }

  return true;
}

// template classes
template class EleMag::UTILS::ScatraCloneStrategy<Core::FE::ShapeFunctionType::polynomial>;
template class EleMag::UTILS::ScatraCloneStrategy<Core::FE::ShapeFunctionType::hdg>;
FOUR_C_NAMESPACE_CLOSE
