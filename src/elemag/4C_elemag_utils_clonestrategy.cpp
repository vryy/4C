// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

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
std::map<std::string, std::string> EleMag::Utils::ScatraCloneStrategy<sft>::conditions_to_copy()
    const
{
  return {{"Dirichlet", "Dirichlet"}};
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <Core::FE::ShapeFunctionType sft>
void EleMag::Utils::ScatraCloneStrategy<sft>::check_material_type(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  Core::Materials::MaterialType mtype =
      Global::Problem::instance()->materials()->parameter_by_id(matid)->type();
  if (mtype != Core::Materials::m_scatra)
    FOUR_C_THROW("Material with ID %d is not admissible for TRANSP elements", matid);
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <Core::FE::ShapeFunctionType sft>
void EleMag::Utils::ScatraCloneStrategy<sft>::set_element_data(
    std::shared_ptr<Core::Elements::Element> newele, Core::Elements::Element* oldele,
    const int matid, const bool nurbsdis)
{
  auto Transport = dynamic_cast<Discret::Elements::Transport*>(newele.get());
  if (Transport != nullptr)
  {
    Transport->set_dis_type(oldele->shape());
    Transport->set_material(0, Mat::factory(matid));
    if (sft == Core::FE::ShapeFunctionType::hdg)
    {
      auto scatraele = dynamic_cast<Discret::Elements::ScaTraHDG*>(Transport);
      scatraele->set_impl_type(Inpar::ScaTra::impltype_std_hdg);
      scatraele->set_degree(oldele->degree());
      scatraele->set_complete_polynomial_space(false);
    }
    else
    {
      Transport->set_impl_type(Inpar::ScaTra::impltype_std);
    }
  }
  else
    FOUR_C_THROW(
        "unsupported ale element type '%s'", Core::Utils::get_dynamic_type_name(*newele).c_str());

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
template <Core::FE::ShapeFunctionType sft>
bool EleMag::Utils::ScatraCloneStrategy<sft>::determine_ele_type(
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
template class EleMag::Utils::ScatraCloneStrategy<Core::FE::ShapeFunctionType::polynomial>;
template class EleMag::Utils::ScatraCloneStrategy<Core::FE::ShapeFunctionType::hdg>;
FOUR_C_NAMESPACE_CLOSE
