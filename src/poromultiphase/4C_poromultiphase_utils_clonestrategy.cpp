// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_poromultiphase_utils_clonestrategy.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_porofluidmultiphase_ele.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | define conditions to copy to the cloned discretization    vuong 08/16 |
 *----------------------------------------------------------------------*/
std::map<std::string, std::string>
POROMULTIPHASE::Utils::PoroFluidMultiPhaseCloneStrategy::conditions_to_copy() const
{
  return {{"PoroDirichlet", "Dirichlet"}, {"PoroPointNeumann", "PointNeumann"},
      {"PoroLineNeumann", "LineNeumann"}, {"PoroSurfaceNeumann", "SurfaceNeumann"},
      {"PoroVolumeNeumann", "VolumeNeumann"}, {"Initfield", "Initfield"},
      {"ArtPorofluidCouplConNodebased", "ArtPorofluidCouplConNodebased"},
      {"ArtPorofluidCouplConNodeToPoint", "ArtPorofluidCouplConNodeToPoint"}};
}


/*----------------------------------------------------------------------*
 | check for correct material                               vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::Utils::PoroFluidMultiPhaseCloneStrategy::check_material_type(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  Core::Materials::MaterialType mtype =
      Global::Problem::instance()->materials()->parameter_by_id(matid)->type();
  if ((mtype != Core::Materials::m_fluidporo_multiphase) and
      (mtype != Core::Materials::m_fluidporo_multiphase_reactions))
    FOUR_C_THROW("Material with ID %d is not admissible for porofluid multiphase elements", matid);
}


/*----------------------------------------------------------------------*
 | set element-specific data (material etc.)                 vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::Utils::PoroFluidMultiPhaseCloneStrategy::set_element_data(
    std::shared_ptr<Core::Elements::Element> newele, Core::Elements::Element* oldele,
    const int matid, const bool isnurbsdis)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: set_material() was reimplemented by the PoroFluidMultiPhase element!
  Discret::Elements::PoroFluidMultiPhase* porofluidele =
      dynamic_cast<Discret::Elements::PoroFluidMultiPhase*>(newele.get());
  if (porofluidele != nullptr)
  {
    porofluidele->set_material(0, Mat::factory(matid));
    porofluidele->set_dis_type(oldele->shape());  // set distype as well!
  }
  else
  {
    FOUR_C_THROW(
        "unsupported element type '%s'", Core::Utils::get_dynamic_type_name(*newele).c_str());
  }
  return;
}


/*----------------------------------------------------------------------*
 | determine whether element is copied or not               vuong 08/16 |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASE::Utils::PoroFluidMultiPhaseCloneStrategy::determine_ele_type(
    Core::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // note: ismyele, actele remain unused here! Used only for ALE creation

  // we only support POROFLUIDMULTIPHASE elements here
  eletype.push_back("POROFLUIDMULTIPHASE");

  // all elements are copied
  return true;
}

FOUR_C_NAMESPACE_CLOSE
