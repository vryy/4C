// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_art_net_utils.hpp"

#include "4C_art_net_artery.hpp"
#include "4C_art_net_explicitintegration.hpp"
#include "4C_art_net_impl_stationary.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_scatra_ele.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | create algorithm                                                      |
 *----------------------------------------------------------------------*/
std::shared_ptr<Adapter::ArtNet> Arteries::Utils::create_algorithm(
    Inpar::ArtDyn::TimeIntegrationScheme timintscheme,
    std::shared_ptr<Core::FE::Discretization> dis, const int linsolvernumber,
    const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& artparams,
    Core::IO::DiscretizationWriter& output)
{
  // Creation of Coupled Problem algorithm.
  std::shared_ptr<Adapter::ArtNet> algo = nullptr;

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------

  switch (timintscheme)
  {
    case Inpar::ArtDyn::TimeIntegrationScheme::tay_gal:
    {
      // create algorithm
      algo = std::make_shared<Arteries::ArtNetExplicitTimeInt>(
          dis, linsolvernumber, probparams, artparams, output);
      break;
    }
    case Inpar::ArtDyn::TimeIntegrationScheme::stationary:
    {
      // create algorithm
      algo = std::make_shared<Arteries::ArtNetImplStationary>(
          dis, linsolvernumber, probparams, artparams, output);
      break;
    }
    default:
      FOUR_C_THROW("Unknown time-integration scheme for artery network problem");
      break;
  }

  return algo;
}

/*----------------------------------------------------------------------*
 | exchange material pointers of both discretizations  kremheller 03/18 |
 *----------------------------------------------------------------------*/
void Arteries::Utils::assign_material_pointers(
    const std::string& artery_disname, const std::string& scatra_disname)
{
  Global::Problem* problem = Global::Problem::instance();

  std::shared_ptr<Core::FE::Discretization> arterydis = problem->get_dis(artery_disname);
  std::shared_ptr<Core::FE::Discretization> scatradis = problem->get_dis(scatra_disname);

  set_material_pointers_matching_grid(*arterydis, *scatradis);
}

/*----------------------------------------------------------------------*
 | reset Material pointers after redistribution        kremheller 03/18 |
 *----------------------------------------------------------------------*/
void Arteries::Utils::set_material_pointers_matching_grid(
    const Core::FE::Discretization& sourcedis, const Core::FE::Discretization& targetdis)
{
  const int numelements = targetdis.num_my_col_elements();

  for (int i = 0; i < numelements; ++i)
  {
    Core::Elements::Element* targetele = targetdis.l_col_element(i);
    const int gid = targetele->id();

    Core::Elements::Element* sourceele = sourcedis.g_element(gid);

    // for coupling we add the source material to the target element and vice versa
    targetele->add_material(sourceele->material());
    sourceele->add_material(targetele->material());
  }
}

/*----------------------------------------------------------------------*
 |                                                     kremheller 03/18 |
 *----------------------------------------------------------------------*/
bool Arteries::ArteryScatraCloneStrategy::determine_ele_type(
    Core::Elements::Element* actele, const bool ismyele, std::vector<std::string>& eletype)
{
  // clone the element
  Discret::Elements::Artery* myele = static_cast<Discret::Elements::Artery*>(actele);
  // only the pressure based artery supports this function so far
  if (myele->impl_type() == Inpar::ArtDyn::impltype_pressure_based)
  {
    // we only support transport elements here
    eletype.push_back("TRANSP");
    return true;
  }

  return false;
}


/*----------------------------------------------------------------------*
 | set the element data (protected)                    kremheller 03/18 |
 *----------------------------------------------------------------------*/
void Arteries::ArteryScatraCloneStrategy::set_element_data(
    std::shared_ptr<Core::Elements::Element> newele, Core::Elements::Element* oldele,
    const int matid, const bool isnurbs)
{
  // We need to set material and possibly other things to complete element setup.
  // This is again really ugly as we have to extract the actual
  // element type in order to access the material property

  // note: set_material() was reimplemented by the transport element!
  Discret::Elements::Transport* trans = dynamic_cast<Discret::Elements::Transport*>(newele.get());
  if (trans != nullptr)
  {
    // set material
    trans->set_material(matid, oldele);
    // set distype as well!
    trans->set_dis_type(oldele->shape());

    // we only have one possible impltype
    Inpar::ScaTra::ImplType impltype = Inpar::ScaTra::impltype_one_d_artery;
    trans->set_impl_type(impltype);
  }
  else
  {
    FOUR_C_THROW(
        "unsupported element type '%s'", Core::Utils::get_dynamic_type_name(*newele).c_str());
  }
  return;
}

/*----------------------------------------------------------------------*
 | check if material type is admissible                kremheller 03/18 |
 *----------------------------------------------------------------------*/
void Arteries::ArteryScatraCloneStrategy::check_material_type(const int matid)
{
  // We take the material with the ID specified by the user
  // Here we check first, whether this material is of admissible type
  Core::Materials::MaterialType mtype =
      Global::Problem::instance()->materials()->parameter_by_id(matid)->type();
  if ((mtype != Core::Materials::m_scatra) && (mtype != Core::Materials::m_matlist))
    FOUR_C_THROW("Material with ID %d is not admissible for scalar transport elements", matid);
}


/*----------------------------------------------------------------------*
 |                                                     kremheller 03/18 |
 *----------------------------------------------------------------------*/
std::map<std::string, std::string> Arteries::ArteryScatraCloneStrategy::conditions_to_copy() const
{
  return {{"TransportDirichlet", "Dirichlet"}, {"TransportPointNeumann", "PointNeumann"},
      {"Initfield", "Initfield"}, {"ArtScatraCouplConNodebased", "ArtScatraCouplConNodebased"},
      {"PoroMultiphaseScatraOxyPartPressCalcCond", "PoroMultiphaseScatraOxyPartPressCalcCond"},
      {"ArtScatraCouplConNodeToPoint", "ArtScatraCouplConNodeToPoint"}};
}

FOUR_C_NAMESPACE_CLOSE
