// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_beaminteraction_beam_to_beam_contact_condition.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_pair.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_fem_condition.hpp"
#include "4C_fem_discretization.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BeamInteraction::BeamToBeamContactCondition::BeamToBeamContactCondition(
    const std::shared_ptr<const Core::Conditions::Condition>& condition_line_1,
    const std::shared_ptr<const Core::Conditions::Condition>& condition_line_2)
    : BeamInteractionConditionBase(condition_line_1),
      condition_other_(condition_line_2),
      condition_contact_pairs_(),
      other_line_ids_()
{
}

/**
 *
 */
void BeamInteraction::BeamToBeamContactCondition::build_id_sets(
    const std::shared_ptr<const Core::FE::Discretization>& discretization)
{
  // Call the parent method to build the line maps.
  BeamInteractionConditionBase::build_id_sets(discretization);

  // Build the other line map.
  std::vector<int> line_ids;
  condition_to_element_ids(*condition_other_, line_ids);
  other_line_ids_ = std::set<int>(line_ids.begin(), line_ids.end());
}

/**
 *
 */
bool BeamInteraction::BeamToBeamContactCondition::ids_in_condition(
    const int id_line, const int id_other) const
{
  if (id_is_in_condition(line_ids_, id_line) and id_is_in_condition(other_line_ids_, id_other))
    return true;
  if (id_is_in_condition(line_ids_, id_other) and id_is_in_condition(other_line_ids_, id_line))
    return true;
  return false;
}

/**
 *
 */
void BeamInteraction::BeamToBeamContactCondition::clear()
{
  BeamInteractionConditionBase::clear();
  condition_contact_pairs_.clear();
}

/**
 *
 */
std::shared_ptr<BeamInteraction::BeamContactPair>
BeamInteraction::BeamToBeamContactCondition::create_contact_pair(
    const std::vector<Core::Elements::Element const*>& ele_ptrs)
{
  // Check if the given elements are in this condition.
  if (!ids_in_condition(ele_ptrs[0]->id(), ele_ptrs[1]->id())) return nullptr;

  // note: numnodes is to be interpreted as number of nodes used for centerline interpolation.
  // numnodalvalues = 1: only positions as primary nodal DoFs ==> Lagrange interpolation
  // numnodalvalues = 2: positions AND tangents ==> Hermite interpolation

  const Discret::Elements::Beam3Base* beamele1 =
      dynamic_cast<const Discret::Elements::Beam3Base*>(ele_ptrs[0]);

  const unsigned int numnodes_centerline = beamele1->num_centerline_nodes();
  const unsigned int numnodalvalues = beamele1->hermite_centerline_interpolation() ? 2 : 1;

  switch (numnodalvalues)
  {
    case 1:
    {
      switch (numnodes_centerline)
      {
        case 2:
        {
          return std::make_shared<BeamInteraction::BeamToBeamContactPair<2, 1>>();
        }
        case 3:
        {
          return std::make_shared<BeamInteraction::BeamToBeamContactPair<3, 1>>();
        }
        case 4:
        {
          return std::make_shared<BeamInteraction::BeamToBeamContactPair<4, 1>>();
        }
        case 5:
        {
          return std::make_shared<BeamInteraction::BeamToBeamContactPair<5, 1>>();
        }
        default:
        {
          FOUR_C_THROW(
              "%d and %d is no valid template parameter combination for the "
              "number of nodes and number of types of nodal DoFs used for centerline "
              "interpolation!",
              numnodes_centerline, numnodalvalues);
          break;
        }
      }
      break;
    }
    case 2:
    {
      switch (numnodes_centerline)
      {
        case 2:
        {
          return std::make_shared<BeamInteraction::BeamToBeamContactPair<2, 2>>();
        }
        default:
        {
          FOUR_C_THROW(
              "%d and %d is no valid template parameter combination for the "
              "number of nodes and number of types of nodal DoFs used for centerline "
              "interpolation!",
              numnodes_centerline, numnodalvalues);
          break;
        }
      }
      break;
    }
    default:
    {
      FOUR_C_THROW(
          "%d and %d is no valid template parameter combination for the "
          "number of nodes and number of types of nodal DoFs used for centerline "
          "interpolation!",
          numnodes_centerline, numnodalvalues);
      break;
    }
  }

  return nullptr;
}

FOUR_C_NAMESPACE_CLOSE
