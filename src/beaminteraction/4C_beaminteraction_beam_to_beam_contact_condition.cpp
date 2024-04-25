/*----------------------------------------------------------------------*/
/*! \file

\brief Object to hold and manage a beam-to-beam condition.

\level 3
*/


#include "4C_beaminteraction_beam_to_beam_contact_condition.hpp"

#include "4C_beam3_base.hpp"
#include "4C_beaminteraction_beam_to_beam_contact_pair.hpp"
#include "4C_beaminteraction_contact_pair.hpp"
#include "4C_beaminteraction_contact_params.hpp"
#include "4C_lib_condition.hpp"
#include "4C_lib_discret.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
BEAMINTERACTION::BeamToBeamContactCondition::BeamToBeamContactCondition(
    const Teuchos::RCP<const DRT::Condition>& condition_line_1,
    const Teuchos::RCP<const DRT::Condition>& condition_line_2)
    : BeamInteractionConditionBase(condition_line_1),
      condition_other_(condition_line_2),
      condition_contact_pairs_(),
      other_line_ids_()
{
}

/**
 *
 */
void BEAMINTERACTION::BeamToBeamContactCondition::BuildIdSets(
    const Teuchos::RCP<const DRT::Discretization>& discretization)
{
  // Call the parent method to build the line maps.
  BeamInteractionConditionBase::BuildIdSets(discretization);

  // Build the other line map.
  std::vector<int> line_ids;
  ConditionToElementIds(condition_other_, line_ids);
  other_line_ids_ = std::set<int>(line_ids.begin(), line_ids.end());
}

/**
 *
 */
bool BEAMINTERACTION::BeamToBeamContactCondition::IdsInCondition(
    const int id_line, const int id_other) const
{
  if (IdIsInCondition(line_ids_, id_line) and IdIsInCondition(other_line_ids_, id_other))
    return true;
  if (IdIsInCondition(line_ids_, id_other) and IdIsInCondition(other_line_ids_, id_line))
    return true;
  return false;
}

/**
 *
 */
void BEAMINTERACTION::BeamToBeamContactCondition::Clear()
{
  BeamInteractionConditionBase::Clear();
  condition_contact_pairs_.clear();
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToBeamContactCondition::CreateContactPair(
    const std::vector<DRT::Element const*>& ele_ptrs)
{
  // Check if the given elements are in this condition.
  if (!IdsInCondition(ele_ptrs[0]->Id(), ele_ptrs[1]->Id())) return Teuchos::null;

  // note: numnodes is to be interpreted as number of nodes used for centerline interpolation.
  // numnodalvalues = 1: only positions as primary nodal DoFs ==> Lagrange interpolation
  // numnodalvalues = 2: positions AND tangents ==> Hermite interpolation

  const DRT::ELEMENTS::Beam3Base* beamele1 =
      dynamic_cast<const DRT::ELEMENTS::Beam3Base*>(ele_ptrs[0]);

  const unsigned int numnodes_centerline = beamele1->NumCenterlineNodes();
  const unsigned int numnodalvalues = beamele1->HermiteCenterlineInterpolation() ? 2 : 1;

  switch (numnodalvalues)
  {
    case 1:
    {
      switch (numnodes_centerline)
      {
        case 2:
        {
          return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamContactPair<2, 1>());
        }
        case 3:
        {
          return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamContactPair<3, 1>());
        }
        case 4:
        {
          return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamContactPair<4, 1>());
        }
        case 5:
        {
          return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamContactPair<5, 1>());
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
          return Teuchos::rcp(new BEAMINTERACTION::BeamToBeamContactPair<2, 2>());
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

  return Teuchos::null;
}

FOUR_C_NAMESPACE_CLOSE
