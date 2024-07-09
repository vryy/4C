/*----------------------------------------------------------------------*/
/*! \file

\brief Class to manage beam-to-beam point couplings.

\level 3
*/


#include "4C_beaminteraction_beam_to_beam_point_coupling_pair_condition.hpp"

#include "4C_beaminteraction_beam_to_beam_point_coupling_pair.hpp"
#include "4C_fem_condition.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_geometry_pair_element.hpp"

FOUR_C_NAMESPACE_OPEN


/**
 *
 */
bool BEAMINTERACTION::BeamToBeamPointCouplingCondition::ids_in_condition(
    const int id_line, const int id_other) const
{
  if (line_ids_.find(id_line) != line_ids_.end() and line_ids_.find(id_other) != line_ids_.end())
  {
    return true;
  }
  return false;
}

/**
 *
 */
void BEAMINTERACTION::BeamToBeamPointCouplingCondition::clear() {}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToBeamPointCouplingCondition::create_contact_pair(
    const std::vector<Core::Elements::Element const*>& ele_ptrs)
{
  {
    // Check if the given elements are in this condition.
    if (!ids_in_condition(ele_ptrs[0]->id(), ele_ptrs[1]->id())) return Teuchos::null;

    // Create the beam contact pair.
    Teuchos::RCP<BEAMINTERACTION::BeamContactPair> contact_pair = Teuchos::rcp(
        new BeamToBeamPointCouplingPair<GEOMETRYPAIR::t_hermite>(rotational_penalty_parameter_,
            positional_penalty_parameter_, local_parameter_coordinates_));
    // Return the newly created pair.
    return contact_pair;
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamToBeamPointCouplingCondition::build_id_sets(
    const Teuchos::RCP<const Core::FE::Discretization>& discretization)
{
  // Set the IDs of the nodes to be coupled
  const std::vector<int> node_ids = *(condition_line_->get_nodes());

  if (node_ids.size() != 2)
    FOUR_C_THROW(
        "The Penalty Point Coupling Condition can only handle 2 nodes per condition! If you want "
        "to couple multiple nodes, please split them into multiple conditions, each coupling two "
        "of the beam nodes.");

  std::vector<int> element_ids(node_ids.size());
  std::vector<double> position_in_parameter_space(node_ids.size());

  int i = 0;
  for (auto node_id : node_ids)
  {
    i++;
    Core::Nodes::Node* node = discretization->g_node(node_id);
    // This means that the node is not in the column map of this proc and the element pair will thus
    // be created on a different processor
    if (node == nullptr) return;
    Core::Elements::Element* element = node->elements()[0];
    element_ids[i - 1] = element->id();
    if (element->node_ids()[0] == node_id)
      position_in_parameter_space[i - 1] = -1;
    else
      position_in_parameter_space[i - 1] = 1;
  }

  line_ids_ = std::set<int>(element_ids.begin(), element_ids.end());
  local_parameter_coordinates_ = {position_in_parameter_space[0], position_in_parameter_space[1]};
}

FOUR_C_NAMESPACE_CLOSE
