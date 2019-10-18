/*----------------------------------------------------------------------*/
/*! \file

\brief Base classes to manage the beam interactions defined by conditions.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beaminteraction_conditions.H"

#include "beam_to_solid_conditions.H"
#include "beam_contact_params.H"
#include "beam_contact_pair.H"

#include "../drt_inpar/inpar_beam_to_solid.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_geometry_pair/geometry_pair_element.H"
#include "../drt_geometry_pair/geometry_pair_evaluation_data_base.H"
#include "../drt_so3/so_base.H"


/**
 *
 */
BEAMINTERACTION::BeamInteractionConditionBase::BeamInteractionConditionBase(
    const Teuchos::RCP<const DRT::Condition>& condition_line)
    : condition_line_(condition_line), line_ids_(), geometry_evaluation_data_(Teuchos::null)
{
}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditionBase::BuildIdSets()
{
  // Set the IDs of the line elements.
  std::vector<int> line_ids;
  ConditionToElementIds(condition_line_, line_ids);
  line_ids_ = std::set<int>(line_ids.begin(), line_ids.end());
}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditionBase::Reset()
{
  // Reset the geometry evaluation tracker.
  geometry_evaluation_data_->Reset();
}

/**
 *
 */
BEAMINTERACTION::BeamInteractionConditions::BeamInteractionConditions() {}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditions::SetBeamInteractionConditions(
    const Teuchos::RCP<const DRT::Discretization>& discret)
{
  condition_map_.clear();

  // Get all available interaction types.
  std::vector<INPAR::BEAMINTERACTION::BeamInteractionConditions> interaction_types;
  INPAR::BEAMINTERACTION::BeamInteractionConditionsGetAll(interaction_types);

  // Loop over interaction types.
  for (const auto& interaction_type : interaction_types)
  {
    // Add all beam-to-solid contitions.
    if (interaction_type ==
            INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_volume_meshtying ||
        interaction_type ==
            INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_meshtying)
    {
      std::vector<Teuchos::RCP<BeamInteractionConditionBase>>& interaction_vector =
          condition_map_[interaction_type];

      // Get the names for the conditions of this type.
      std::array<std::string, 2> condition_names;
      INPAR::BEAMTOSOLID::BeamToSolidInteractionGetString(interaction_type, condition_names);

      // Get the conditions from the discretization.
      std::vector<Teuchos::RCP<DRT::Condition>> condition_line;
      std::vector<Teuchos::RCP<DRT::Condition>> condition_other;
      discret->GetCondition(condition_names[0], condition_line);
      discret->GetCondition(condition_names[1], condition_other);

      // There has to be an equal number of sections for lines and surfaces / volumes.
      if (condition_line.size() != condition_other.size())
        dserror("There are %d %s sections and %d %s sections. The numbers have to match!",
            condition_line.size(), condition_names[0].c_str(), condition_other.size(),
            condition_names[1].c_str());

      // Match the coupling IDs from the input line.
      std::map<int,
          std::pair<Teuchos::RCP<const DRT::Condition>, Teuchos::RCP<const DRT::Condition>>>
          coupling_id_map;
      for (const auto& condition : condition_line)
        coupling_id_map[condition->GetInt("COUPLING_ID")].first = condition;
      for (const auto& condition : condition_other)
        coupling_id_map[condition->GetInt("COUPLING_ID")].second = condition;
      for (const auto& map_item : coupling_id_map)
      {
        if (map_item.second.first != Teuchos::null && map_item.second.second != Teuchos::null)
        {
          // We found the matching conditions, now create the beam-to-solid condition objects.
          Teuchos::RCP<BeamInteractionConditionBase> new_condition;
          if (interaction_type ==
              INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_volume_meshtying)
            new_condition = Teuchos::rcp(new BEAMINTERACTION::BeamToSolidConditionVolumeMeshtying(
                map_item.second.first, map_item.second.second));
          else if (interaction_type == INPAR::BEAMINTERACTION::BeamInteractionConditions::
                                           beam_to_solid_surface_meshtying)
            new_condition = Teuchos::rcp(new BEAMINTERACTION::BeamToSolidConditionSurfaceMeshtying(
                map_item.second.first, map_item.second.second));
          else
            dserror("Got unexpected interaction type.");
          interaction_vector.push_back(new_condition);
        }
        else
          dserror("Could not find both conditions (%s, %s) for the COUPLING_ID %d",
              condition_names[0].c_str(), condition_names[1].c_str(), map_item.first);
      }

      // Check that all conditions were added, i.e. that there are no double definitions of
      // COUPLING_ID.
      if (interaction_vector.size() != condition_line.size())
        dserror("There are multiple definitions of the same COUPLING_ID for %s and %s",
            condition_names[0].c_str(), condition_names[1].c_str());
    }
    else
      dserror("Got unexpected interaction type.");
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditions::BuildIdSets()
{
  for (auto const& map_pair : condition_map_)
    for (auto const& condition : map_pair.second) condition->BuildIdSets();
}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditions::Reset()
{
  for (auto const& map_pair : condition_map_)
    for (auto const& condition : map_pair.second) condition->Reset();
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamInteractionConditions::CreateContactPair(
    const std::vector<DRT::Element const*>& ele_ptrs,
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr)
{
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> new_pair;
  for (auto& map_pair : condition_map_)
  {
    for (auto& condition : map_pair.second)
    {
      new_pair = condition->CreateContactPair(ele_ptrs, params_ptr);
      if (new_pair != Teuchos::null) return new_pair;
    }
  }

  // Default return value, i.e. the pair was not found in any of the conditions.
  return Teuchos::null;
}

/**
 *
 */
void BEAMINTERACTION::ConditionToElementIds(
    const Teuchos::RCP<const DRT::Condition>& condition, std::vector<int>& element_ids)
{
  // Loop over the elements in the condition and get the "real" element by comparing the node IDs.
  element_ids.clear();
  element_ids.reserve(condition->Geometry().size());
  for (const auto& item : condition->Geometry())
  {
    int n_nodes = item.second->NumNode();

    // Create the node sets and store the node IDs from the condition element in it.
    std::set<int> nodes_condition;
    std::set<int> nodes_element;
    for (int i = 0; i < n_nodes; i++) nodes_condition.insert(item.second->Nodes()[i]->Id());

    // Loop over all elements connected to a node and check if the nodal IDs are the same. Use the
    // last node, since if there are nodes connected to fewer elements, those are usually at the
    // end of the list.
    int local_node_id = n_nodes - 1;
    DRT::Element** elements = item.second->Nodes()[local_node_id]->Elements();
    for (int i_element = 0; i_element < item.second->Nodes()[local_node_id]->NumElement();
         i_element++)
    {
      if (elements[i_element]->NumNode() != n_nodes) continue;

      // Fill up the node ID map.
      nodes_element.clear();
      for (int i_nodes = 0; i_nodes < n_nodes; i_nodes++)
        nodes_element.insert(elements[i_element]->Nodes()[i_nodes]->Id());

      // Check if the maps are equal.
      if (std::equal(nodes_condition.begin(), nodes_condition.end(), nodes_element.begin()))
        element_ids.push_back(elements[i_element]->Id());
    }
  }

  // Check if all elements were found.
  if (condition->Geometry().size() != element_ids.size())
    dserror("Could not find the IDs of all elements!");
}
