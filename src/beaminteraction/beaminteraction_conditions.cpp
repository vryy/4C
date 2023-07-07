/*----------------------------------------------------------------------*/
/*! \file

\brief Base classes to manage the beam interactions defined by conditions.

\level 3
*/


#include "beaminteraction_conditions.H"

#include "beaminteraction_beam_to_beam_contact_condition.H"
#include "beaminteraction_beam_to_beam_point_coupling_pair.H"
#include "beaminteraction_beam_to_beam_point_coupling_pair_condition.H"
#include "beaminteraction_beam_to_solid_conditions.H"
#include "beaminteraction_contact_params.H"
#include "beaminteraction_contact_pair.H"
#include "beaminteraction_beam_to_solid_volume_meshtying_params.H"
#include "beaminteraction_beam_to_solid_surface_meshtying_params.H"
#include "beaminteraction_beam_to_solid_surface_contact_params.H"

#include "inpar_beam_to_solid.H"
#include "lib_discret.H"
#include "lib_condition.H"
#include "geometry_pair_element.H"
#include "geometry_pair_evaluation_data_base.H"
#include "so3_base.H"


/**
 *
 */
BEAMINTERACTION::BeamInteractionConditionBase::BeamInteractionConditionBase(
    const Teuchos::RCP<const DRT::Condition>& condition_line)
    : condition_line_(condition_line), line_ids_()
{
}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditionBase::BuildIdSets(
    const Teuchos::RCP<const ::DRT::Discretization>& discretization)
{
  // Set the IDs of the line elements.
  std::vector<int> line_ids;
  ConditionToElementIds(condition_line_, line_ids);
  line_ids_ = std::set<int>(line_ids.begin(), line_ids.end());
}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditionBase::Setup(
    const Teuchos::RCP<const ::DRT::Discretization>& discret)
{
}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditionBase::Clear() {}

/**
 *
 */
BEAMINTERACTION::BeamInteractionConditions::BeamInteractionConditions() {}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditions::SetBeamInteractionConditions(
    const Teuchos::RCP<const ::DRT::Discretization>& discret,
    const Teuchos::RCP<const BeamContactParams>& params_ptr)
{
  condition_map_.clear();

  // Get all available interaction types.
  std::vector<INPAR::BEAMINTERACTION::BeamInteractionConditions> interaction_types;
  INPAR::BEAMINTERACTION::BeamInteractionConditionsGetAll(interaction_types);

  // Loop over interaction types.
  for (const auto& interaction_type : interaction_types)
  {
    if (interaction_type == INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_beam_contact)
    {
      // Add all beam-to-beam contitions.
      std::vector<Teuchos::RCP<BeamInteractionConditionBase>>& interaction_vector =
          condition_map_[interaction_type];

      // Get the names for the conditions of this type.
      std::string condition_name = "BeamToBeamContact";

      // Get the line conditions from the discretization.
      std::vector<Teuchos::RCP<DRT::Condition>> condition_lines;
      discret->GetCondition(condition_name, condition_lines);

      // Match the coupling IDs from the input line.
      std::map<int,
          std::pair<Teuchos::RCP<const DRT::Condition>, Teuchos::RCP<const DRT::Condition>>>
          coupling_id_map;
      for (const auto& condition : condition_lines)
      {
        const int coupling_id = condition->GetInt("COUPLING_ID");
        auto& condition_1 = coupling_id_map[coupling_id].first;
        auto& condition_2 = coupling_id_map[coupling_id].second;
        if (condition_1 == Teuchos::null)
          condition_1 = condition;
        else if (condition_2 == Teuchos::null)
          condition_2 = condition;
        else
          dserror("There can not be three different beam-to-beam coupling conditions.");
      }

      for (const auto& map_item : coupling_id_map)
      {
        const auto& condition_1 = map_item.second.first;
        const auto& condition_2 = map_item.second.second;
        if (condition_1 != Teuchos::null && condition_2 != Teuchos::null)
        {
          // We found the matching conditions, now create the beam-to-beam condition objects.
          interaction_vector.push_back(Teuchos::rcp(
              new BEAMINTERACTION::BeamToBeamContactCondition(condition_1, condition_2)));
        }
        else
          dserror("Could not find both conditions (%s) for the COUPLING_ID %d",
              condition_name.c_str(), map_item.first);
      }

      // Check that all conditions were added, i.e. that there are no double definitions of
      // COUPLING_ID.
      if (2 * interaction_vector.size() != condition_lines.size())
        dserror("There are multiple definitions of the same COUPLING_ID for %s",
            condition_name.c_str());
    }
    else if (interaction_type == INPAR::BEAMINTERACTION::BeamInteractionConditions::
                                     beam_to_solid_volume_meshtying or
             interaction_type == INPAR::BEAMINTERACTION::BeamInteractionConditions::
                                     beam_to_solid_surface_meshtying or
             interaction_type ==
                 INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_contact)
    {
      // Add all beam-to-solid conditions.
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
            new_condition = Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidConditionVolumeMeshtying(map_item.second.first,
                    map_item.second.second, params_ptr->BeamToSolidVolumeMeshtyingParams()));
          else if (interaction_type == INPAR::BEAMINTERACTION::BeamInteractionConditions::
                                           beam_to_solid_surface_meshtying)
            new_condition =
                Teuchos::rcp(new BEAMINTERACTION::BeamToSolidConditionSurface(map_item.second.first,
                    map_item.second.second, params_ptr->BeamToSolidSurfaceMeshtyingParams(), true));
          else if (interaction_type ==
                   INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_solid_surface_contact)
            new_condition =
                Teuchos::rcp(new BEAMINTERACTION::BeamToSolidConditionSurface(map_item.second.first,
                    map_item.second.second, params_ptr->BeamToSolidSurfaceContactParams(), false));
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
    else if (interaction_type ==
             INPAR::BEAMINTERACTION::BeamInteractionConditions::beam_to_beam_point_coupling)
    {
      std::vector<Teuchos::RCP<BeamInteractionConditionBase>>& interaction_vector =
          condition_map_[interaction_type];

      // Get the conditions from the discretization.
      std::vector<Teuchos::RCP<DRT::Condition>> coupling_conditions;
      discret->GetCondition("PenaltyPointCouplingCondition", coupling_conditions);
      for (const auto& condition : coupling_conditions)
      {
        // We found the matching conditions, now create the beam-to-beam coupling condition object
        Teuchos::RCP<BeamInteractionConditionBase> new_condition;

        new_condition = Teuchos::rcp(new BEAMINTERACTION::BeamToBeamPointCouplingCondition(
            condition, condition->GetDouble("POSITIONAL_PENALTY_PARAMETER"),
            condition->GetDouble("ROTATIONAL_PENALTY_PARAMETER")));

        interaction_vector.push_back(new_condition);
      }
    }
    else
      dserror("Got unexpected interaction type.");
  }
}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditions::BuildIdSets(
    Teuchos::RCP<::DRT::Discretization> discretization)
{
  for (auto const& map_pair : condition_map_)
    for (auto const& condition : map_pair.second) condition->BuildIdSets(discretization);
}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditions::SetState(
    const Teuchos::RCP<const ::DRT::Discretization>& discret,
    const Teuchos::RCP<const STR::MODELEVALUATOR::BeamInteractionDataState>&
        beaminteraction_data_state)
{
  for (auto const& map_pair : condition_map_)
    for (auto const& condition : map_pair.second)
      condition->SetState(discret, beaminteraction_data_state);
}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditions::Setup(
    const Teuchos::RCP<const ::DRT::Discretization>& discret)
{
  for (auto const& map_pair : condition_map_)
    for (auto const& condition : map_pair.second) condition->Setup(discret);
}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditions::Clear()
{
  for (auto const& map_pair : condition_map_)
    for (auto const& condition : map_pair.second) condition->Clear();
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamInteractionConditions::CreateContactPair(
    const std::vector<DRT::Element const*>& ele_ptrs)
{
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> new_pair;
  for (auto& map_pair : condition_map_)
  {
    for (auto& condition : map_pair.second)
    {
      new_pair = condition->CreateContactPair(ele_ptrs);
      if (new_pair != Teuchos::null) return new_pair;
    }
  }

  // Default return value, i.e. the pair was not found in any of the conditions.
  return Teuchos::null;
}

/**
 *
 */
void BEAMINTERACTION::BeamInteractionConditions::CreateIndirectAssemblyManagers(
    const Teuchos::RCP<const ::DRT::Discretization>& discret,
    std::vector<Teuchos::RCP<SUBMODELEVALUATOR::BeamContactAssemblyManager>>& assembly_managers)
{
  Teuchos::RCP<BEAMINTERACTION::SUBMODELEVALUATOR::BeamContactAssemblyManager>
      condition_assembly_manager = Teuchos::null;
  for (auto& map_pair : condition_map_)
  {
    for (auto& condition : map_pair.second)
    {
      condition_assembly_manager = condition->CreateIndirectAssemblyManager(discret);
      if (not(condition_assembly_manager == Teuchos::null))
        assembly_managers.push_back(condition_assembly_manager);
    }
  }
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
