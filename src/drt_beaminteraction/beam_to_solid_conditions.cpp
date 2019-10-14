/*----------------------------------------------------------------------*/
/*! \file

\brief Manage the beam-to-solid conditions.

\level 3
\maintainer Ivo Steinbrecher
*/


#include "beam_to_solid_conditions.H"

#include "beam_contact_params.H"
#include "beam_contact_pair.H"
#include "beam_to_solid_volume_meshtying_pair_gauss_point.H"
#include "beam_to_solid_volume_meshtying_pair_mortar.H"
#include "beam_to_solid_volume_meshtying_pair_gauss_point_cross_section.H"
#include "beam_to_solid_volume_meshtying_params.H"

#include "../drt_inpar/inpar_beam_to_solid.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_condition.H"
#include "../drt_geometry_pair/geometry_pair_element.H"
#include "../drt_so3/so_base.H"


/**
 *
 */
BEAMINTERACTION::BeamToSolidCondition::BeamToSolidCondition(
    const Teuchos::RCP<const DRT::Condition>& condition_line,
    const Teuchos::RCP<const DRT::Condition>& condition_other)
    : condition_line_(condition_line), condition_other_(condition_other)
{
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidCondition::BuildIdSets()
{
  // Set the IDs of the line elements.
  std::vector<int> line_ids;
  ConditionToElementIds(condition_line_, line_ids);
  line_ids_ = std::set<int>(line_ids.begin(), line_ids.end());
}

/**
 *
 */
bool BEAMINTERACTION::BeamToSolidCondition::IdsInCondition(
    const int id_line, const int id_other) const
{
  if (line_ids_.find(id_line) != line_ids_.end())
    if (IdInOther(id_other)) return true;
  return false;
}

/**
 *
 */
BEAMINTERACTION::BeamToSolidConditionVolumeMeshtying::BeamToSolidConditionVolumeMeshtying(
    const Teuchos::RCP<const DRT::Condition>& condition_line,
    const Teuchos::RCP<const DRT::Condition>& condition_other)
    : BeamToSolidCondition(condition_line, condition_other)
{
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidConditionVolumeMeshtying::CreateBeamToSolidPair(
    const std::vector<DRT::Element const*>& ele_ptrs,
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr)
{
  // Check if the given elements are in this condition.
  if (!IdsInCondition(ele_ptrs[0]->Id(), ele_ptrs[1]->Id())) return Teuchos::null;

  // Cast the solid element.
  DRT::ELEMENTS::So_base const* solidele = dynamic_cast<DRT::ELEMENTS::So_base const*>(ele_ptrs[1]);
  DRT::Element::DiscretizationType shape = solidele->Shape();

  // Get the contact discretization method.
  INPAR::BEAMTOSOLID::BeamToSolidVolumeContactDiscretization contact_discretization =
      params_ptr->BeamToSolidVolumeMeshtyingParams()->GetContactDiscretization();

  // Check which contact discretization is wanted.
  if (contact_discretization ==
      INPAR::BEAMTOSOLID::BeamToSolidVolumeContactDiscretization::gauss_point_to_segment)
  {
    switch (shape)
    {
      case DRT::Element::hex8:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex8>());
      case DRT::Element::hex20:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex20>());
      case DRT::Element::hex27:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_hex27>());
      case DRT::Element::tet4:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_tet4>());
      case DRT::Element::tet10:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_tet10>());
      case DRT::Element::nurbs27:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPoint<GEOMETRYPAIR::t_hermite,
                GEOMETRYPAIR::t_nurbs27>());

      default:
        dserror("Wrong element type for solid element.");
    }
  }
  else if (contact_discretization ==
           INPAR::BEAMTOSOLID::BeamToSolidVolumeContactDiscretization::mortar)
  {
    INPAR::BEAMTOSOLID::BeamToSolidVolumeMortarShapefunctions mortar_shape_function =
        params_ptr->BeamToSolidVolumeMeshtyingParams()->GetMortarShapeFunctionType();

    switch (mortar_shape_function)
    {
      case INPAR::BEAMTOSOLID::BeamToSolidVolumeMortarShapefunctions::line2:
      {
        switch (shape)
        {
          case DRT::Element::hex8:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line2>());
          case DRT::Element::hex20:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line2>());
          case DRT::Element::hex27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line2>());
          case DRT::Element::tet4:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line2>());
          case DRT::Element::tet10:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line2>());
          case DRT::Element::nurbs27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_nurbs27, GEOMETRYPAIR::t_line2>());
          default:
            dserror("Wrong element type for solid element.");
        }
        break;
      }
      case INPAR::BEAMTOSOLID::BeamToSolidVolumeMortarShapefunctions::line3:
      {
        switch (shape)
        {
          case DRT::Element::hex8:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line3>());
          case DRT::Element::hex20:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line3>());
          case DRT::Element::hex27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line3>());
          case DRT::Element::tet4:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line3>());
          case DRT::Element::tet10:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line3>());
          case DRT::Element::nurbs27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_nurbs27, GEOMETRYPAIR::t_line3>());
          default:
            dserror("Wrong element type for solid element.");
        }
        break;
      }
      case INPAR::BEAMTOSOLID::BeamToSolidVolumeMortarShapefunctions::line4:
      {
        switch (shape)
        {
          case DRT::Element::hex8:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex8, GEOMETRYPAIR::t_line4>());
          case DRT::Element::hex20:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex20, GEOMETRYPAIR::t_line4>());
          case DRT::Element::hex27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_hex27, GEOMETRYPAIR::t_line4>());
          case DRT::Element::tet4:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet4, GEOMETRYPAIR::t_line4>());
          case DRT::Element::tet10:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_tet10, GEOMETRYPAIR::t_line4>());
          case DRT::Element::nurbs27:
            return Teuchos::rcp(
                new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairMortar<GEOMETRYPAIR::t_hermite,
                    GEOMETRYPAIR::t_nurbs27, GEOMETRYPAIR::t_line4>());
          default:
            dserror("Wrong element type for solid element.");
        }
        break;
      }
      default:
        dserror("Wrong mortar shape function.");
    }
  }
  if (contact_discretization ==
      INPAR::BEAMTOSOLID::BeamToSolidVolumeContactDiscretization::gauss_point_cross_section)
  {
    switch (shape)
    {
      case DRT::Element::hex8:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex8>());
      case DRT::Element::hex20:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex20>());
      case DRT::Element::hex27:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_hex27>());
      case DRT::Element::tet4:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet4>());
      case DRT::Element::tet10:
        return Teuchos::rcp(
            new BEAMINTERACTION::BeamToSolidVolumeMeshtyingPairGaussPointCrossSection<
                GEOMETRYPAIR::t_hermite, GEOMETRYPAIR::t_tet10>());
      default:
        dserror("Wrong element type for solid element.");
    }
  }

  // Default return value.
  return Teuchos::null;
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidConditionVolumeMeshtying::BuildIdSets()
{
  // Call the parent method to build the line maps.
  BeamToSolidCondition::BuildIdSets();

  // Build the volume map.
  std::vector<int> volume_ids;
  ConditionToElementIds(condition_other_, volume_ids);
  volume_ids_ = std::set<int>(volume_ids.begin(), volume_ids.end());
}

/**
 *
 */
BEAMINTERACTION::BeamToSolidConditionSurfaceMeshtying::BeamToSolidConditionSurfaceMeshtying(
    const Teuchos::RCP<const DRT::Condition>& condition_line,
    const Teuchos::RCP<const DRT::Condition>& condition_other)
    : BeamToSolidCondition(condition_line, condition_other)
{
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair>
BEAMINTERACTION::BeamToSolidConditionSurfaceMeshtying::CreateBeamToSolidPair(
    const std::vector<DRT::Element const*>& ele_ptrs,
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr)
{
  dserror("Not yet implemented");
  return Teuchos::null;
}

/**
 *
 */
BEAMINTERACTION::BeamToSolidConditions::BeamToSolidConditions() {}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidConditions::SetBeamToSolidConditions(
    const Teuchos::RCP<const DRT::Discretization>& discret)
{
  condition_map_.clear();

  // Get all available interaction types.
  std::vector<INPAR::BEAMTOSOLID::BeamToSolidInteraction> interaction_types;
  INPAR::BEAMTOSOLID::BeamToSolidInteractionGetAll(interaction_types);

  // Loop over interaction types.
  for (const auto& interaction_type : interaction_types)
  {
    std::vector<Teuchos::RCP<BeamToSolidCondition>>& interaction_vector =
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
    std::map<int, std::pair<Teuchos::RCP<const DRT::Condition>, Teuchos::RCP<const DRT::Condition>>>
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
        Teuchos::RCP<BeamToSolidCondition> new_condition;
        if (interaction_type ==
            INPAR::BEAMTOSOLID::BeamToSolidInteraction::beam_to_solid_volume_meshtying)
          new_condition = Teuchos::rcp(new BeamToSolidConditionVolumeMeshtying(
              map_item.second.first, map_item.second.second));
        else if (interaction_type ==
                 INPAR::BEAMTOSOLID::BeamToSolidInteraction::beam_to_solid_surface_meshtying)
          new_condition = Teuchos::rcp(new BeamToSolidConditionSurfaceMeshtying(
              map_item.second.first, map_item.second.second));
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
}

/**
 *
 */
void BEAMINTERACTION::BeamToSolidConditions::BuildIdSets()
{
  for (auto const& map_pair : condition_map_)
    for (auto const& condition : map_pair.second) condition->BuildIdSets();
}

/**
 *
 */
Teuchos::RCP<BEAMINTERACTION::BeamContactPair> BEAMINTERACTION::BeamToSolidConditions::CreatePair(
    const std::vector<DRT::Element const*>& ele_ptrs,
    const Teuchos::RCP<BEAMINTERACTION::BeamContactParams>& params_ptr)
{
  Teuchos::RCP<BEAMINTERACTION::BeamContactPair> new_pair;
  for (auto& map_pair : condition_map_)
  {
    for (auto& condition : map_pair.second)
    {
      new_pair = condition->CreateBeamToSolidPair(ele_ptrs, params_ptr);
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
