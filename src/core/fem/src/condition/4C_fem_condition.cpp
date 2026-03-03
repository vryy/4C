// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_fem_condition.hpp"

#include "4C_fem_general_element.hpp"

FOUR_C_NAMESPACE_OPEN


Core::Conditions::Condition::Condition(const int id, const Core::Conditions::ConditionType type,
    const bool buildgeometry, const Core::Conditions::GeometryType gtype,
    const EntityType entity_type, const std::vector<int>& nodes,
    IO::InputParameterContainer parameters, const std::optional<std::string>& node_set_name)
    : id_(id),
      buildgeometry_(buildgeometry),
      type_(type),
      gtype_(gtype),
      entity_type_(entity_type),
      node_set_name_(node_set_name),
      container_(std::move(parameters))
{
  if (entity_type_ == EntityType::node_set_name)
  {
    FOUR_C_ASSERT_ALWAYS(node_set_name_.has_value(),
        "Condition {} of ENTITY_TYPE 'node_set_name' requires a NODE_SET_NAME to be provided.",
        id_);
  }

  set_nodes(nodes);
}

std::ostream& operator<<(std::ostream& os, const Core::Conditions::Condition& cond)
{
  cond.print(os);
  return os;
}


void Core::Conditions::Condition::print(std::ostream& os) const
{
  os << "Condition " << (node_set_name_.has_value() ? node_set_name_.value() : std::to_string(id_))
     << " " << to_string(type_) << ": ";
  container_.print(os);
  os << std::endl;
  if (nodes_.size() != 0)
  {
    os << "Nodes of this condition:";
    for (const auto& node_gid : nodes_) os << " " << node_gid;
    os << std::endl;
  }
  if (!geometry_.empty())
  {
    os << "Elements of this condition:";
    for (const auto& [ele_id, ele] : geometry_) os << " " << ele_id;
    os << std::endl;
  }
}

std::unique_ptr<Core::Conditions::Condition> Core::Conditions::Condition::copy_without_geometry()
    const
{
  std::unique_ptr<Core::Conditions::Condition> copy(new Condition(*this));
  copy->clear_geometry();
  return copy;
}

const std::string& Core::Conditions::Condition::node_set_name() const
{
  FOUR_C_ASSERT(node_set_name_.has_value(),
      "Condition of ENTITY_TYPE '{}' has no node_set_name assigned!", entity_type_);
  return node_set_name_.value();
}

std::unique_ptr<Core::Conditions::Condition> Core::Conditions::make_condition(
    const Core::Conditions::ConditionSpec& condition_spec,
    const Core::Conditions::InputNodeSets& input_node_sets)
{
  switch (condition_spec.entity_type)
  {
    case Core::Conditions::EntityType::legacy_id:
    {
      FOUR_C_ASSERT_ALWAYS(condition_spec.id.has_value(), "Legacy ID condition must have an ID.");
      const int id = condition_spec.id.value() -
                     1;  // legacy IDs in the input file are 1-based, we convert to 0-based here

      if (input_node_sets.dnode_fenode.size() == 0 && input_node_sets.dline_fenode.size() == 0 &&
          input_node_sets.dsurf_fenode.size() == 0 && input_node_sets.dvol_fenode.size() == 0)
      {
        FOUR_C_THROW(
            "{} condition {} uses legacy_id entity type but no legacy entities were defined in "
            "the input file.\n"
            "This is probably because the geometry is handled in an external file.\n"
            "If this is the case, you must specify a specific entity type (node_set_id or "
            "element_block_id) or identify the node set via its name.\n",
            condition_spec.entity_type, id);
      }
      switch (condition_spec.geometry_type)
      {
        case Core::Conditions::geometry_type_point:
          if (id < 0 or static_cast<unsigned>(id) >= input_node_sets.dnode_fenode.size())
          {
            FOUR_C_THROW(
                "DPoint {} not in range [0:{}[\n"
                "DPoint condition on non existent DPoint?"
                "Could not read set from entity type.",
                id, input_node_sets.dnode_fenode.size());
          }
          return std::make_unique<Core::Conditions::Condition>(id, condition_spec.condition_type,
              condition_spec.build_geometry, condition_spec.geometry_type,
              condition_spec.entity_type, input_node_sets.dnode_fenode[id],
              condition_spec.condition_data);
        case Core::Conditions::geometry_type_line:
          if (id < 0 or static_cast<unsigned>(id) >= input_node_sets.dline_fenode.size())
          {
            FOUR_C_THROW(
                "DLine {} not in range [0:{}[\n"
                "DLine condition on non existent DLine?"
                "Could not read set from entity type.",
                id, input_node_sets.dline_fenode.size());
          }
          return std::make_unique<Core::Conditions::Condition>(id, condition_spec.condition_type,
              condition_spec.build_geometry, condition_spec.geometry_type,
              condition_spec.entity_type, input_node_sets.dline_fenode[id],
              condition_spec.condition_data);
        case Core::Conditions::geometry_type_surface:
          if (id < 0 or static_cast<unsigned>(id) >= input_node_sets.dsurf_fenode.size())
          {
            FOUR_C_THROW(
                "DSurface {} not in range [0:{}[\n"
                "DSurface condition on non existent DSurface?"
                "Could not read set from entity type.",
                id, input_node_sets.dsurf_fenode.size());
          }
          return std::make_unique<Core::Conditions::Condition>(id, condition_spec.condition_type,
              condition_spec.build_geometry, condition_spec.geometry_type,
              condition_spec.entity_type, input_node_sets.dsurf_fenode[id],
              condition_spec.condition_data);
        case Core::Conditions::geometry_type_volume:
          if (id < 0 or static_cast<unsigned>(id) >= input_node_sets.dvol_fenode.size())
          {
            FOUR_C_THROW(
                "DVolume {} not in range [0:{}[\n"
                "DVolume condition on non existent DVolume?",
                id, input_node_sets.dvol_fenode.size());
          }
          return std::make_unique<Core::Conditions::Condition>(id, condition_spec.condition_type,
              condition_spec.build_geometry, condition_spec.geometry_type,
              condition_spec.entity_type, input_node_sets.dvol_fenode[id],
              condition_spec.condition_data);
          break;
        default:
          FOUR_C_THROW("geometry type unspecified");
      }
    }
    case Core::Conditions::EntityType::node_set_id:
    {
      FOUR_C_ASSERT_ALWAYS(condition_spec.id.has_value(), "Node set ID condition must have an ID.");
      const int id = condition_spec.id.value();
      FOUR_C_ASSERT_ALWAYS(input_node_sets.node_sets.contains(id),
          "Cannot apply condition '{}' to node set {} which is not specified in the mesh file.",
          condition_spec.condition_type, id);
      return std::make_unique<Core::Conditions::Condition>(id, condition_spec.condition_type,
          condition_spec.build_geometry, condition_spec.geometry_type, condition_spec.entity_type,
          input_node_sets.node_sets.at(id), condition_spec.condition_data);
    }
    case Core::Conditions::EntityType::node_set_name:
    {
      FOUR_C_ASSERT_ALWAYS(condition_spec.node_set_name.has_value(),
          "Node set name condition must have a node set name.");
      const std::string& node_set_name = condition_spec.node_set_name.value();
      FOUR_C_ASSERT_ALWAYS(input_node_sets.node_sets_names.contains(node_set_name),
          "NODE_SET_NAME '{}' could not be found in the meshfile.", node_set_name);

      const auto& ids = input_node_sets.node_sets_names.at(node_set_name);
      FOUR_C_ASSERT_ALWAYS(ids.size() == 1,
          "NODE_SET_NAME '{}' is not unique in the meshfile ({} occurrences).", node_set_name,
          ids.size());

      return std::make_unique<Core::Conditions::Condition>(ids[0], condition_spec.condition_type,
          condition_spec.build_geometry, condition_spec.geometry_type, condition_spec.entity_type,
          input_node_sets.node_sets.at(ids[0]), condition_spec.condition_data, node_set_name);
    }
    case Core::Conditions::EntityType::element_block_id:
    {
      FOUR_C_ASSERT_ALWAYS(
          condition_spec.id.has_value(), "Element block ID condition must have an ID.");
      const int id = condition_spec.id.value();
      FOUR_C_ASSERT_ALWAYS(input_node_sets.element_block_nodes.contains(id),
          "Cannot apply condition '{}' to element block {} which is not specified in the mesh "
          "file.",
          condition_spec.condition_type, id);
      return std::make_unique<Core::Conditions::Condition>(id, condition_spec.condition_type,
          condition_spec.build_geometry, condition_spec.geometry_type, condition_spec.entity_type,
          input_node_sets.element_block_nodes.at(id), condition_spec.condition_data);
    }
    default:
      FOUR_C_THROW("Unsupported entity type '{}'", condition_spec.entity_type);
  }
}

FOUR_C_NAMESPACE_CLOSE
