// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_FEM_CONDITION_DEFINITION_HPP
#define FOUR_C_FEM_CONDITION_DEFINITION_HPP

#include "4C_config.hpp"

#include "4C_io_input_parameter_container.hpp"
#include "4C_io_input_spec.hpp"
#include "4C_legacy_enum_definitions_conditions.hpp"

#include <optional>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

namespace Core::FE
{
  class Discretization;
}  // namespace Core::FE

namespace Core::IO
{
  class InputFile;
}

namespace Core::Conditions
{

  /**
   * Which type of entity is the ID of a condition referring to?
   */
  enum class EntityType
  {
    /**
     * Refers to a numbering of nodes specified in the input file directly. Note that an ID is only
     * unique within a GeometryType. The same ID may be used across different GeometryTypes.
     */
    legacy_id,

    /**
     * Refers to a node set ID.
     */
    node_set_id,

    /**
     * Refers to an element block ID.
     */
    element_block_id,

    /**
     * Refers to an external name specified in the mesh file.
     */
    node_set_name
  };

  /// Struct storing the necessary information to create a condition, parsed from the input file.
  struct ConditionSpec
  {
    const EntityType entity_type;
    const std::optional<int>
        id;  //< The ID of the condition, referring to the entity specified by the EntityType. Not
             // set if the condition is identified via a node set name.
    const std::optional<std::string>
        node_set_name;  //< The name of the node set in the mesh file, if the condition is
                        // identified via a node set name. Not set if the condition is identified
                        // via an ID.
    const ConditionType
        condition_type;                //< The type of the condition (e.g. dirichlet, neumann, etc.)
    const GeometryType geometry_type;  //< The type of geometry the condition lives on (e.g. point,
                                       // line, surface, volume)
    const bool build_geometry;  //< Whether this condition requires an explicit geometry description
                                //(elements) to be built.
    const IO::InputParameterContainer
        condition_data;  //< The parameters parsed from the input file for this condition.
  };

  /**
   * @brief Definition of a condition.
   *
   * This class groups all data that is needed to create a Condition. It contains the InputSpec
   * defining the parameters of the condition and information about the associated geometry.
   */
  class ConditionDefinition
  {
   public:
    /// construction of a condition definition
    /*!
      \param sectionname name of input file section
      \param conditionname name of conditions in Core::FE::Discretization
      \param description description of condition type
      \param condtype type of conditions to be build
      \param buildgeometry whether we need conditions elements
      \param gtype type of geometry the condition lives on
     */
    ConditionDefinition(std::string sectionname, std::string conditionname, std::string description,
        Core::Conditions::ConditionType condtype, bool buildgeometry,
        Core::Conditions::GeometryType gtype);

    /**
     * Add an InputSpec @p spec as another component of this condition. The ordering of the
     * components is irrelevant.
     */
    void add_component(Core::IO::InputSpec&& spec);
    void add_component(const Core::IO::InputSpec& spec);

    /// read all conditions from my input file section
    /*!
      \param input the input file
      \param condition_specs vector of the validated condition specifications in the input to be
      filled.
     */
    void read(Core::IO::InputFile& input, std::vector<ConditionSpec>& condition_specs) const;

    /// name of my section in input file
    std::string section_name() const { return sectionname_; }

    /// my condition name
    std::string name() const { return conditionname_; }

    /// my GeometryType
    Core::Conditions::GeometryType geometry_type() const { return gtype_; }

    /// Get the InputSpec for this ConditionDefinition
    [[nodiscard]] Core::IO::InputSpec spec() const;

   private:
    std::string sectionname_;
    std::string conditionname_;
    std::string description_;
    Core::Conditions::ConditionType condtype_;
    bool buildgeometry_;
    Core::Conditions::GeometryType gtype_;

    std::vector<Core::IO::InputSpec> specs_;

    /**
     * @brief Parse the entity specification for a single condition based on the given input
     * container into a valid ConditionSpec. This checks the consistency of the given entity
     * specification and throws an error if the specification is invalid.
     *
     * @param condition_data Input parameter container of a single condition read from the input
     * file.
     * @return ConditionSpec
     */
    [[nodiscard]] ConditionSpec read_condition_data(
        const Core::IO::InputParameterContainer& condition_data) const;
  };

}  // namespace Core::Conditions


FOUR_C_NAMESPACE_CLOSE

#endif
