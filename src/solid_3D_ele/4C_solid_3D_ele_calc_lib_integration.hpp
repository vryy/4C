// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SOLID_3D_ELE_CALC_LIB_INTEGRATION_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_LIB_INTEGRATION_HPP

#include "4C_config.hpp"

#include "4C_fem_general_cell_type.hpp"
#include "4C_fem_general_cell_type_traits.hpp"
#include "4C_fem_general_element_integration_select.hpp"
#include "4C_fem_general_utils_gausspoints.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_utils_parameter_list.fwd.hpp"

FOUR_C_NAMESPACE_OPEN
namespace Discret::Elements
{
  /*!
   * @brief Compare two Gauss integration rules for equality
   */
  inline bool compare_gauss_integration(const Core::FE::GaussIntegration& integration_a,
      const Core::FE::GaussIntegration& integration_b)
  {
    // currently this simple check is sufficient as we only use the same type of gauss integrations.
    return integration_a.num_points() == integration_b.num_points();
  }

  /*!
   * @brief Get the default Gauss integration rules for different Cell types.
   *
   * @note It follows the rules defined in Discret::Elements::DisTypeToOptGaussRule<celltype>::rule,
   * except for the stiffness matrix of tetrahedral elements.
   *
   */
  /// @{

  template <Core::FE::CellType celltype>
  constexpr auto get_gauss_rule_mass_matrix()
  {
    return Discret::Elements::DisTypeToOptGaussRule<celltype>::rule;
  }

  template <Core::FE::CellType celltype>
  constexpr auto get_gauss_rule_stiffness_matrix()
  {
    return Discret::Elements::DisTypeToOptGaussRule<celltype>::rule;
  }

  template <>
  constexpr auto get_gauss_rule_stiffness_matrix<Core::FE::CellType::tet10>()
  {
    return Core::FE::GaussRule3D::tet_4point;
  }

  template <>
  constexpr auto get_gauss_rule_stiffness_matrix<Core::FE::CellType::tet4>()
  {
    return Core::FE::GaussRule3D::tet_1point;
  }
  /// @}

  namespace Internal
  {
    template <Core::FE::CellType celltype>
    struct ApplicableIntegrationRules;

    template <Core::FE::CellType celltype>
      requires Core::FE::is_hex<celltype>
    struct ApplicableIntegrationRules<celltype>
    {
      static constexpr std::array value = {
          Core::FE::GaussRule3D::hex_1point,
          Core::FE::GaussRule3D::hex_8point,
          Core::FE::GaussRule3D::hex_18point,
          Core::FE::GaussRule3D::hex_27point,
          Core::FE::GaussRule3D::hex_64point,
          Core::FE::GaussRule3D::hex_125point,
          Core::FE::GaussRule3D::hex_216point,
          Core::FE::GaussRule3D::hex_343point,
          Core::FE::GaussRule3D::hex_512point,
          Core::FE::GaussRule3D::hex_729point,
          Core::FE::GaussRule3D::hex_1000point,
      };
    };

    template <Core::FE::CellType celltype>
      requires Core::FE::is_tet<celltype>
    struct ApplicableIntegrationRules<celltype>
    {
      static constexpr std::array value = {
          Core::FE::GaussRule3D::tet_1point,
          Core::FE::GaussRule3D::tet_4point,
          Core::FE::GaussRule3D::tet_5point,
          Core::FE::GaussRule3D::tet_11point,
          Core::FE::GaussRule3D::tet_15point,
          Core::FE::GaussRule3D::tet_24point,
          Core::FE::GaussRule3D::tet_45point,
      };
    };

    template <Core::FE::CellType celltype>
      requires Core::FE::is_pyramid<celltype>
    struct ApplicableIntegrationRules<celltype>
    {
      static constexpr std::array value = {
          Core::FE::GaussRule3D::pyramid_1point,
          Core::FE::GaussRule3D::pyramid_8point,
      };
    };

    template <Core::FE::CellType celltype>
      requires Core::FE::is_wedge<celltype>
    struct ApplicableIntegrationRules<celltype>
    {
      static constexpr std::array value = {
          Core::FE::GaussRule3D::wedge_1point,
          Core::FE::GaussRule3D::wedge_6point,
          Core::FE::GaussRule3D::wedge_9point,
      };
    };

    template <Core::FE::CellType celltype>
      requires(celltype == Core::FE::CellType::nurbs27)
    struct ApplicableIntegrationRules<celltype>
    {
      static constexpr std::array value =
          ApplicableIntegrationRules<Core::FE::CellType::hex8>::value;
    };
  }  // namespace Internal

  template <Core::FE::CellType celltype>
  static constexpr std::array applicable_integration_rules =
      Internal::ApplicableIntegrationRules<celltype>::value;

  template <unsigned dim>
  struct SolidIntegrationRules;

  template <>
  struct SolidIntegrationRules<3>
  {
    Core::FE::GaussRule3D rule_residuum;
    Core::FE::GaussRule3D rule_mass;
  };

  template <Core::FE::CellType celltype>
  SolidIntegrationRules<3> make_default_solid_integration_rules()
  {
    return {
        .rule_residuum = get_gauss_rule_stiffness_matrix<celltype>(),
        .rule_mass = get_gauss_rule_mass_matrix<celltype>(),
    };
  }
}  // namespace Discret::Elements
FOUR_C_NAMESPACE_CLOSE

#endif