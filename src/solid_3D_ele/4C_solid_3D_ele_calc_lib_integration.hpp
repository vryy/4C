/*! \file

\brief A library of free functions for a default solid element

\level 1
*/

#ifndef FOUR_C_SOLID_3D_ELE_CALC_LIB_INTEGRATION_HPP
#define FOUR_C_SOLID_3D_ELE_CALC_LIB_INTEGRATION_HPP

#include "4C_config.hpp"

#include "4C_discretization_fem_general_cell_type.hpp"
#include "4C_discretization_fem_general_element_integration_select.hpp"
#include "4C_discretization_fem_general_utils_gausspoints.hpp"

#include <Teuchos_ParameterList.hpp>

FOUR_C_NAMESPACE_OPEN
namespace Discret::ELEMENTS
{
  /*!
   * @brief Compare two Gauss integration rules for equality
   */
  inline bool CompareGaussIntegration(const Core::FE::GaussIntegration& integration_a,
      const Core::FE::GaussIntegration& integration_b)
  {
    // currently this simple check is sufficient as we only use the same type of gauss integrations.
    return integration_a.NumPoints() == integration_b.NumPoints();
  }

  /*!
   * @brief Get the default Gauss integration rules for different Cell types.
   *
   * @note It follows the rules defined in Discret::ELEMENTS::DisTypeToOptGaussRule<celltype>::rule,
   * except for the stiffness matrix of tetrahedral elements.
   *
   */
  /// @{

  template <Core::FE::CellType celltype>
  constexpr auto GetGaussRuleMassMatrix()
  {
    return Discret::ELEMENTS::DisTypeToOptGaussRule<celltype>::rule;
  }

  template <Core::FE::CellType celltype>
  constexpr auto GetGaussRuleStiffnessMatrix()
  {
    return Discret::ELEMENTS::DisTypeToOptGaussRule<celltype>::rule;
  }

  template <>
  constexpr auto GetGaussRuleStiffnessMatrix<Core::FE::CellType::tet10>()
  {
    return Core::FE::GaussRule3D::tet_4point;
  }

  template <>
  constexpr auto GetGaussRuleStiffnessMatrix<Core::FE::CellType::tet4>()
  {
    return Core::FE::GaussRule3D::tet_1point;
  }

  /*!
   * @brief Create a Gauss integration interface from a given Gauss rule type
   */
  template <Core::FE::CellType celltype, typename GaussRuleType>
  Core::FE::GaussIntegration CreateGaussIntegration(GaussRuleType rule)
  {
    // setup default integration
    Core::FE::IntPointsAndWeights<Core::FE::dim<celltype>> intpoints(rule);

    // format as Discret::UTILS::GaussIntegration
    Teuchos::RCP<Core::FE::CollectedGaussPoints> gp =
        Teuchos::rcp(new Core::FE::CollectedGaussPoints);

    std::array<double, 3> xi = {0., 0., 0.};
    for (int i = 0; i < intpoints.IP().nquad; ++i)
    {
      for (int d = 0; d < Core::FE::dim<celltype>; ++d) xi[d] = intpoints.IP().qxg[i][d];
      gp->Append(xi[0], xi[1], xi[2], intpoints.IP().qwgt[i]);
    }

    return Core::FE::GaussIntegration(gp);
  }
  /// @}
}  // namespace Discret::ELEMENTS
FOUR_C_NAMESPACE_CLOSE

#endif