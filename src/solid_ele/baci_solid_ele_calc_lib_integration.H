/*! \file

\brief A library of free functions for a default solid element

\level 1
*/

#ifndef BACI_SOLID_ELE_CALC_LIB_INTEGRATION_H
#define BACI_SOLID_ELE_CALC_LIB_INTEGRATION_H

#include "baci_discretization_fem_general_cell_type.H"
#include "baci_discretization_fem_general_utils_gausspoints.H"
#include "baci_lib_element_integration_select.H"

#include <Teuchos_ParameterList.hpp>

BACI_NAMESPACE_OPEN
namespace DRT::ELEMENTS
{
  /*!
   * @brief Compare two Gauss integration rules for equality
   */
  inline bool CompareGaussIntegration(const CORE::FE::GaussIntegration& integration_a,
      const CORE::FE::GaussIntegration& integration_b)
  {
    // currently this simple check is sufficient as we only use the same type of gauss integrations.
    return integration_a.NumPoints() == integration_b.NumPoints();
  }

  /*!
   * @brief Get the default Gauss integration rules for different Cell types.
   *
   * @note It follows the rules defined in DRT::ELEMENTS::DisTypeToOptGaussRule<celltype>::rule,
   * except for the stiffness matrix of tetrahedral elements.
   *
   */
  /// @{

  template <CORE::FE::CellType celltype>
  constexpr auto GetGaussRuleMassMatrix()
  {
    return DRT::ELEMENTS::DisTypeToOptGaussRule<celltype>::rule;
  }

  template <CORE::FE::CellType celltype>
  constexpr auto GetGaussRuleStiffnessMatrix()
  {
    return DRT::ELEMENTS::DisTypeToOptGaussRule<celltype>::rule;
  }

  template <>
  constexpr auto GetGaussRuleStiffnessMatrix<CORE::FE::CellType::tet10>()
  {
    return CORE::FE::GaussRule3D::tet_4point;
  }

  template <>
  constexpr auto GetGaussRuleStiffnessMatrix<CORE::FE::CellType::tet4>()
  {
    return CORE::FE::GaussRule3D::tet_1point;
  }

  /*!
   * @brief Create a Gauss integration interface from a given Gauss rule type
   */
  template <CORE::FE::CellType celltype, typename GaussRuleType>
  CORE::FE::GaussIntegration CreateGaussIntegration(GaussRuleType rule)
  {
    // setup default integration
    CORE::FE::IntPointsAndWeights<CORE::FE::dim<celltype>> intpoints(rule);

    // format as DRT::UTILS::GaussIntegration
    Teuchos::RCP<CORE::FE::CollectedGaussPoints> gp =
        Teuchos::rcp(new CORE::FE::CollectedGaussPoints);

    std::array<double, 3> xi = {0., 0., 0.};
    for (int i = 0; i < intpoints.IP().nquad; ++i)
    {
      for (int d = 0; d < CORE::FE::dim<celltype>; ++d) xi[d] = intpoints.IP().qxg[i][d];
      gp->Append(xi[0], xi[1], xi[2], intpoints.IP().qwgt[i]);
    }

    return CORE::FE::GaussIntegration(gp);
  }
  /// @}
}  // namespace DRT::ELEMENTS
BACI_NAMESPACE_CLOSE

#endif  // BACI_SOLID_ELE_CALC_LIB_INTEGRATION_H