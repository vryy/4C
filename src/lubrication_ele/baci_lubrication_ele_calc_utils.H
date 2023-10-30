/*--------------------------------------------------------------------------*/
/*! \file

\brief integration rules

\level 3

*/
/*--------------------------------------------------------------------------*/
#ifndef BACI_LUBRICATION_ELE_CALC_UTILS_H
#define BACI_LUBRICATION_ELE_CALC_UTILS_H

#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_lib_element.H"

namespace LUBRICATION
{
  //! Template Meta Programming version of switch over discretization type
  template <CORE::FE::CellType DISTYPE>
  struct DisTypeToOptGaussRule
  {
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::quad4>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::quad8>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_25point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::quad9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_25point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::tri3>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_7point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::tri6>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_16point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::line2>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_3point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::line3>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_5point;
  };

  //! Template Meta Programming version of switch over discretization type
  template <CORE::FE::CellType DISTYPE>
  struct DisTypeToGaussRuleForExactSol
  {
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::quad4>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::quad8>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_25point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::quad9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_25point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tri3>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_7point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tri6>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_16point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::line2>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_3point;
  };  // not tested
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::line3>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_5point;
  };  // not tested


}  // namespace LUBRICATION


#endif  // LUBRICATION_ELE_CALC_UTILS_H
