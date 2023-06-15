/*--------------------------------------------------------------------------*/
/*! \file

\brief integration rules

\level 3

*/
/*--------------------------------------------------------------------------*/
#ifndef LUBRICATION_ELE_CALC_UTILS_H
#define LUBRICATION_ELE_CALC_UTILS_H

#include "lib_element.H"
#include "discretization_fem_general_utils_integration.H"

namespace LUBRICATION
{
  //! Template Meta Programming version of switch over discretization type
  template <DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToOptGaussRule
  {
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::quad4>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::quad8>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_25point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::quad9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_25point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::tri3>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_7point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::tri6>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_16point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::line2>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_3point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::line3>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_5point;
  };

  //! Template Meta Programming version of switch over discretization type
  template <DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToGaussRuleForExactSol
  {
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::quad4>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::quad8>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_25point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::quad9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_25point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::tri3>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_7point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::tri6>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_16point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::line2>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_3point;
  };  // not tested
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::line3>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_5point;
  };  // not tested


}  // namespace LUBRICATION


#endif  // LUBRICATION_ELE_CALC_UTILS_H
