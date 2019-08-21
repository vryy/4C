/*--------------------------------------------------------------------------*/
/*! \file

\brief integration rules

\level 3

\maintainer Mostafa Faraji
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_CALC_UTILS_H_
#define SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_CALC_UTILS_H_


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
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::quad8>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_25point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::quad9>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_25point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::tri3>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_tri_7point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::tri6>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_tri_16point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::line2>
  {
    static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule_line_3point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::line3>
  {
    static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule_line_5point;
  };

  //! Template Meta Programming version of switch over discretization type
  template <DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToGaussRuleForExactSol
  {
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::quad4>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::quad8>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_25point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::quad9>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_25point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::tri3>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_tri_7point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::tri6>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_tri_16point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::line2>
  {
    static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule_line_3point;
  };  // not tested
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::line3>
  {
    static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule_line_5point;
  };  // not tested


}  // namespace LUBRICATION


#endif /* SRC_DRT_LUBRICATION_ELE_LUBRICATION_ELE_CALC_UTILS_H_ */
