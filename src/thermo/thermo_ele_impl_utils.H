/*----------------------------------------------------------------------*/
/*! \file

\level 1

*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 08/09 |
 *----------------------------------------------------------------------*/
#ifndef THERMO_ELE_IMPL_UTILS_H
#define THERMO_ELE_IMPL_UTILS_H

/*----------------------------------------------------------------------*
 | headers                                                   dano 08/09 |
 *----------------------------------------------------------------------*/
#include "discretization_fem_general_utils_integration.H"
#include "lib_element.H"

/*----------------------------------------------------------------------*
 |                                                           dano 08/09 |
 *----------------------------------------------------------------------*/
namespace THR
{
  //! Template Meta Programming version of switch over discretization type
  template <DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToOptGaussRule
  {
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::hex8>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::hex_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::hex18>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::hex_18point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::hex20>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::hex27>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::nurbs27>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::tet4>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::tet_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::tet10>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::tet_5point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::wedge6>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::wedge_6point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::pyramid5>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::pyramid_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::quad4>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::GaussRule2D::quad_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::quad8>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::quad9>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::nurbs9>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::tri3>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::GaussRule2D::tri_3point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::tri6>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::GaussRule2D::tri_6point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::line2>
  {
    static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::GaussRule1D::line_2point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::line3>
  {
    static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::GaussRule1D::line_3point;
  };

  //! Template Meta Programming version of switch over discretization type
  template <DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToGaussRuleForExactSol
  {
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::hex8>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::hex20>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::hex27>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs27>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::tet4>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::tet10>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::wedge6>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::pyramid5>
  {
    static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::quad4>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::quad8>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::quad9>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::tri3>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::GaussRule2D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::tri6>
  {
    static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::GaussRule2D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::line2>
  {
    static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::GaussRule1D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::line3>
  {
    static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::GaussRule1D::undefined;
  };

  //! Template Meta Programming version of switch over discretization type
  template <DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToNumGaussPoints
  {
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::hex8>
  {
    static const int nquad = 8;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::hex20>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::hex27>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::nurbs27>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::tet4>
  {
    static const int nquad = 4;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::tet10>
  {
    static const int nquad = 5;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::wedge6>
  {
    static const int nquad = 6;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::pyramid5>
  {
    static const int nquad = 8;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::quad4>
  {
    static const int nquad = 4;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::quad8>
  {
    static const int nquad = 9;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::quad9>
  {
    static const int nquad = 9;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::nurbs9>
  {
    static const int nquad = 9;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::tri3>
  {
    static const int nquad = 3;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::tri6>
  {
    static const int nquad = 6;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::line2>
  {
    static const int nquad = 2;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::line3>
  {
    static const int nquad = 3;
  };

  //! Template Meta Programming version of switch over discretization type
  template <DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToSTRNumGaussPoints
  {
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::hex8>
  {
    static const int nquad = 8;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::tet4>
  {
    static const int nquad = 5;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::tet10>
  {
    static const int nquad = 11;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::hex27>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::hex20>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::hex18>
  {
    static const int nquad = 18;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::nurbs27>
  {
    static const int nquad = 27;
  };

}  // namespace THR

/*----------------------------------------------------------------------*/
#endif  // THERMO_ELE_IMPL_UTILS_H
