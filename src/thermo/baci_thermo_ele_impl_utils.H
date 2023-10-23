/*----------------------------------------------------------------------*/
/*! \file

\level 1

*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 08/09 |
 *----------------------------------------------------------------------*/
#ifndef BACI_THERMO_ELE_IMPL_UTILS_H
#define BACI_THERMO_ELE_IMPL_UTILS_H

/*----------------------------------------------------------------------*
 | headers                                                   dano 08/09 |
 *----------------------------------------------------------------------*/
#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_lib_element.H"

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
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::hex8>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::hex18>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_18point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::hex20>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::hex27>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::nurbs27>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::tet4>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::tet_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::tet10>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::tet_5point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::wedge6>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::wedge_6point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::pyramid5>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::pyramid_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::quad4>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::quad8>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::quad9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::nurbs9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::tri3>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_3point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::tri6>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_6point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::line2>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_2point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::line3>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_3point;
  };

  //! Template Meta Programming version of switch over discretization type
  template <DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToGaussRuleForExactSol
  {
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::hex8>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::hex20>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::hex27>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::nurbs27>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::tet4>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::tet10>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::wedge6>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::pyramid5>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::quad4>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::quad8>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::quad9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::tri3>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::tri6>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::line2>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::line3>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
  };

  //! Template Meta Programming version of switch over discretization type
  template <DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToNumGaussPoints
  {
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::hex8>
  {
    static const int nquad = 8;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::hex20>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::hex27>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::nurbs27>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::tet4>
  {
    static const int nquad = 4;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::tet10>
  {
    static const int nquad = 5;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::wedge6>
  {
    static const int nquad = 6;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::pyramid5>
  {
    static const int nquad = 8;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::quad4>
  {
    static const int nquad = 4;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::quad8>
  {
    static const int nquad = 9;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::quad9>
  {
    static const int nquad = 9;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::nurbs9>
  {
    static const int nquad = 9;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::tri3>
  {
    static const int nquad = 3;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::tri6>
  {
    static const int nquad = 6;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::line2>
  {
    static const int nquad = 2;
  };
  template <>
  struct DisTypeToNumGaussPoints<DRT::Element::DiscretizationType::line3>
  {
    static const int nquad = 3;
  };

  //! Template Meta Programming version of switch over discretization type
  template <DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToSTRNumGaussPoints
  {
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::DiscretizationType::hex8>
  {
    static const int nquad = 8;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::DiscretizationType::tet4>
  {
    static const int nquad = 5;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::DiscretizationType::tet10>
  {
    static const int nquad = 11;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::DiscretizationType::hex27>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::DiscretizationType::hex20>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::DiscretizationType::hex18>
  {
    static const int nquad = 18;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<DRT::Element::DiscretizationType::nurbs27>
  {
    static const int nquad = 27;
  };

}  // namespace THR

/*----------------------------------------------------------------------*/
#endif  // THERMO_ELE_IMPL_UTILS_H
