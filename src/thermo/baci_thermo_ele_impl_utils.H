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
  template <CORE::FE::CellType DISTYPE>
  struct DisTypeToOptGaussRule
  {
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::hex8>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::hex18>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_18point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::hex20>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::hex27>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs27>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::tet4>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::tet_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::tet10>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::tet_5point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::wedge6>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::wedge_6point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::pyramid5>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::pyramid_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::quad4>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::quad8>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::quad9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::tri3>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_3point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::tri6>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_6point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::line2>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_2point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::line3>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_3point;
  };

  //! Template Meta Programming version of switch over discretization type
  template <CORE::FE::CellType DISTYPE>
  struct DisTypeToGaussRuleForExactSol
  {
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::hex8>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::hex20>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::hex27>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::nurbs27>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tet4>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tet10>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::wedge6>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::pyramid5>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::quad4>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::quad8>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::quad9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tri3>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tri6>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::line2>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::line3>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
  };

  //! Template Meta Programming version of switch over discretization type
  template <CORE::FE::CellType DISTYPE>
  struct DisTypeToNumGaussPoints
  {
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::hex8>
  {
    static const int nquad = 8;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::hex20>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::hex27>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::nurbs27>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::tet4>
  {
    static const int nquad = 4;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::tet10>
  {
    static const int nquad = 5;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::wedge6>
  {
    static const int nquad = 6;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::pyramid5>
  {
    static const int nquad = 8;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::quad4>
  {
    static const int nquad = 4;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::quad8>
  {
    static const int nquad = 9;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::quad9>
  {
    static const int nquad = 9;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::nurbs9>
  {
    static const int nquad = 9;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::tri3>
  {
    static const int nquad = 3;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::tri6>
  {
    static const int nquad = 6;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::line2>
  {
    static const int nquad = 2;
  };
  template <>
  struct DisTypeToNumGaussPoints<CORE::FE::CellType::line3>
  {
    static const int nquad = 3;
  };

  //! Template Meta Programming version of switch over discretization type
  template <CORE::FE::CellType DISTYPE>
  struct DisTypeToSTRNumGaussPoints
  {
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<CORE::FE::CellType::hex8>
  {
    static const int nquad = 8;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<CORE::FE::CellType::tet4>
  {
    static const int nquad = 5;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<CORE::FE::CellType::tet10>
  {
    static const int nquad = 11;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<CORE::FE::CellType::hex27>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<CORE::FE::CellType::hex20>
  {
    static const int nquad = 27;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<CORE::FE::CellType::hex18>
  {
    static const int nquad = 18;
  };
  template <>
  struct DisTypeToSTRNumGaussPoints<CORE::FE::CellType::nurbs27>
  {
    static const int nquad = 27;
  };

}  // namespace THR

/*----------------------------------------------------------------------*/
#endif  // THERMO_ELE_IMPL_UTILS_H
