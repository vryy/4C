/*----------------------------------------------------------------------*/
/*! \file

\brief Utility methods for scatra

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_CALC_UTILS_HPP
#define FOUR_C_SCATRA_ELE_CALC_UTILS_HPP


#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_integration.hpp"
#include "baci_inpar_scatra.hpp"
#include "baci_lib_element.hpp"

FOUR_C_NAMESPACE_OPEN


namespace SCATRA
{
  /*!
  \brief Decide, whether second derivatives are needed  (template version)
   *  In convection-diffusion problems, ONLY N,xx , N,yy and N,zz are needed
   *  to evaluate the laplacian operator for the residual-based stabilization.
   *  Hence, unlike to the Navier-Stokes equations, hex8, wedge6 and pyramid5
   *  return false although they have non-zero MIXED second derivatives.*/
  template <CORE::FE::CellType DISTYPE>
  struct Use2ndDerivs
  {
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::hex8>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::tet4>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::wedge6>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::pyramid5>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::nurbs8>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::quad4>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::nurbs4>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::tri3>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::line2>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::nurbs2>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::hex20>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::hex27>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::nurbs27>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::tet10>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::quad8>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::quad9>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::nurbs9>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::tri6>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::line3>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<CORE::FE::CellType::nurbs3>
  {
    static constexpr bool use = true;
  };

  //! Template Meta Programming version of switch over discretization type
  template <CORE::FE::CellType DISTYPE>
  struct DisTypeToOptGaussRule
  {
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::hex8>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::hex20>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::hex27>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::tet4>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::tet_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::tet10>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::tet_11point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::wedge6>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::wedge_6point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::pyramid5>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::pyramid_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs8>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs27>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::quad4>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::quad8>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::quad9>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::tri3>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::tri_3point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::tri6>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::tri_6point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs4>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs9>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::line2>
  {
    static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::line_2point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::line3>
  {
    static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::line_3point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs2>
  {
    static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::line_2point;
  };
  template <>
  struct DisTypeToOptGaussRule<CORE::FE::CellType::nurbs3>
  {
    static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::line_3point;
  };

  //! Template Meta Programming version of switch over discretization type
  template <CORE::FE::CellType DISTYPE>
  struct DisTypeToMatGaussRule
  {
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::hex8>
  {
    static CORE::FE::GaussRule3D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::FE::GaussRule3D::hex_1point;
          break;
        case 2:
        case 3:
          return CORE::FE::GaussRule3D::hex_8point;
          break;
        case 4:
        case 5:
          return CORE::FE::GaussRule3D::hex_27point;
          break;
        case 6:
        case 7:
          return CORE::FE::GaussRule3D::hex_64point;
          break;
        case 8:
        case 9:
          return CORE::FE::GaussRule3D::hex_125point;
          break;
        case 10:
        case 11:
          return CORE::FE::GaussRule3D::hex_216point;
          break;
        case 12:
        case 13:
          return CORE::FE::GaussRule3D::hex_343point;
          break;
        case 14:
        case 15:
          return CORE::FE::GaussRule3D::hex_512point;
          break;
        case 16:
        case 17:
          return CORE::FE::GaussRule3D::hex_729point;
          break;
        case 18:
        case 19:
          return CORE::FE::GaussRule3D::hex_1000point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for HEX elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::FE::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::hex20>
  {
    static CORE::FE::GaussRule3D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::FE::GaussRule3D::hex_1point;
          break;
        case 2:
        case 3:
          return CORE::FE::GaussRule3D::hex_8point;
          break;
        case 4:
        case 5:
          return CORE::FE::GaussRule3D::hex_27point;
          break;
        case 6:
        case 7:
          return CORE::FE::GaussRule3D::hex_64point;
          break;
        case 8:
        case 9:
          return CORE::FE::GaussRule3D::hex_125point;
          break;
        case 10:
        case 11:
          return CORE::FE::GaussRule3D::hex_216point;
          break;
        case 12:
        case 13:
          return CORE::FE::GaussRule3D::hex_343point;
          break;
        case 14:
        case 15:
          return CORE::FE::GaussRule3D::hex_512point;
          break;
        case 16:
        case 17:
          return CORE::FE::GaussRule3D::hex_729point;
          break;
        case 18:
        case 19:
          return CORE::FE::GaussRule3D::hex_1000point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for HEX elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::FE::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::hex27>
  {
    static CORE::FE::GaussRule3D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::FE::GaussRule3D::hex_1point;
          break;
        case 2:
        case 3:
          return CORE::FE::GaussRule3D::hex_8point;
          break;
        case 4:
        case 5:
          return CORE::FE::GaussRule3D::hex_27point;
          break;
        case 6:
        case 7:
          return CORE::FE::GaussRule3D::hex_64point;
          break;
        case 8:
        case 9:
          return CORE::FE::GaussRule3D::hex_125point;
          break;
        case 10:
        case 11:
          return CORE::FE::GaussRule3D::hex_216point;
          break;
        case 12:
        case 13:
          return CORE::FE::GaussRule3D::hex_343point;
          break;
        case 14:
        case 15:
          return CORE::FE::GaussRule3D::hex_512point;
          break;
        case 16:
        case 17:
          return CORE::FE::GaussRule3D::hex_729point;
          break;
        case 18:
        case 19:
          return CORE::FE::GaussRule3D::hex_1000point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for HEX elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::FE::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::tet4>
  {
    static CORE::FE::GaussRule3D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::FE::GaussRule3D::tet_1point;
          break;
        case 2:
          return CORE::FE::GaussRule3D::tet_4point;
          break;
        case 3:
          return CORE::FE::GaussRule3D::tet_5point;
          break;
        case 4:
          return CORE::FE::GaussRule3D::tet_11point;
          break;
        case 5:
          return CORE::FE::GaussRule3D::tet_15point;
          break;
        case 6:
          return CORE::FE::GaussRule3D::tet_24point;
          break;
        case 7:
        case 8:
          return CORE::FE::GaussRule3D::tet_45point;
          break;
        case 9:
          return CORE::FE::GaussRule3D::tet_125point_peano;
          break;
        case 10:
        case 11:
        case 12:
        case 13:
          return CORE::FE::GaussRule3D::tet_343point_peano;
          break;
        case 14:
        case 15:
        case 16:
          return CORE::FE::GaussRule3D::tet_729point_peano;
          break;
        default:
          dserror(
              "Integration rule only until degree 16 for TET elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::FE::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::tet10>
  {
    static CORE::FE::GaussRule3D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::FE::GaussRule3D::tet_1point;
          break;
        case 2:
          return CORE::FE::GaussRule3D::tet_4point;
          break;
        case 3:
          return CORE::FE::GaussRule3D::tet_5point;
          break;
        case 4:
          return CORE::FE::GaussRule3D::tet_11point;
          break;
        case 5:
          return CORE::FE::GaussRule3D::tet_15point;
          break;
        case 6:
          return CORE::FE::GaussRule3D::tet_24point;
          break;
        case 7:
        case 8:
          return CORE::FE::GaussRule3D::tet_45point;
          break;
        case 9:
          return CORE::FE::GaussRule3D::tet_125point_peano;
          break;
        case 10:
        case 11:
        case 12:
        case 13:
          return CORE::FE::GaussRule3D::tet_343point_peano;
          break;
        case 14:
        case 15:
        case 16:
          return CORE::FE::GaussRule3D::tet_729point_peano;
          break;
        default:
          dserror(
              "Integration rule only until degree 16 for TET elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::FE::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::wedge6>
  {
    static CORE::FE::GaussRule3D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for WEDGE elements defined. Feel free to add the Gauss rule.");
      return CORE::FE::GaussRule3D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::pyramid5>
  {
    static CORE::FE::GaussRule3D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for PYRAMID elements defined. Feel free to add the Gauss rule.");
      return CORE::FE::GaussRule3D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::nurbs8>
  {
    static CORE::FE::GaussRule3D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return CORE::FE::GaussRule3D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::nurbs27>
  {
    static CORE::FE::GaussRule3D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return CORE::FE::GaussRule3D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::quad4>
  {
    static CORE::FE::GaussRule2D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::FE::GaussRule2D::quad_1point;
          break;
        case 2:
        case 3:
          return CORE::FE::GaussRule2D::quad_4point;
          break;
        case 4:
        case 5:
          return CORE::FE::GaussRule2D::quad_9point;
          break;
        case 6:
        case 7:
          return CORE::FE::GaussRule2D::quad_16point;
          break;
        case 8:
        case 9:
          return CORE::FE::GaussRule2D::quad_25point;
          break;
        case 10:
        case 11:
          return CORE::FE::GaussRule2D::quad_36point;
          break;
        case 12:
        case 13:
          return CORE::FE::GaussRule2D::quad_49point;
          break;
        case 14:
        case 15:
          return CORE::FE::GaussRule2D::quad_64point;
          break;
        case 16:
        case 17:
          return CORE::FE::GaussRule2D::quad_81point;
          break;
        case 18:
        case 19:
          return CORE::FE::GaussRule2D::quad_100point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for QUAD elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::FE::GaussRule2D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::quad8>
  {
    static CORE::FE::GaussRule2D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::FE::GaussRule2D::quad_1point;
          break;
        case 2:
        case 3:
          return CORE::FE::GaussRule2D::quad_4point;
          break;
        case 4:
        case 5:
          return CORE::FE::GaussRule2D::quad_9point;
          break;
        case 6:
        case 7:
          return CORE::FE::GaussRule2D::quad_16point;
          break;
        case 8:
        case 9:
          return CORE::FE::GaussRule2D::quad_25point;
          break;
        case 10:
        case 11:
          return CORE::FE::GaussRule2D::quad_36point;
          break;
        case 12:
        case 13:
          return CORE::FE::GaussRule2D::quad_49point;
          break;
        case 14:
        case 15:
          return CORE::FE::GaussRule2D::quad_64point;
          break;
        case 16:
        case 17:
          return CORE::FE::GaussRule2D::quad_81point;
          break;
        case 18:
        case 19:
          return CORE::FE::GaussRule2D::quad_100point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for QUAD elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::FE::GaussRule2D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::quad9>
  {
    static CORE::FE::GaussRule2D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::FE::GaussRule2D::quad_1point;
          break;
        case 2:
        case 3:
          return CORE::FE::GaussRule2D::quad_4point;
          break;
        case 4:
        case 5:
          return CORE::FE::GaussRule2D::quad_9point;
          break;
        case 6:
        case 7:
          return CORE::FE::GaussRule2D::quad_16point;
          break;
        case 8:
        case 9:
          return CORE::FE::GaussRule2D::quad_25point;
          break;
        case 10:
        case 11:
          return CORE::FE::GaussRule2D::quad_36point;
          break;
        case 12:
        case 13:
          return CORE::FE::GaussRule2D::quad_49point;
          break;
        case 14:
        case 15:
          return CORE::FE::GaussRule2D::quad_64point;
          break;
        case 16:
        case 17:
          return CORE::FE::GaussRule2D::quad_81point;
          break;
        case 18:
        case 19:
          return CORE::FE::GaussRule2D::quad_100point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for QUAD elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::FE::GaussRule2D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::tri3>
  {
    static CORE::FE::GaussRule2D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::FE::GaussRule2D::tri_1point;
          break;
        case 2:
          return CORE::FE::GaussRule2D::tri_3point;
          break;
        case 3:
          return CORE::FE::GaussRule2D::tri_4point;
          break;
        case 4:
          return CORE::FE::GaussRule2D::tri_6point;
          break;
        case 5:
          return CORE::FE::GaussRule2D::tri_7point;
          break;
        case 6:
          return CORE::FE::GaussRule2D::tri_12point;
          break;
        case 7:
        case 8:
          return CORE::FE::GaussRule2D::tri_16point;
          break;
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
          return CORE::FE::GaussRule2D::tri_37point;
          break;
        case 14:
        case 15:
          return CORE::FE::GaussRule2D::tri_64point;
          break;
        default:
          dserror(
              "Integration rule only until degree 15 for TRI elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::FE::GaussRule2D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::tri6>
  {
    static CORE::FE::GaussRule2D GetGaussRule(int degree)
    {
      return CORE::FE::GaussRule2D::tri_16point;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::nurbs4>
  {
    static CORE::FE::GaussRule2D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return CORE::FE::GaussRule2D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::nurbs9>
  {
    static CORE::FE::GaussRule2D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return CORE::FE::GaussRule2D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::line2>
  {
    static CORE::FE::GaussRule1D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::FE::GaussRule1D::line_1point;
          break;
        case 2:
        case 3:
          return CORE::FE::GaussRule1D::line_2point;
          break;
        case 4:
        case 5:
          return CORE::FE::GaussRule1D::line_3point;
          break;
        case 6:
        case 7:
          return CORE::FE::GaussRule1D::line_4point;
          break;
        case 8:
        case 9:
          return CORE::FE::GaussRule1D::line_5point;
          break;
        case 10:
        case 11:
          return CORE::FE::GaussRule1D::line_6point;
          break;
        case 12:
        case 13:
          return CORE::FE::GaussRule1D::line_7point;
          break;
        case 14:
        case 15:
          return CORE::FE::GaussRule1D::line_8point;
          break;
        case 16:
        case 17:
          return CORE::FE::GaussRule1D::line_9point;
          break;
        case 18:
        case 19:
          return CORE::FE::GaussRule1D::line_10point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for LINE elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::FE::GaussRule1D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::line3>
  {
    static CORE::FE::GaussRule1D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for LINE elements defined. Feel free to add the Gauss rule.");
      return CORE::FE::GaussRule1D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::nurbs2>
  {
    static CORE::FE::GaussRule1D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return CORE::FE::GaussRule1D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<CORE::FE::CellType::nurbs3>
  {
    static CORE::FE::GaussRule1D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return CORE::FE::GaussRule1D::undefined;
    };
  };


  //! Template Meta Programming version of switch over discretization type
  template <CORE::FE::CellType DISTYPE>
  struct DisTypeToGaussRuleForExactSol
  {
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::hex8>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::hex20>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::hex27>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tet4>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::tet_5point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tet10>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::wedge6>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::pyramid5>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::nurbs8>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::nurbs27>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::quad4>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::quad8>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::quad9>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tri3>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::tri6>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::nurbs4>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_4point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::nurbs9>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::line2>
  {
    static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::line_2point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::line3>
  {
    static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::line_8point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::nurbs2>
  {
    static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::nurbs3>
  {
    static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::undefined;
  };

  //! for each distype provide the m_k needed for stabilization computation
  template <CORE::FE::CellType DISTYPE>
  inline double MK()
  {
    switch (DISTYPE)
    {
      case CORE::FE::CellType::tet4:
      case CORE::FE::CellType::pyramid5:
      case CORE::FE::CellType::hex8:
      case CORE::FE::CellType::wedge6:
      case CORE::FE::CellType::nurbs8:
        return 0.333333333333333333333;
        break;
      case CORE::FE::CellType::hex20:
      case CORE::FE::CellType::hex27:
      case CORE::FE::CellType::tet10:
      case CORE::FE::CellType::wedge15:
      case CORE::FE::CellType::nurbs27:
        return 0.083333333333333333333;
        break;
        // do the 2D case after the more usual 3D case
      case CORE::FE::CellType::tri3:
      case CORE::FE::CellType::quad4:
      case CORE::FE::CellType::nurbs4:
        return 0.333333333333333333333;
        break;
      case CORE::FE::CellType::tri6:
      case CORE::FE::CellType::quad8:
      case CORE::FE::CellType::quad9:
      case CORE::FE::CellType::nurbs9:
        return 0.083333333333333333333;
        break;
        // finally, do the 1D case
      case CORE::FE::CellType::line2:
      case CORE::FE::CellType::nurbs2:
        return 0.333333333333333333333;
        break;
      case CORE::FE::CellType::line3:
      case CORE::FE::CellType::nurbs3:
        return 0.083333333333333333333;
        break;
      default:
        dserror("Element shape not supported.");
        break;
    }
    return -1.0;
  }

  //! Template Meta Programming version of switch over discretization type
  template <CORE::FE::CellType DISTYPE>
  struct DisTypeToStabGaussRule
  {
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::hex8>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::hex20>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::hex27>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::tet4>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::tet_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::tet10>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::tet_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::wedge6>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::wedge_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::pyramid5>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::pyramid_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::nurbs8>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::nurbs27>
  {
    static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::quad4>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::quad8>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::quad9>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::nurbs4>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::nurbs9>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::tri3>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::tri_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::tri6>
  {
    static constexpr CORE::FE::GaussRule2D rule = CORE::FE::GaussRule2D::tri_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::line2>
  {
    static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::line_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::line3>
  {
    static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::line_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::nurbs2>
  {
    static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::line_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<CORE::FE::CellType::nurbs3>
  {
    static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::line_1point;
  };


  //! determine whether there are only two charged species present or not
  bool IsBinaryElectrolyte(const std::vector<double>& valence);

  //! determine indices of the two charged species in case of an binary electrolyte
  std::vector<int> GetIndicesBinaryElectrolyte(const std::vector<double>& valence);

  //! calculate resultant diffusion coefficient for stabilization of binary electrolyte systems
  double CalResDiffCoeff(const std::vector<double>& valence,  ///< valences
      const std::vector<double>& diffus,                      ///< diffusivities
      const std::vector<int>& indices                         ///< indices
  );

  //! identify elements of inflow section
  bool InflowElement(const DRT::Element* ele);

  //! convert implementation type of scalar transport elements into corresponding string for output
  //! purposes
  std::string ImplTypeToString(const INPAR::SCATRA::ImplType impltype);
}  // namespace SCATRA

FOUR_C_NAMESPACE_CLOSE

#endif
