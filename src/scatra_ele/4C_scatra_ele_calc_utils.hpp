/*----------------------------------------------------------------------*/
/*! \file

\brief Utility methods for scatra

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_ELE_CALC_UTILS_HPP
#define FOUR_C_SCATRA_ELE_CALC_UTILS_HPP


#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_integration.hpp"
#include "4C_inpar_scatra.hpp"

FOUR_C_NAMESPACE_OPEN


namespace ScaTra
{
  /*!
  \brief Decide, whether second derivatives are needed  (template version)
   *  In convection-diffusion problems, ONLY N,xx , N,yy and N,zz are needed
   *  to evaluate the laplacian operator for the residual-based stabilization.
   *  Hence, unlike to the Navier-Stokes equations, hex8, wedge6 and pyramid5
   *  return false although they have non-zero MIXED second derivatives.*/
  template <Core::FE::CellType DISTYPE>
  struct Use2ndDerivs
  {
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::hex8>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::tet4>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::wedge6>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::pyramid5>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::nurbs8>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::quad4>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::nurbs4>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::tri3>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::line2>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::nurbs2>
  {
    static constexpr bool use = false;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::hex20>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::hex27>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::nurbs27>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::tet10>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::quad8>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::quad9>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::nurbs9>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::tri6>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::line3>
  {
    static constexpr bool use = true;
  };
  template <>
  struct Use2ndDerivs<Core::FE::CellType::nurbs3>
  {
    static constexpr bool use = true;
  };

  //! Template Meta Programming version of switch over discretization type
  template <Core::FE::CellType DISTYPE>
  struct DisTypeToOptGaussRule
  {
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::hex8>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::hex20>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::hex27>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::tet4>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::tet_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::tet10>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::tet_11point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::wedge6>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::wedge_6point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::pyramid5>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::pyramid_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::nurbs8>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::nurbs27>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::quad4>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::quad8>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::quad9>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::tri3>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::tri_3point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::tri6>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::tri_6point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::nurbs4>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::nurbs9>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::line2>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_2point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::line3>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_3point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::nurbs2>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_2point;
  };
  template <>
  struct DisTypeToOptGaussRule<Core::FE::CellType::nurbs3>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_3point;
  };

  //! Template Meta Programming version of switch over discretization type
  template <Core::FE::CellType DISTYPE>
  struct DisTypeToMatGaussRule
  {
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::hex8>
  {
    static Core::FE::GaussRule3D get_gauss_rule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return Core::FE::GaussRule3D::hex_1point;
          break;
        case 2:
        case 3:
          return Core::FE::GaussRule3D::hex_8point;
          break;
        case 4:
        case 5:
          return Core::FE::GaussRule3D::hex_27point;
          break;
        case 6:
        case 7:
          return Core::FE::GaussRule3D::hex_64point;
          break;
        case 8:
        case 9:
          return Core::FE::GaussRule3D::hex_125point;
          break;
        case 10:
        case 11:
          return Core::FE::GaussRule3D::hex_216point;
          break;
        case 12:
        case 13:
          return Core::FE::GaussRule3D::hex_343point;
          break;
        case 14:
        case 15:
          return Core::FE::GaussRule3D::hex_512point;
          break;
        case 16:
        case 17:
          return Core::FE::GaussRule3D::hex_729point;
          break;
        case 18:
        case 19:
          return Core::FE::GaussRule3D::hex_1000point;
          break;
        default:
          FOUR_C_THROW(
              "Integration rule only until degree 19 for HEX elements defined. You used a degree "
              "of %d",
              degree);
          return Core::FE::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::hex20>
  {
    static Core::FE::GaussRule3D get_gauss_rule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return Core::FE::GaussRule3D::hex_1point;
          break;
        case 2:
        case 3:
          return Core::FE::GaussRule3D::hex_8point;
          break;
        case 4:
        case 5:
          return Core::FE::GaussRule3D::hex_27point;
          break;
        case 6:
        case 7:
          return Core::FE::GaussRule3D::hex_64point;
          break;
        case 8:
        case 9:
          return Core::FE::GaussRule3D::hex_125point;
          break;
        case 10:
        case 11:
          return Core::FE::GaussRule3D::hex_216point;
          break;
        case 12:
        case 13:
          return Core::FE::GaussRule3D::hex_343point;
          break;
        case 14:
        case 15:
          return Core::FE::GaussRule3D::hex_512point;
          break;
        case 16:
        case 17:
          return Core::FE::GaussRule3D::hex_729point;
          break;
        case 18:
        case 19:
          return Core::FE::GaussRule3D::hex_1000point;
          break;
        default:
          FOUR_C_THROW(
              "Integration rule only until degree 19 for HEX elements defined. You used a degree "
              "of %d",
              degree);
          return Core::FE::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::hex27>
  {
    static Core::FE::GaussRule3D get_gauss_rule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return Core::FE::GaussRule3D::hex_1point;
          break;
        case 2:
        case 3:
          return Core::FE::GaussRule3D::hex_8point;
          break;
        case 4:
        case 5:
          return Core::FE::GaussRule3D::hex_27point;
          break;
        case 6:
        case 7:
          return Core::FE::GaussRule3D::hex_64point;
          break;
        case 8:
        case 9:
          return Core::FE::GaussRule3D::hex_125point;
          break;
        case 10:
        case 11:
          return Core::FE::GaussRule3D::hex_216point;
          break;
        case 12:
        case 13:
          return Core::FE::GaussRule3D::hex_343point;
          break;
        case 14:
        case 15:
          return Core::FE::GaussRule3D::hex_512point;
          break;
        case 16:
        case 17:
          return Core::FE::GaussRule3D::hex_729point;
          break;
        case 18:
        case 19:
          return Core::FE::GaussRule3D::hex_1000point;
          break;
        default:
          FOUR_C_THROW(
              "Integration rule only until degree 19 for HEX elements defined. You used a degree "
              "of %d",
              degree);
          return Core::FE::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::tet4>
  {
    static Core::FE::GaussRule3D get_gauss_rule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return Core::FE::GaussRule3D::tet_1point;
          break;
        case 2:
          return Core::FE::GaussRule3D::tet_4point;
          break;
        case 3:
          return Core::FE::GaussRule3D::tet_5point;
          break;
        case 4:
          return Core::FE::GaussRule3D::tet_11point;
          break;
        case 5:
          return Core::FE::GaussRule3D::tet_15point;
          break;
        case 6:
          return Core::FE::GaussRule3D::tet_24point;
          break;
        case 7:
        case 8:
          return Core::FE::GaussRule3D::tet_45point;
          break;
        case 9:
          return Core::FE::GaussRule3D::tet_125point_peano;
          break;
        case 10:
        case 11:
        case 12:
        case 13:
          return Core::FE::GaussRule3D::tet_343point_peano;
          break;
        case 14:
        case 15:
        case 16:
          return Core::FE::GaussRule3D::tet_729point_peano;
          break;
        default:
          FOUR_C_THROW(
              "Integration rule only until degree 16 for TET elements defined. You used a degree "
              "of %d",
              degree);
          return Core::FE::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::tet10>
  {
    static Core::FE::GaussRule3D get_gauss_rule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return Core::FE::GaussRule3D::tet_1point;
          break;
        case 2:
          return Core::FE::GaussRule3D::tet_4point;
          break;
        case 3:
          return Core::FE::GaussRule3D::tet_5point;
          break;
        case 4:
          return Core::FE::GaussRule3D::tet_11point;
          break;
        case 5:
          return Core::FE::GaussRule3D::tet_15point;
          break;
        case 6:
          return Core::FE::GaussRule3D::tet_24point;
          break;
        case 7:
        case 8:
          return Core::FE::GaussRule3D::tet_45point;
          break;
        case 9:
          return Core::FE::GaussRule3D::tet_125point_peano;
          break;
        case 10:
        case 11:
        case 12:
        case 13:
          return Core::FE::GaussRule3D::tet_343point_peano;
          break;
        case 14:
        case 15:
        case 16:
          return Core::FE::GaussRule3D::tet_729point_peano;
          break;
        default:
          FOUR_C_THROW(
              "Integration rule only until degree 16 for TET elements defined. You used a degree "
              "of %d",
              degree);
          return Core::FE::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::wedge6>
  {
    static Core::FE::GaussRule3D get_gauss_rule(int degree)
    {
      FOUR_C_THROW("Gauss rule not for WEDGE elements defined. Feel free to add the Gauss rule.");
      return Core::FE::GaussRule3D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::pyramid5>
  {
    static Core::FE::GaussRule3D get_gauss_rule(int degree)
    {
      FOUR_C_THROW("Gauss rule not for PYRAMID elements defined. Feel free to add the Gauss rule.");
      return Core::FE::GaussRule3D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::nurbs8>
  {
    static Core::FE::GaussRule3D get_gauss_rule(int degree)
    {
      FOUR_C_THROW("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return Core::FE::GaussRule3D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::nurbs27>
  {
    static Core::FE::GaussRule3D get_gauss_rule(int degree)
    {
      FOUR_C_THROW("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return Core::FE::GaussRule3D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::quad4>
  {
    static Core::FE::GaussRule2D get_gauss_rule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return Core::FE::GaussRule2D::quad_1point;
          break;
        case 2:
        case 3:
          return Core::FE::GaussRule2D::quad_4point;
          break;
        case 4:
        case 5:
          return Core::FE::GaussRule2D::quad_9point;
          break;
        case 6:
        case 7:
          return Core::FE::GaussRule2D::quad_16point;
          break;
        case 8:
        case 9:
          return Core::FE::GaussRule2D::quad_25point;
          break;
        case 10:
        case 11:
          return Core::FE::GaussRule2D::quad_36point;
          break;
        case 12:
        case 13:
          return Core::FE::GaussRule2D::quad_49point;
          break;
        case 14:
        case 15:
          return Core::FE::GaussRule2D::quad_64point;
          break;
        case 16:
        case 17:
          return Core::FE::GaussRule2D::quad_81point;
          break;
        case 18:
        case 19:
          return Core::FE::GaussRule2D::quad_100point;
          break;
        default:
          FOUR_C_THROW(
              "Integration rule only until degree 19 for QUAD elements defined. You used a degree "
              "of %d",
              degree);
          return Core::FE::GaussRule2D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::quad8>
  {
    static Core::FE::GaussRule2D get_gauss_rule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return Core::FE::GaussRule2D::quad_1point;
          break;
        case 2:
        case 3:
          return Core::FE::GaussRule2D::quad_4point;
          break;
        case 4:
        case 5:
          return Core::FE::GaussRule2D::quad_9point;
          break;
        case 6:
        case 7:
          return Core::FE::GaussRule2D::quad_16point;
          break;
        case 8:
        case 9:
          return Core::FE::GaussRule2D::quad_25point;
          break;
        case 10:
        case 11:
          return Core::FE::GaussRule2D::quad_36point;
          break;
        case 12:
        case 13:
          return Core::FE::GaussRule2D::quad_49point;
          break;
        case 14:
        case 15:
          return Core::FE::GaussRule2D::quad_64point;
          break;
        case 16:
        case 17:
          return Core::FE::GaussRule2D::quad_81point;
          break;
        case 18:
        case 19:
          return Core::FE::GaussRule2D::quad_100point;
          break;
        default:
          FOUR_C_THROW(
              "Integration rule only until degree 19 for QUAD elements defined. You used a degree "
              "of %d",
              degree);
          return Core::FE::GaussRule2D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::quad9>
  {
    static Core::FE::GaussRule2D get_gauss_rule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return Core::FE::GaussRule2D::quad_1point;
          break;
        case 2:
        case 3:
          return Core::FE::GaussRule2D::quad_4point;
          break;
        case 4:
        case 5:
          return Core::FE::GaussRule2D::quad_9point;
          break;
        case 6:
        case 7:
          return Core::FE::GaussRule2D::quad_16point;
          break;
        case 8:
        case 9:
          return Core::FE::GaussRule2D::quad_25point;
          break;
        case 10:
        case 11:
          return Core::FE::GaussRule2D::quad_36point;
          break;
        case 12:
        case 13:
          return Core::FE::GaussRule2D::quad_49point;
          break;
        case 14:
        case 15:
          return Core::FE::GaussRule2D::quad_64point;
          break;
        case 16:
        case 17:
          return Core::FE::GaussRule2D::quad_81point;
          break;
        case 18:
        case 19:
          return Core::FE::GaussRule2D::quad_100point;
          break;
        default:
          FOUR_C_THROW(
              "Integration rule only until degree 19 for QUAD elements defined. You used a degree "
              "of %d",
              degree);
          return Core::FE::GaussRule2D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::tri3>
  {
    static Core::FE::GaussRule2D get_gauss_rule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return Core::FE::GaussRule2D::tri_1point;
          break;
        case 2:
          return Core::FE::GaussRule2D::tri_3point;
          break;
        case 3:
          return Core::FE::GaussRule2D::tri_4point;
          break;
        case 4:
          return Core::FE::GaussRule2D::tri_6point;
          break;
        case 5:
          return Core::FE::GaussRule2D::tri_7point;
          break;
        case 6:
          return Core::FE::GaussRule2D::tri_12point;
          break;
        case 7:
        case 8:
          return Core::FE::GaussRule2D::tri_16point;
          break;
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
          return Core::FE::GaussRule2D::tri_37point;
          break;
        case 14:
        case 15:
          return Core::FE::GaussRule2D::tri_64point;
          break;
        default:
          FOUR_C_THROW(
              "Integration rule only until degree 15 for TRI elements defined. You used a degree "
              "of %d",
              degree);
          return Core::FE::GaussRule2D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::tri6>
  {
    static Core::FE::GaussRule2D get_gauss_rule(int degree)
    {
      return Core::FE::GaussRule2D::tri_16point;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::nurbs4>
  {
    static Core::FE::GaussRule2D get_gauss_rule(int degree)
    {
      FOUR_C_THROW("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return Core::FE::GaussRule2D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::nurbs9>
  {
    static Core::FE::GaussRule2D get_gauss_rule(int degree)
    {
      FOUR_C_THROW("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return Core::FE::GaussRule2D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::line2>
  {
    static Core::FE::GaussRule1D get_gauss_rule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return Core::FE::GaussRule1D::line_1point;
          break;
        case 2:
        case 3:
          return Core::FE::GaussRule1D::line_2point;
          break;
        case 4:
        case 5:
          return Core::FE::GaussRule1D::line_3point;
          break;
        case 6:
        case 7:
          return Core::FE::GaussRule1D::line_4point;
          break;
        case 8:
        case 9:
          return Core::FE::GaussRule1D::line_5point;
          break;
        case 10:
        case 11:
          return Core::FE::GaussRule1D::line_6point;
          break;
        case 12:
        case 13:
          return Core::FE::GaussRule1D::line_7point;
          break;
        case 14:
        case 15:
          return Core::FE::GaussRule1D::line_8point;
          break;
        case 16:
        case 17:
          return Core::FE::GaussRule1D::line_9point;
          break;
        case 18:
        case 19:
          return Core::FE::GaussRule1D::line_10point;
          break;
        default:
          FOUR_C_THROW(
              "Integration rule only until degree 19 for LINE elements defined. You used a degree "
              "of %d",
              degree);
          return Core::FE::GaussRule1D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::line3>
  {
    static Core::FE::GaussRule1D get_gauss_rule(int degree)
    {
      FOUR_C_THROW("Gauss rule not for LINE elements defined. Feel free to add the Gauss rule.");
      return Core::FE::GaussRule1D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::nurbs2>
  {
    static Core::FE::GaussRule1D get_gauss_rule(int degree)
    {
      FOUR_C_THROW("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return Core::FE::GaussRule1D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<Core::FE::CellType::nurbs3>
  {
    static Core::FE::GaussRule1D get_gauss_rule(int degree)
    {
      FOUR_C_THROW("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return Core::FE::GaussRule1D::undefined;
    };
  };


  //! Template Meta Programming version of switch over discretization type
  template <Core::FE::CellType DISTYPE>
  struct DisTypeToGaussRuleForExactSol
  {
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::hex8>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::hex20>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::hex27>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::tet4>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::tet_5point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::tet10>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::wedge6>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::pyramid5>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::nurbs8>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::nurbs27>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::quad4>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::quad8>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::quad9>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::tri3>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::tri6>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::nurbs4>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_4point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::nurbs9>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::line2>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_2point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::line3>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_8point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::nurbs2>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::nurbs3>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::undefined;
  };

  //! for each distype provide the m_k needed for stabilization computation
  template <Core::FE::CellType DISTYPE>
  inline double MK()
  {
    switch (DISTYPE)
    {
      case Core::FE::CellType::tet4:
      case Core::FE::CellType::pyramid5:
      case Core::FE::CellType::hex8:
      case Core::FE::CellType::wedge6:
      case Core::FE::CellType::nurbs8:
        return 0.333333333333333333333;
        break;
      case Core::FE::CellType::hex20:
      case Core::FE::CellType::hex27:
      case Core::FE::CellType::tet10:
      case Core::FE::CellType::wedge15:
      case Core::FE::CellType::nurbs27:
        return 0.083333333333333333333;
        break;
        // do the 2D case after the more usual 3D case
      case Core::FE::CellType::tri3:
      case Core::FE::CellType::quad4:
      case Core::FE::CellType::nurbs4:
        return 0.333333333333333333333;
        break;
      case Core::FE::CellType::tri6:
      case Core::FE::CellType::quad8:
      case Core::FE::CellType::quad9:
      case Core::FE::CellType::nurbs9:
        return 0.083333333333333333333;
        break;
        // finally, do the 1D case
      case Core::FE::CellType::line2:
      case Core::FE::CellType::nurbs2:
        return 0.333333333333333333333;
        break;
      case Core::FE::CellType::line3:
      case Core::FE::CellType::nurbs3:
        return 0.083333333333333333333;
        break;
      default:
        FOUR_C_THROW("Element shape not supported.");
        break;
    }
    return -1.0;
  }

  //! Template Meta Programming version of switch over discretization type
  template <Core::FE::CellType DISTYPE>
  struct DisTypeToStabGaussRule
  {
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::hex8>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::hex20>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::hex27>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::tet4>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::tet_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::tet10>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::tet_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::wedge6>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::wedge_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::pyramid5>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::pyramid_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::nurbs8>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::nurbs27>
  {
    static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::quad4>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::quad8>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::quad9>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::nurbs4>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::nurbs9>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::tri3>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::tri_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::tri6>
  {
    static constexpr Core::FE::GaussRule2D rule = Core::FE::GaussRule2D::tri_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::line2>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::line3>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::nurbs2>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<Core::FE::CellType::nurbs3>
  {
    static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_1point;
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
  bool inflow_element(const Core::Elements::Element* ele);

  //! convert implementation type of scalar transport elements into corresponding string for output
  //! purposes
  std::string ImplTypeToString(const Inpar::ScaTra::ImplType impltype);
}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
