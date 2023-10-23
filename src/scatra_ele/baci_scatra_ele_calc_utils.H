/*----------------------------------------------------------------------*/
/*! \file

\brief Utility methods for scatra

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_ELE_CALC_UTILS_H
#define BACI_SCATRA_ELE_CALC_UTILS_H


#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_inpar_scatra.H"
#include "baci_lib_element.H"


namespace SCATRA
{
  /*!
  \brief Decide, whether second derivatives are needed  (template version)
   *  In convection-diffusion problems, ONLY N,xx , N,yy and N,zz are needed
   *  to evaluate the laplacian operator for the residual-based stabilization.
   *  Hence, unlike to the Navier-Stokes equations, hex8, wedge6 and pyramid5
   *  return false although they have non-zero MIXED second derivatives.*/
  template <DRT::Element::DiscretizationType DISTYPE>
  struct Use2ndDerivs
  {
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::hex8>
  {
    static const bool use = false;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::tet4>
  {
    static const bool use = false;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::wedge6>
  {
    static const bool use = false;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::pyramid5>
  {
    static const bool use = false;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::nurbs8>
  {
    static const bool use = false;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::quad4>
  {
    static const bool use = false;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::nurbs4>
  {
    static const bool use = false;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::tri3>
  {
    static const bool use = false;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::line2>
  {
    static const bool use = false;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::nurbs2>
  {
    static const bool use = false;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::hex20>
  {
    static const bool use = true;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::hex27>
  {
    static const bool use = true;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::nurbs27>
  {
    static const bool use = true;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::tet10>
  {
    static const bool use = true;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::quad8>
  {
    static const bool use = true;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::quad9>
  {
    static const bool use = true;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::nurbs9>
  {
    static const bool use = true;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::tri6>
  {
    static const bool use = true;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::line3>
  {
    static const bool use = true;
  };
  template <>
  struct Use2ndDerivs<DRT::Element::DiscretizationType::nurbs3>
  {
    static const bool use = true;
  };

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
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::tet4>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::tet_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::tet10>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::tet_11point;
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
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::nurbs8>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_8point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::nurbs27>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
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
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::nurbs4>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_4point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::nurbs9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
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
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::nurbs2>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_2point;
  };
  template <>
  struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::nurbs3>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_3point;
  };

  //! Template Meta Programming version of switch over discretization type
  template <DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToMatGaussRule
  {
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::hex8>
  {
    static CORE::DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::DRT::UTILS::GaussRule3D::hex_1point;
          break;
        case 2:
        case 3:
          return CORE::DRT::UTILS::GaussRule3D::hex_8point;
          break;
        case 4:
        case 5:
          return CORE::DRT::UTILS::GaussRule3D::hex_27point;
          break;
        case 6:
        case 7:
          return CORE::DRT::UTILS::GaussRule3D::hex_64point;
          break;
        case 8:
        case 9:
          return CORE::DRT::UTILS::GaussRule3D::hex_125point;
          break;
        case 10:
        case 11:
          return CORE::DRT::UTILS::GaussRule3D::hex_216point;
          break;
        case 12:
        case 13:
          return CORE::DRT::UTILS::GaussRule3D::hex_343point;
          break;
        case 14:
        case 15:
          return CORE::DRT::UTILS::GaussRule3D::hex_512point;
          break;
        case 16:
        case 17:
          return CORE::DRT::UTILS::GaussRule3D::hex_729point;
          break;
        case 18:
        case 19:
          return CORE::DRT::UTILS::GaussRule3D::hex_1000point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for HEX elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::DRT::UTILS::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::hex20>
  {
    static CORE::DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::DRT::UTILS::GaussRule3D::hex_1point;
          break;
        case 2:
        case 3:
          return CORE::DRT::UTILS::GaussRule3D::hex_8point;
          break;
        case 4:
        case 5:
          return CORE::DRT::UTILS::GaussRule3D::hex_27point;
          break;
        case 6:
        case 7:
          return CORE::DRT::UTILS::GaussRule3D::hex_64point;
          break;
        case 8:
        case 9:
          return CORE::DRT::UTILS::GaussRule3D::hex_125point;
          break;
        case 10:
        case 11:
          return CORE::DRT::UTILS::GaussRule3D::hex_216point;
          break;
        case 12:
        case 13:
          return CORE::DRT::UTILS::GaussRule3D::hex_343point;
          break;
        case 14:
        case 15:
          return CORE::DRT::UTILS::GaussRule3D::hex_512point;
          break;
        case 16:
        case 17:
          return CORE::DRT::UTILS::GaussRule3D::hex_729point;
          break;
        case 18:
        case 19:
          return CORE::DRT::UTILS::GaussRule3D::hex_1000point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for HEX elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::DRT::UTILS::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::hex27>
  {
    static CORE::DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::DRT::UTILS::GaussRule3D::hex_1point;
          break;
        case 2:
        case 3:
          return CORE::DRT::UTILS::GaussRule3D::hex_8point;
          break;
        case 4:
        case 5:
          return CORE::DRT::UTILS::GaussRule3D::hex_27point;
          break;
        case 6:
        case 7:
          return CORE::DRT::UTILS::GaussRule3D::hex_64point;
          break;
        case 8:
        case 9:
          return CORE::DRT::UTILS::GaussRule3D::hex_125point;
          break;
        case 10:
        case 11:
          return CORE::DRT::UTILS::GaussRule3D::hex_216point;
          break;
        case 12:
        case 13:
          return CORE::DRT::UTILS::GaussRule3D::hex_343point;
          break;
        case 14:
        case 15:
          return CORE::DRT::UTILS::GaussRule3D::hex_512point;
          break;
        case 16:
        case 17:
          return CORE::DRT::UTILS::GaussRule3D::hex_729point;
          break;
        case 18:
        case 19:
          return CORE::DRT::UTILS::GaussRule3D::hex_1000point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for HEX elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::DRT::UTILS::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::tet4>
  {
    static CORE::DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::DRT::UTILS::GaussRule3D::tet_1point;
          break;
        case 2:
          return CORE::DRT::UTILS::GaussRule3D::tet_4point;
          break;
        case 3:
          return CORE::DRT::UTILS::GaussRule3D::tet_5point;
          break;
        case 4:
          return CORE::DRT::UTILS::GaussRule3D::tet_11point;
          break;
        case 5:
          return CORE::DRT::UTILS::GaussRule3D::tet_15point;
          break;
        case 6:
          return CORE::DRT::UTILS::GaussRule3D::tet_24point;
          break;
        case 7:
        case 8:
          return CORE::DRT::UTILS::GaussRule3D::tet_45point;
          break;
        case 9:
          return CORE::DRT::UTILS::GaussRule3D::tet_125point_peano;
          break;
        case 10:
        case 11:
        case 12:
        case 13:
          return CORE::DRT::UTILS::GaussRule3D::tet_343point_peano;
          break;
        case 14:
        case 15:
        case 16:
          return CORE::DRT::UTILS::GaussRule3D::tet_729point_peano;
          break;
        default:
          dserror(
              "Integration rule only until degree 16 for TET elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::DRT::UTILS::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::tet10>
  {
    static CORE::DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::DRT::UTILS::GaussRule3D::tet_1point;
          break;
        case 2:
          return CORE::DRT::UTILS::GaussRule3D::tet_4point;
          break;
        case 3:
          return CORE::DRT::UTILS::GaussRule3D::tet_5point;
          break;
        case 4:
          return CORE::DRT::UTILS::GaussRule3D::tet_11point;
          break;
        case 5:
          return CORE::DRT::UTILS::GaussRule3D::tet_15point;
          break;
        case 6:
          return CORE::DRT::UTILS::GaussRule3D::tet_24point;
          break;
        case 7:
        case 8:
          return CORE::DRT::UTILS::GaussRule3D::tet_45point;
          break;
        case 9:
          return CORE::DRT::UTILS::GaussRule3D::tet_125point_peano;
          break;
        case 10:
        case 11:
        case 12:
        case 13:
          return CORE::DRT::UTILS::GaussRule3D::tet_343point_peano;
          break;
        case 14:
        case 15:
        case 16:
          return CORE::DRT::UTILS::GaussRule3D::tet_729point_peano;
          break;
        default:
          dserror(
              "Integration rule only until degree 16 for TET elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::DRT::UTILS::GaussRule3D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::wedge6>
  {
    static CORE::DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for WEDGE elements defined. Feel free to add the Gauss rule.");
      return CORE::DRT::UTILS::GaussRule3D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::pyramid5>
  {
    static CORE::DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for PYRAMID elements defined. Feel free to add the Gauss rule.");
      return CORE::DRT::UTILS::GaussRule3D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::nurbs8>
  {
    static CORE::DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return CORE::DRT::UTILS::GaussRule3D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::nurbs27>
  {
    static CORE::DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return CORE::DRT::UTILS::GaussRule3D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::quad4>
  {
    static CORE::DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::DRT::UTILS::GaussRule2D::quad_1point;
          break;
        case 2:
        case 3:
          return CORE::DRT::UTILS::GaussRule2D::quad_4point;
          break;
        case 4:
        case 5:
          return CORE::DRT::UTILS::GaussRule2D::quad_9point;
          break;
        case 6:
        case 7:
          return CORE::DRT::UTILS::GaussRule2D::quad_16point;
          break;
        case 8:
        case 9:
          return CORE::DRT::UTILS::GaussRule2D::quad_25point;
          break;
        case 10:
        case 11:
          return CORE::DRT::UTILS::GaussRule2D::quad_36point;
          break;
        case 12:
        case 13:
          return CORE::DRT::UTILS::GaussRule2D::quad_49point;
          break;
        case 14:
        case 15:
          return CORE::DRT::UTILS::GaussRule2D::quad_64point;
          break;
        case 16:
        case 17:
          return CORE::DRT::UTILS::GaussRule2D::quad_81point;
          break;
        case 18:
        case 19:
          return CORE::DRT::UTILS::GaussRule2D::quad_100point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for QUAD elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::DRT::UTILS::GaussRule2D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::quad8>
  {
    static CORE::DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::DRT::UTILS::GaussRule2D::quad_1point;
          break;
        case 2:
        case 3:
          return CORE::DRT::UTILS::GaussRule2D::quad_4point;
          break;
        case 4:
        case 5:
          return CORE::DRT::UTILS::GaussRule2D::quad_9point;
          break;
        case 6:
        case 7:
          return CORE::DRT::UTILS::GaussRule2D::quad_16point;
          break;
        case 8:
        case 9:
          return CORE::DRT::UTILS::GaussRule2D::quad_25point;
          break;
        case 10:
        case 11:
          return CORE::DRT::UTILS::GaussRule2D::quad_36point;
          break;
        case 12:
        case 13:
          return CORE::DRT::UTILS::GaussRule2D::quad_49point;
          break;
        case 14:
        case 15:
          return CORE::DRT::UTILS::GaussRule2D::quad_64point;
          break;
        case 16:
        case 17:
          return CORE::DRT::UTILS::GaussRule2D::quad_81point;
          break;
        case 18:
        case 19:
          return CORE::DRT::UTILS::GaussRule2D::quad_100point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for QUAD elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::DRT::UTILS::GaussRule2D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::quad9>
  {
    static CORE::DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::DRT::UTILS::GaussRule2D::quad_1point;
          break;
        case 2:
        case 3:
          return CORE::DRT::UTILS::GaussRule2D::quad_4point;
          break;
        case 4:
        case 5:
          return CORE::DRT::UTILS::GaussRule2D::quad_9point;
          break;
        case 6:
        case 7:
          return CORE::DRT::UTILS::GaussRule2D::quad_16point;
          break;
        case 8:
        case 9:
          return CORE::DRT::UTILS::GaussRule2D::quad_25point;
          break;
        case 10:
        case 11:
          return CORE::DRT::UTILS::GaussRule2D::quad_36point;
          break;
        case 12:
        case 13:
          return CORE::DRT::UTILS::GaussRule2D::quad_49point;
          break;
        case 14:
        case 15:
          return CORE::DRT::UTILS::GaussRule2D::quad_64point;
          break;
        case 16:
        case 17:
          return CORE::DRT::UTILS::GaussRule2D::quad_81point;
          break;
        case 18:
        case 19:
          return CORE::DRT::UTILS::GaussRule2D::quad_100point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for QUAD elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::DRT::UTILS::GaussRule2D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::tri3>
  {
    static CORE::DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::DRT::UTILS::GaussRule2D::tri_1point;
          break;
        case 2:
          return CORE::DRT::UTILS::GaussRule2D::tri_3point;
          break;
        case 3:
          return CORE::DRT::UTILS::GaussRule2D::tri_4point;
          break;
        case 4:
          return CORE::DRT::UTILS::GaussRule2D::tri_6point;
          break;
        case 5:
          return CORE::DRT::UTILS::GaussRule2D::tri_7point;
          break;
        case 6:
          return CORE::DRT::UTILS::GaussRule2D::tri_12point;
          break;
        case 7:
        case 8:
          return CORE::DRT::UTILS::GaussRule2D::tri_16point;
          break;
        case 9:
        case 10:
        case 11:
        case 12:
        case 13:
          return CORE::DRT::UTILS::GaussRule2D::tri_37point;
          break;
        case 14:
        case 15:
          return CORE::DRT::UTILS::GaussRule2D::tri_64point;
          break;
        default:
          dserror(
              "Integration rule only until degree 15 for TRI elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::DRT::UTILS::GaussRule2D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::tri6>
  {
    static CORE::DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      return CORE::DRT::UTILS::GaussRule2D::tri_16point;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::nurbs4>
  {
    static CORE::DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return CORE::DRT::UTILS::GaussRule2D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::nurbs9>
  {
    static CORE::DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return CORE::DRT::UTILS::GaussRule2D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::line2>
  {
    static CORE::DRT::UTILS::GaussRule1D GetGaussRule(int degree)
    {
      switch (degree)
      {
        case 0:
        case 1:
          return CORE::DRT::UTILS::GaussRule1D::line_1point;
          break;
        case 2:
        case 3:
          return CORE::DRT::UTILS::GaussRule1D::line_2point;
          break;
        case 4:
        case 5:
          return CORE::DRT::UTILS::GaussRule1D::line_3point;
          break;
        case 6:
        case 7:
          return CORE::DRT::UTILS::GaussRule1D::line_4point;
          break;
        case 8:
        case 9:
          return CORE::DRT::UTILS::GaussRule1D::line_5point;
          break;
        case 10:
        case 11:
          return CORE::DRT::UTILS::GaussRule1D::line_6point;
          break;
        case 12:
        case 13:
          return CORE::DRT::UTILS::GaussRule1D::line_7point;
          break;
        case 14:
        case 15:
          return CORE::DRT::UTILS::GaussRule1D::line_8point;
          break;
        case 16:
        case 17:
          return CORE::DRT::UTILS::GaussRule1D::line_9point;
          break;
        case 18:
        case 19:
          return CORE::DRT::UTILS::GaussRule1D::line_10point;
          break;
        default:
          dserror(
              "Integration rule only until degree 19 for LINE elements defined. You used a degree "
              "of %d",
              degree);
          return CORE::DRT::UTILS::GaussRule1D::undefined;
          break;
      }
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::line3>
  {
    static CORE::DRT::UTILS::GaussRule1D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for LINE elements defined. Feel free to add the Gauss rule.");
      return CORE::DRT::UTILS::GaussRule1D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::nurbs2>
  {
    static CORE::DRT::UTILS::GaussRule1D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return CORE::DRT::UTILS::GaussRule1D::undefined;
    };
  };
  template <>
  struct DisTypeToMatGaussRule<DRT::Element::DiscretizationType::nurbs3>
  {
    static CORE::DRT::UTILS::GaussRule1D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return CORE::DRT::UTILS::GaussRule1D::undefined;
    };
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
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::tet4>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::tet_5point;
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
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::nurbs8>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::nurbs27>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_27point;
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
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::nurbs4>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_4point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::nurbs9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_9point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::line2>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_2point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::line3>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_8point;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::nurbs2>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
  };
  template <>
  struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::nurbs3>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
  };

  //! for each distype provide the m_k needed for stabilization computation
  template <DRT::Element::DiscretizationType DISTYPE>
  inline double MK()
  {
    switch (DISTYPE)
    {
      case DRT::Element::DiscretizationType::tet4:
      case DRT::Element::DiscretizationType::pyramid5:
      case DRT::Element::DiscretizationType::hex8:
      case DRT::Element::DiscretizationType::wedge6:
      case DRT::Element::DiscretizationType::nurbs8:
        return 0.333333333333333333333;
        break;
      case DRT::Element::DiscretizationType::hex20:
      case DRT::Element::DiscretizationType::hex27:
      case DRT::Element::DiscretizationType::tet10:
      case DRT::Element::DiscretizationType::wedge15:
      case DRT::Element::DiscretizationType::nurbs27:
        return 0.083333333333333333333;
        break;
        // do the 2D case after the more usual 3D case
      case DRT::Element::DiscretizationType::tri3:
      case DRT::Element::DiscretizationType::quad4:
      case DRT::Element::DiscretizationType::nurbs4:
        return 0.333333333333333333333;
        break;
      case DRT::Element::DiscretizationType::tri6:
      case DRT::Element::DiscretizationType::quad8:
      case DRT::Element::DiscretizationType::quad9:
      case DRT::Element::DiscretizationType::nurbs9:
        return 0.083333333333333333333;
        break;
        // finally, do the 1D case
      case DRT::Element::DiscretizationType::line2:
      case DRT::Element::DiscretizationType::nurbs2:
        return 0.333333333333333333333;
        break;
      case DRT::Element::DiscretizationType::line3:
      case DRT::Element::DiscretizationType::nurbs3:
        return 0.083333333333333333333;
        break;
      default:
        dserror("Element shape not supported.");
        break;
    }
    return -1.0;
  }

  //! Template Meta Programming version of switch over discretization type
  template <DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToStabGaussRule
  {
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::hex8>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::hex20>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::hex27>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::tet4>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::tet_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::tet10>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::tet_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::wedge6>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::wedge_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::pyramid5>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::pyramid_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::nurbs8>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::nurbs27>
  {
    static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::quad4>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::quad8>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::quad9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::nurbs4>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::nurbs9>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::quad_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::tri3>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::tri6>
  {
    static const CORE::DRT::UTILS::GaussRule2D rule = CORE::DRT::UTILS::GaussRule2D::tri_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::line2>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::line3>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::nurbs2>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_1point;
  };
  template <>
  struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::nurbs3>
  {
    static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_1point;
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
  const std::string ImplTypeToString(const INPAR::SCATRA::ImplType impltype);
}  // namespace SCATRA

#endif
