/*---------------------------------------------------------------------*/
/*! \file

\brief templated version for selecting integration points from element type

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef BACI_LIB_ELEMENT_INTEGRATION_SELECT_H
#define BACI_LIB_ELEMENT_INTEGRATION_SELECT_H


#include "baci_discretization_fem_general_utils_integration.H"
#include "baci_lib_element.H"


namespace DRT
{
  namespace ELEMENTS
  {
    //! decide, whether second derivatives are needed  (template version)
    /*  Hence, unlike to the Navier-Stokes equations, hex8, wedge6 and pyramid5
     *  return false although they have non-zero MIXED second derivatives.*/

    // tet10 -> integrationrule_tet5  ???
    // wedge15 -> integrationrule_wedge9  ???

    //! Template Meta Programming version of switch over discretization type
    //! get optimal Gauss rule
    template <DRT::Element::DiscretizationType distype>
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
    struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::wedge15>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::wedge_9point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::DiscretizationType::pyramid5>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule =
          CORE::DRT::UTILS::GaussRule3D::pyramid_8point;
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

    //! one-point integration rule should be exact for all elements with straight lines and surfaces
    //! Otherwise, it is assumed as an accurate enough approximation

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
    struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::hex16>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::hex_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::hex18>
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
    struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::wedge15>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule = CORE::DRT::UTILS::GaussRule3D::wedge_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::pyramid5>
    {
      static const CORE::DRT::UTILS::GaussRule3D rule =
          CORE::DRT::UTILS::GaussRule3D::pyramid_1point;
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
    struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::line2>
    {
      static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::DiscretizationType::line3>
    {
      static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::line_1point;
    };

    //! wedge6 and pyramid5 have nonzero 2nd mixed derivatives, but they are neglected
    //! quad4 and hex8 are pseudo higherOrder elements
    //! Template Meta Programming version of switch over discretization type
    template <DRT::Element::DiscretizationType distype>
    struct IsHigherOrder
    {
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::hex8>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::hex20>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::hex27>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::tet4>
    {
      static const bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::tet10>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::wedge6>
    {
      static const bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::wedge15>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::pyramid5>
    {
      static const bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::quad4>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::quad8>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::quad9>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::tri3>
    {
      static const bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::tri6>
    {
      static const bool ishigherorder = true;
    };
    // template<> struct IsHigherOrder<DRT::Element::DiscretizationType::line2>   {static const bool
    // ishigherorder = false; }; template<> struct
    // IsHigherOrder<DRT::Element::DiscretizationType::line3>   {static const bool ishigherorder =
    // false; }; template<> struct IsHigherOrder<DRT::Element::DiscretizationType::nurbs2> {static
    // const bool ishigherorder = true; }; template<> struct
    // IsHigherOrder<DRT::Element::DiscretizationType::nurbs3> {static const bool ishigherorder =
    // true; };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::nurbs4>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::nurbs9>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::nurbs8>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::DiscretizationType::nurbs27>
    {
      static const bool ishigherorder = true;
    };


    //! Template Meta Programming version of switch over discretization type
    template <DRT::Element::DiscretizationType distype>
    struct IsNurbs
    {
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::hex8>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::hex20>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::hex27>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::tet4>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::tet10>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::wedge6>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::wedge15>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::pyramid5>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::quad4>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::quad8>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::quad9>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::tri3>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::tri6>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::line2>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::line3>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::nurbs2>
    {
      static const bool isnurbs = true;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::nurbs3>
    {
      static const bool isnurbs = true;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::nurbs4>
    {
      static const bool isnurbs = true;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::nurbs9>
    {
      static const bool isnurbs = true;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::nurbs8>
    {
      static const bool isnurbs = true;
    };
    template <>
    struct IsNurbs<DRT::Element::DiscretizationType::nurbs27>
    {
      static const bool isnurbs = true;
    };


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
      static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::DiscretizationType::line3>
    {
      static const CORE::DRT::UTILS::GaussRule1D rule = CORE::DRT::UTILS::GaussRule1D::undefined;
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
    //! provide parameter m_k for each distype required for certain stabilization parameter
    //! definitions
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

  }  // namespace ELEMENTS

}  // namespace DRT

#endif
