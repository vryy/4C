/*!
\file drt_element_integration_select.H

\brief templated version for selecting integration points from element type

<pre>
\level 1
\maintainer Martin Kronbichler
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef DRT_ELEMENT_INTEGRATION_SELECT_H
#define DRT_ELEMENT_INTEGRATION_SELECT_H


#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_integration.H"


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
    struct DisTypeToOptGaussRule<DRT::Element::hex8>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_8point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::hex20>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_27point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::hex27>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_27point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::tet4>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_tet_4point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::tet10>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_tet_11point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::wedge6>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_wedge_6point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::wedge15>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_wedge_9point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::pyramid5>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_pyramid_8point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::nurbs8>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_8point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::nurbs27>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_27point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::quad4>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_4point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::quad8>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_9point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::quad9>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_9point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::tri3>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_tri_3point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::tri6>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_tri_6point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::nurbs4>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_4point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::nurbs9>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_9point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::line2>
    {
      static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule_line_2point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::line3>
    {
      static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule_line_3point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::nurbs2>
    {
      static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule_line_2point;
    };
    template <>
    struct DisTypeToOptGaussRule<DRT::Element::nurbs3>
    {
      static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule_line_3point;
    };

    //! one-point integration rule should be exact for all elements with straight lines and surfaces
    //! Otherwise, it is assumed as an accurate enough approximation

    //! Template Meta Programming version of switch over discretization type
    template <DRT::Element::DiscretizationType DISTYPE>
    struct DisTypeToStabGaussRule
    {
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::hex8>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::hex20>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::hex27>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::hex16>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::hex18>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::tet4>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_tet_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::tet10>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_tet_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::wedge6>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_wedge_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::wedge15>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_wedge_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::pyramid5>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_pyramid_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::nurbs8>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::nurbs27>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::quad4>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::quad8>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::quad9>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::tri3>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_tri_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::tri6>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_tri_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::nurbs4>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::nurbs9>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::line2>
    {
      static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule_line_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<DRT::Element::line3>
    {
      static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule_line_1point;
    };

    //! wedge6 and pyramid5 have nonzero 2nd mixed derivatives, but they are neglected
    //! quad4 and hex8 are pseudo higherOrder elements
    //! Template Meta Programming version of switch over discretization type
    template <DRT::Element::DiscretizationType distype>
    struct IsHigherOrder
    {
    };
    template <>
    struct IsHigherOrder<DRT::Element::hex8>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::hex20>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::hex27>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::tet4>
    {
      static const bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<DRT::Element::tet10>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::wedge6>
    {
      static const bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<DRT::Element::wedge15>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::pyramid5>
    {
      static const bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<DRT::Element::quad4>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::quad8>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::quad9>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::tri3>
    {
      static const bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<DRT::Element::tri6>
    {
      static const bool ishigherorder = true;
    };
    // template<> struct IsHigherOrder<DRT::Element::line2>   {static const bool ishigherorder =
    // false; }; template<> struct IsHigherOrder<DRT::Element::line3>   {static const bool
    // ishigherorder = false; }; template<> struct IsHigherOrder<DRT::Element::nurbs2>  {static
    // const bool ishigherorder = true; }; template<> struct IsHigherOrder<DRT::Element::nurbs3>
    // {static const bool ishigherorder = true; };
    template <>
    struct IsHigherOrder<DRT::Element::nurbs4>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::nurbs9>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::nurbs8>
    {
      static const bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<DRT::Element::nurbs27>
    {
      static const bool ishigherorder = true;
    };


    //! Template Meta Programming version of switch over discretization type
    template <DRT::Element::DiscretizationType distype>
    struct IsNurbs
    {
    };
    template <>
    struct IsNurbs<DRT::Element::hex8>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::hex20>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::hex27>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::tet4>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::tet10>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::wedge6>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::wedge15>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::pyramid5>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::quad4>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::quad8>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::quad9>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::tri3>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::tri6>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::line2>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::line3>
    {
      static const bool isnurbs = false;
    };
    template <>
    struct IsNurbs<DRT::Element::nurbs2>
    {
      static const bool isnurbs = true;
    };
    template <>
    struct IsNurbs<DRT::Element::nurbs3>
    {
      static const bool isnurbs = true;
    };
    template <>
    struct IsNurbs<DRT::Element::nurbs4>
    {
      static const bool isnurbs = true;
    };
    template <>
    struct IsNurbs<DRT::Element::nurbs9>
    {
      static const bool isnurbs = true;
    };
    template <>
    struct IsNurbs<DRT::Element::nurbs8>
    {
      static const bool isnurbs = true;
    };
    template <>
    struct IsNurbs<DRT::Element::nurbs27>
    {
      static const bool isnurbs = true;
    };


    template <DRT::Element::DiscretizationType DISTYPE>
    struct DisTypeToGaussRuleForExactSol
    {
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::hex8>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_27point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::hex20>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_27point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::hex27>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_27point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::tet4>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule3D_undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::tet10>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule3D_undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::wedge6>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule3D_undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::pyramid5>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule3D_undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs8>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_27point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs27>
    {
      static const DRT::UTILS::GaussRule3D rule = DRT::UTILS::intrule_hex_27point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::quad4>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_9point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::quad8>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_9point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::quad9>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_9point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::tri3>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule2D_undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::tri6>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule2D_undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs4>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_4point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs9>
    {
      static const DRT::UTILS::GaussRule2D rule = DRT::UTILS::intrule_quad_9point;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::line2>
    {
      static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule1D_undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::line3>
    {
      static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule1D_undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs2>
    {
      static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule1D_undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs3>
    {
      static const DRT::UTILS::GaussRule1D rule = DRT::UTILS::intrule1D_undefined;
    };
    /*
      //! Template Meta Programming version of switch over discretization type
      template<DRT::Element::DiscretizationType DISTYPE>
      struct DisTypeToGaussRuleForExactSol {};
      template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::hex8>    {static const
      DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;}; template<> struct
      DisTypeToGaussRuleForExactSol<DRT::Element::hex20>   {static const DRT::UTILS::GaussRule3D
      rule =  DRT::UTILS::intrule_hex_27point;}; template<> struct
      DisTypeToGaussRuleForExactSol<DRT::Element::hex27>   {static const DRT::UTILS::GaussRule3D
      rule =  DRT::UTILS::intrule_hex_27point;}; template<> struct
      DisTypeToGaussRuleForExactSol<DRT::Element::tet4>    {static const DRT::UTILS::GaussRule3D
      rule =  DRT::UTILS::intrule3D_undefined;}; template<> struct
      DisTypeToGaussRuleForExactSol<DRT::Element::tet10>   {static const DRT::UTILS::GaussRule3D
      rule =  DRT::UTILS::intrule3D_undefined;}; template<> struct
      DisTypeToGaussRuleForExactSol<DRT::Element::wedge6>  {static const DRT::UTILS::GaussRule3D
      rule =  DRT::UTILS::intrule3D_undefined;}; template<> struct
      DisTypeToGaussRuleForExactSol<DRT::Element::pyramid5>{static const DRT::UTILS::GaussRule3D
      rule =  DRT::UTILS::intrule3D_undefined;}; template<> struct
      DisTypeToGaussRuleForExactSol<DRT::Element::quad4>   {static const DRT::UTILS::GaussRule2D
      rule =  DRT::UTILS::intrule_quad_9point;}; template<> struct
      DisTypeToGaussRuleForExactSol<DRT::Element::quad8>   {static const DRT::UTILS::GaussRule2D
      rule =  DRT::UTILS::intrule_quad_9point;}; template<> struct
      DisTypeToGaussRuleForExactSol<DRT::Element::quad9>   {static const DRT::UTILS::GaussRule2D
      rule =  DRT::UTILS::intrule_quad_9point;}; template<> struct
      DisTypeToGaussRuleForExactSol<DRT::Element::tri3>    {static const DRT::UTILS::GaussRule2D
      rule =  DRT::UTILS::intrule2D_undefined;}; template<> struct
      DisTypeToGaussRuleForExactSol<DRT::Element::tri6>    {static const DRT::UTILS::GaussRule2D
      rule =  DRT::UTILS::intrule2D_undefined;}; template<> struct
      DisTypeToGaussRuleForExactSol<DRT::Element::line2>   {static const DRT::UTILS::GaussRule1D
      rule =  DRT::UTILS::intrule1D_undefined; }; template<> struct
      DisTypeToGaussRuleForExactSol<DRT::Element::line3>   {static const DRT::UTILS::GaussRule1D
      rule =  DRT::UTILS::intrule1D_undefined; };
    */
    //! provide parameter m_k for each distype required for certain stabilization parameter
    //! definitions
    template <DRT::Element::DiscretizationType DISTYPE>
    inline double MK()
    {
      switch (DISTYPE)
      {
        case DRT::Element::tet4:
        case DRT::Element::pyramid5:
        case DRT::Element::hex8:
        case DRT::Element::wedge6:
        case DRT::Element::nurbs8:
          return 0.333333333333333333333;
          break;
        case DRT::Element::hex20:
        case DRT::Element::hex27:
        case DRT::Element::tet10:
        case DRT::Element::wedge15:
        case DRT::Element::nurbs27:
          return 0.083333333333333333333;
          break;
          // do the 2D case after the more usual 3D case
        case DRT::Element::tri3:
        case DRT::Element::quad4:
        case DRT::Element::nurbs4:
          return 0.333333333333333333333;
          break;
        case DRT::Element::tri6:
        case DRT::Element::quad8:
        case DRT::Element::quad9:
        case DRT::Element::nurbs9:
          return 0.083333333333333333333;
          break;
          // finally, do the 1D case
        case DRT::Element::line2:
        case DRT::Element::nurbs2:
          return 0.333333333333333333333;
          break;
        case DRT::Element::line3:
        case DRT::Element::nurbs3:
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
