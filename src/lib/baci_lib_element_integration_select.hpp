/*---------------------------------------------------------------------*/
/*! \file

\brief templated version for selecting integration points from element type

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_LIB_ELEMENT_INTEGRATION_SELECT_HPP
#define FOUR_C_LIB_ELEMENT_INTEGRATION_SELECT_HPP


#include "baci_config.hpp"

#include "baci_discretization_fem_general_utils_integration.hpp"
#include "baci_lib_element.hpp"

FOUR_C_NAMESPACE_OPEN


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
    template <CORE::FE::CellType distype>
    struct DisTypeToOptGaussRule
    {
    };
    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::hex8>
    {
      static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_8point;
    };

    template <>
    struct DisTypeToOptGaussRule<CORE::FE::CellType::hex18>
    {
      static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_18point;
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
    struct DisTypeToOptGaussRule<CORE::FE::CellType::wedge15>
    {
      static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::wedge_9point;
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

    //! one-point integration rule should be exact for all elements with straight lines and surfaces
    //! Otherwise, it is assumed as an accurate enough approximation

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
    struct DisTypeToStabGaussRule<CORE::FE::CellType::hex16>
    {
      static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::hex_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<CORE::FE::CellType::hex18>
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
    struct DisTypeToStabGaussRule<CORE::FE::CellType::wedge15>
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
    struct DisTypeToStabGaussRule<CORE::FE::CellType::line2>
    {
      static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::line_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<CORE::FE::CellType::line3>
    {
      static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::line_1point;
    };

    //! wedge6 and pyramid5 have nonzero 2nd mixed derivatives, but they are neglected
    //! quad4 and hex8 are pseudo higherOrder elements
    //! Template Meta Programming version of switch over discretization type
    template <CORE::FE::CellType distype>
    struct IsHigherOrder
    {
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::hex8>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::hex20>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::hex27>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::tet4>
    {
      static constexpr bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::tet10>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::wedge6>
    {
      static constexpr bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::wedge15>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::pyramid5>
    {
      static constexpr bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::quad4>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::quad8>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::quad9>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::tri3>
    {
      static constexpr bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::tri6>
    {
      static constexpr bool ishigherorder = true;
    };
    // template<> struct IsHigherOrder<CORE::FE::CellType::line2>   {static const bool
    // ishigherorder = false; }; template<> struct
    // IsHigherOrder<CORE::FE::CellType::line3>   {static constexpr bool ishigherorder =
    // false; }; template<> struct IsHigherOrder<CORE::FE::CellType::nurbs2> {static
    // const bool ishigherorder = true; }; template<> struct
    // IsHigherOrder<CORE::FE::CellType::nurbs3> {static constexpr bool ishigherorder =
    // true; };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::nurbs4>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::nurbs9>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::nurbs8>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<CORE::FE::CellType::nurbs27>
    {
      static constexpr bool ishigherorder = true;
    };


    //! Template Meta Programming version of switch over discretization type
    template <CORE::FE::CellType distype>
    struct IsNurbs
    {
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::hex8>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::hex20>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::hex27>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::tet4>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::tet10>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::wedge6>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::wedge15>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::pyramid5>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::quad4>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::quad8>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::quad9>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::tri3>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::tri6>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::line2>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::line3>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::nurbs2>
    {
      static constexpr bool isnurbs = true;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::nurbs3>
    {
      static constexpr bool isnurbs = true;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::nurbs4>
    {
      static constexpr bool isnurbs = true;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::nurbs9>
    {
      static constexpr bool isnurbs = true;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::nurbs8>
    {
      static constexpr bool isnurbs = true;
    };
    template <>
    struct IsNurbs<CORE::FE::CellType::nurbs27>
    {
      static constexpr bool isnurbs = true;
    };


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
      static constexpr CORE::FE::GaussRule3D rule = CORE::FE::GaussRule3D::undefined;
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
      static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<CORE::FE::CellType::line3>
    {
      static constexpr CORE::FE::GaussRule1D rule = CORE::FE::GaussRule1D::undefined;
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
    //! provide parameter m_k for each distype required for certain stabilization parameter
    //! definitions
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
          FOUR_C_THROW("Element shape not supported.");
          break;
      }
      return -1.0;
    }

  }  // namespace ELEMENTS

}  // namespace DRT

FOUR_C_NAMESPACE_CLOSE

#endif
