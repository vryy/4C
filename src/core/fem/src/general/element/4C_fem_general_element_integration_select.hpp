/*---------------------------------------------------------------------*/
/*! \file

\brief templated version for selecting integration points from element type

\level 1


*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FEM_GENERAL_ELEMENT_INTEGRATION_SELECT_HPP
#define FOUR_C_FEM_GENERAL_ELEMENT_INTEGRATION_SELECT_HPP


#include "4C_config.hpp"

#include "4C_fem_general_element.hpp"
#include "4C_fem_general_utils_integration.hpp"

FOUR_C_NAMESPACE_OPEN


namespace Discret
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
    template <Core::FE::CellType distype>
    struct DisTypeToOptGaussRule
    {
    };
    template <>
    struct DisTypeToOptGaussRule<Core::FE::CellType::hex8>
    {
      static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_8point;
    };

    template <>
    struct DisTypeToOptGaussRule<Core::FE::CellType::hex18>
    {
      static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_18point;
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
    struct DisTypeToOptGaussRule<Core::FE::CellType::wedge15>
    {
      static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::wedge_9point;
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

    //! one-point integration rule should be exact for all elements with straight lines and surfaces
    //! Otherwise, it is assumed as an accurate enough approximation

    //! Template Meta Programming version of switch over discretization type
    template <Core::FE::CellType distype>
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
    struct DisTypeToStabGaussRule<Core::FE::CellType::hex16>
    {
      static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::hex_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<Core::FE::CellType::hex18>
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
    struct DisTypeToStabGaussRule<Core::FE::CellType::wedge15>
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
    struct DisTypeToStabGaussRule<Core::FE::CellType::line2>
    {
      static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_1point;
    };
    template <>
    struct DisTypeToStabGaussRule<Core::FE::CellType::line3>
    {
      static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::line_1point;
    };

    //! wedge6 and pyramid5 have nonzero 2nd mixed derivatives, but they are neglected
    //! quad4 and hex8 are pseudo higherOrder elements
    //! Template Meta Programming version of switch over discretization type
    template <Core::FE::CellType distype>
    struct IsHigherOrder
    {
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::hex8>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::hex20>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::hex27>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::tet4>
    {
      static constexpr bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::tet10>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::wedge6>
    {
      static constexpr bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::wedge15>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::pyramid5>
    {
      static constexpr bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::quad4>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::quad8>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::quad9>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::tri3>
    {
      static constexpr bool ishigherorder = false;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::tri6>
    {
      static constexpr bool ishigherorder = true;
    };
    // template<> struct IsHigherOrder<Core::FE::CellType::line2>   {static const bool
    // ishigherorder = false; }; template<> struct
    // IsHigherOrder<Core::FE::CellType::line3>   {static constexpr bool ishigherorder =
    // false; }; template<> struct IsHigherOrder<Core::FE::CellType::nurbs2> {static
    // const bool ishigherorder = true; }; template<> struct
    // IsHigherOrder<Core::FE::CellType::nurbs3> {static constexpr bool ishigherorder =
    // true; };
    template <>
    struct IsHigherOrder<Core::FE::CellType::nurbs4>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::nurbs9>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::nurbs8>
    {
      static constexpr bool ishigherorder = true;
    };
    template <>
    struct IsHigherOrder<Core::FE::CellType::nurbs27>
    {
      static constexpr bool ishigherorder = true;
    };


    //! Template Meta Programming version of switch over discretization type
    template <Core::FE::CellType distype>
    struct IsNurbs
    {
    };
    template <>
    struct IsNurbs<Core::FE::CellType::hex8>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::hex20>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::hex27>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::tet4>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::tet10>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::wedge6>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::wedge15>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::pyramid5>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::quad4>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::quad8>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::quad9>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::tri3>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::tri6>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::line2>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::line3>
    {
      static constexpr bool isnurbs = false;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::nurbs2>
    {
      static constexpr bool isnurbs = true;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::nurbs3>
    {
      static constexpr bool isnurbs = true;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::nurbs4>
    {
      static constexpr bool isnurbs = true;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::nurbs9>
    {
      static constexpr bool isnurbs = true;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::nurbs8>
    {
      static constexpr bool isnurbs = true;
    };
    template <>
    struct IsNurbs<Core::FE::CellType::nurbs27>
    {
      static constexpr bool isnurbs = true;
    };


    template <Core::FE::CellType distype>
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
      static constexpr Core::FE::GaussRule3D rule = Core::FE::GaussRule3D::undefined;
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
      static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::undefined;
    };
    template <>
    struct DisTypeToGaussRuleForExactSol<Core::FE::CellType::line3>
    {
      static constexpr Core::FE::GaussRule1D rule = Core::FE::GaussRule1D::undefined;
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
    //! provide parameter m_k for each distype required for certain stabilization parameter
    //! definitions
    template <Core::FE::CellType distype>
    inline double MK()
    {
      switch (distype)
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

  }  // namespace ELEMENTS

}  // namespace Discret

FOUR_C_NAMESPACE_CLOSE

#endif
