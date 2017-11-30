/*!
\file scatra_ele_calc_utils.H

\brief Utility methods for scatra

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef SCATRA_ELE_CALC_UTILS_H
#define SCATRA_ELE_CALC_UTILS_H


#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_element.H"
#include "../drt_inpar/inpar_scatra.H"


namespace SCATRA
{
  /*!
  \brief Decide, whether second derivatives are needed  (template version)
   *  In convection-diffusion problems, ONLY N,xx , N,yy and N,zz are needed
   *  to evaluate the laplacian operator for the residual-based stabilization.
   *  Hence, unlike to the Navier-Stokes equations, hex8, wedge6 and pyramid5
   *  return false although they have non-zero MIXED second derivatives.*/
  template<DRT::Element::DiscretizationType DISTYPE>
  struct Use2ndDerivs {};
  template<> struct Use2ndDerivs<DRT::Element::hex8>    {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::tet4>    {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::wedge6>  {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::pyramid5>{static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::nurbs8>  {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::quad4>   {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::nurbs4>  {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::tri3>    {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::line2>   {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::nurbs2>  {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::hex20>   {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::hex27>   {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::nurbs27> {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::tet10>   {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::quad8>   {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::quad9>   {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::nurbs9>  {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::tri6>    {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::line3>   {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::nurbs3>  {static const bool use = true;};

  //! Template Meta Programming version of switch over discretization type
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToOptGaussRule {};
  template<> struct DisTypeToOptGaussRule<DRT::Element::hex8>    {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_8point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::hex20>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::hex27>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::tet4>    {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_tet_4point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::tet10>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_tet_5point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::wedge6>  {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_wedge_6point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::pyramid5>{static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_pyramid_8point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::nurbs8>  {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_8point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::nurbs27> {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::quad4>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_4point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::quad8>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::quad9>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::tri3>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_tri_3point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::tri6>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_tri_6point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::nurbs4>  {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_4point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::nurbs9>  {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::line2>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_2point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::line3>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_3point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::nurbs2>  {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_2point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::nurbs3>  {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_3point; };

  //! Template Meta Programming version of switch over discretization type
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToMatGaussRule {};
  template<> struct DisTypeToMatGaussRule<DRT::Element::hex8>
  {
    static DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      switch(degree)
      {
      case 0:
      case 1:
        return DRT::UTILS::intrule_hex_1point;
        break;
      case 2:
      case 3:
        return DRT::UTILS::intrule_hex_8point;
        break;
      case 4:
      case 5:
        return DRT::UTILS::intrule_hex_27point;
        break;
      case 6:
      case 7:
        return DRT::UTILS::intrule_hex_64point;
        break;
      case 8:
      case 9:
        return DRT::UTILS::intrule_hex_125point;
        break;
      case 10:
      case 11:
        return DRT::UTILS::intrule_hex_216point;
        break;
      case 12:
      case 13:
        return DRT::UTILS::intrule_hex_343point;
        break;
      case 14:
      case 15:
        return DRT::UTILS::intrule_hex_512point;
        break;
      case 16:
      case 17:
        return DRT::UTILS::intrule_hex_729point;
        break;
      case 18:
      case 19:
        return DRT::UTILS::intrule_hex_1000point;
        break;
      default:
        dserror("Integration rule only until degree 19 for HEX elements defined. You used a degree of %d",degree);
        return DRT::UTILS::intrule3D_undefined;
        break;
      }
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::hex20>
  {
    static DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      switch(degree)
      {
      case 0:
      case 1:
        return DRT::UTILS::intrule_hex_1point;
        break;
      case 2:
      case 3:
        return DRT::UTILS::intrule_hex_8point;
        break;
      case 4:
      case 5:
        return DRT::UTILS::intrule_hex_27point;
        break;
      case 6:
      case 7:
        return DRT::UTILS::intrule_hex_64point;
        break;
      case 8:
      case 9:
        return DRT::UTILS::intrule_hex_125point;
        break;
      case 10:
      case 11:
        return DRT::UTILS::intrule_hex_216point;
        break;
      case 12:
      case 13:
        return DRT::UTILS::intrule_hex_343point;
        break;
      case 14:
      case 15:
        return DRT::UTILS::intrule_hex_512point;
        break;
      case 16:
      case 17:
        return DRT::UTILS::intrule_hex_729point;
        break;
      case 18:
      case 19:
        return DRT::UTILS::intrule_hex_1000point;
        break;
      default:
        dserror("Integration rule only until degree 19 for HEX elements defined. You used a degree of %d",degree);
        return DRT::UTILS::intrule3D_undefined;
        break;
      }
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::hex27>
  {
    static DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      switch(degree)
      {
      case 0:
      case 1:
        return DRT::UTILS::intrule_hex_1point;
        break;
      case 2:
      case 3:
        return DRT::UTILS::intrule_hex_8point;
        break;
      case 4:
      case 5:
        return DRT::UTILS::intrule_hex_27point;
        break;
      case 6:
      case 7:
        return DRT::UTILS::intrule_hex_64point;
        break;
      case 8:
      case 9:
        return DRT::UTILS::intrule_hex_125point;
        break;
      case 10:
      case 11:
        return DRT::UTILS::intrule_hex_216point;
        break;
      case 12:
      case 13:
        return DRT::UTILS::intrule_hex_343point;
        break;
      case 14:
      case 15:
        return DRT::UTILS::intrule_hex_512point;
        break;
      case 16:
      case 17:
        return DRT::UTILS::intrule_hex_729point;
        break;
      case 18:
      case 19:
        return DRT::UTILS::intrule_hex_1000point;
        break;
      default:
        dserror("Integration rule only until degree 19 for HEX elements defined. You used a degree of %d",degree);
        return DRT::UTILS::intrule3D_undefined;
        break;
      }
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::tet4>
  {
    static DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      switch(degree)
      {
      case 0:
      case 1:
        return DRT::UTILS::intrule_tet_1point;
        break;
      case 2:
        return DRT::UTILS::intrule_tet_4point;
        break;
      case 3:
        return DRT::UTILS::intrule_tet_5point;
        break;
      case 4:
        return DRT::UTILS::intrule_tet_11point;
        break;
      case 5:
        return DRT::UTILS::intrule_tet_15point;
        break;
      case 6:
        return DRT::UTILS::intrule_tet_24point;
        break;
      case 7:
      case 8:
        return DRT::UTILS::intrule_tet_45point;
        break;
      case 9:
        return DRT::UTILS::intrule_tet_125point_peano;
        break;
      case 10:
      case 11:
      case 12:
      case 13:
        return DRT::UTILS::intrule_tet_343point_peano;
        break;
      case 14:
      case 15:
      case 16:
        return DRT::UTILS::intrule_tet_729point_peano;
        break;
      default:
        dserror("Integration rule only until degree 16 for TET elements defined. You used a degree of %d",degree);
        return DRT::UTILS::intrule3D_undefined;
        break;
      }
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::tet10>
  {
    static DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      switch(degree)
      {
      case 0:
      case 1:
        return DRT::UTILS::intrule_tet_1point;
        break;
      case 2:
        return DRT::UTILS::intrule_tet_4point;
        break;
      case 3:
        return DRT::UTILS::intrule_tet_5point;
        break;
      case 4:
        return DRT::UTILS::intrule_tet_11point;
        break;
      case 5:
        return DRT::UTILS::intrule_tet_15point;
        break;
      case 6:
        return DRT::UTILS::intrule_tet_24point;
        break;
      case 7:
      case 8:
        return DRT::UTILS::intrule_tet_45point;
        break;
      case 9:
        return DRT::UTILS::intrule_tet_125point_peano;
        break;
      case 10:
      case 11:
      case 12:
      case 13:
        return DRT::UTILS::intrule_tet_343point_peano;
        break;
      case 14:
      case 15:
      case 16:
        return DRT::UTILS::intrule_tet_729point_peano;
        break;
      default:
        dserror("Integration rule only until degree 16 for TET elements defined. You used a degree of %d",degree);
        return DRT::UTILS::intrule3D_undefined;
        break;
      }
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::wedge6>
  {
    static DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for WEDGE elements defined. Feel free to add the Gauss rule.");
      return DRT::UTILS::intrule3D_undefined;
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::pyramid5>
  {
    static DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for PYRAMID elements defined. Feel free to add the Gauss rule.");
      return DRT::UTILS::intrule3D_undefined;
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::nurbs8>
  {
    static DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return DRT::UTILS::intrule3D_undefined;
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::nurbs27>
  {
    static DRT::UTILS::GaussRule3D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return DRT::UTILS::intrule3D_undefined;
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::quad4>
  {
    static DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      switch(degree)
      {
      case 0:
      case 1:
        return DRT::UTILS::intrule_quad_1point;
        break;
      case 2:
      case 3:
        return DRT::UTILS::intrule_quad_4point;
        break;
      case 4:
      case 5:
        return DRT::UTILS::intrule_quad_9point;
        break;
      case 6:
      case 7:
        return DRT::UTILS::intrule_quad_16point;
        break;
      case 8:
      case 9:
        return DRT::UTILS::intrule_quad_25point;
        break;
      case 10:
      case 11:
        return DRT::UTILS::intrule_quad_36point;
        break;
      case 12:
      case 13:
        return DRT::UTILS::intrule_quad_49point;
        break;
      case 14:
      case 15:
        return DRT::UTILS::intrule_quad_64point;
        break;
      case 16:
      case 17:
        return DRT::UTILS::intrule_quad_81point;
        break;
      case 18:
      case 19:
        return DRT::UTILS::intrule_quad_100point;
        break;
      default:
        dserror("Integration rule only until degree 19 for QUAD elements defined. You used a degree of %d",degree);
        return DRT::UTILS::intrule2D_undefined;
        break;
      }
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::quad8>
  {
    static DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      switch(degree)
      {
      case 0:
      case 1:
        return DRT::UTILS::intrule_quad_1point;
        break;
      case 2:
      case 3:
        return DRT::UTILS::intrule_quad_4point;
        break;
      case 4:
      case 5:
        return DRT::UTILS::intrule_quad_9point;
        break;
      case 6:
      case 7:
        return DRT::UTILS::intrule_quad_16point;
        break;
      case 8:
      case 9:
        return DRT::UTILS::intrule_quad_25point;
        break;
      case 10:
      case 11:
        return DRT::UTILS::intrule_quad_36point;
        break;
      case 12:
      case 13:
        return DRT::UTILS::intrule_quad_49point;
        break;
      case 14:
      case 15:
        return DRT::UTILS::intrule_quad_64point;
        break;
      case 16:
      case 17:
        return DRT::UTILS::intrule_quad_81point;
        break;
      case 18:
      case 19:
        return DRT::UTILS::intrule_quad_100point;
        break;
      default:
        dserror("Integration rule only until degree 19 for QUAD elements defined. You used a degree of %d",degree);
        return DRT::UTILS::intrule2D_undefined;
        break;
      }
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::quad9>
  {
    static DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      switch(degree)
      {
      case 0:
      case 1:
        return DRT::UTILS::intrule_quad_1point;
        break;
      case 2:
      case 3:
        return DRT::UTILS::intrule_quad_4point;
        break;
      case 4:
      case 5:
        return DRT::UTILS::intrule_quad_9point;
        break;
      case 6:
      case 7:
        return DRT::UTILS::intrule_quad_16point;
        break;
      case 8:
      case 9:
        return DRT::UTILS::intrule_quad_25point;
        break;
      case 10:
      case 11:
        return DRT::UTILS::intrule_quad_36point;
        break;
      case 12:
      case 13:
        return DRT::UTILS::intrule_quad_49point;
        break;
      case 14:
      case 15:
        return DRT::UTILS::intrule_quad_64point;
        break;
      case 16:
      case 17:
        return DRT::UTILS::intrule_quad_81point;
        break;
      case 18:
      case 19:
        return DRT::UTILS::intrule_quad_100point;
        break;
      default:
        dserror("Integration rule only until degree 19 for QUAD elements defined. You used a degree of %d",degree);
        return DRT::UTILS::intrule2D_undefined;
        break;
      }
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::tri3>
  {
    static DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      switch(degree)
      {
      case 0:
      case 1:
        return DRT::UTILS::intrule_tri_1point;
        break;
      case 2:
        return DRT::UTILS::intrule_tri_3point;
        break;
      case 3:
        return DRT::UTILS::intrule_tri_4point;
        break;
      case 4:
        return DRT::UTILS::intrule_tri_6point;
        break;
      case 5:
        return DRT::UTILS::intrule_tri_7point;
        break;
      case 6:
        return DRT::UTILS::intrule_tri_12point;
        break;
      case 7:
      case 8:
        return DRT::UTILS::intrule_tri_16point;
        break;
      case 9:
      case 10:
      case 11:
      case 12:
      case 13:
        return DRT::UTILS::intrule_tri_37point;
        break;
      case 14:
      case 15:
        return DRT::UTILS::intrule_tri_64point;
        break;
      default:
        dserror("Integration rule only until degree 15 for TRI elements defined. You used a degree of %d",degree);
        return DRT::UTILS::intrule2D_undefined;
        break;
      }
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::tri6>
  {
    static DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      return DRT::UTILS::intrule_tri_16point;
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::nurbs4>
  {
    static DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return DRT::UTILS::intrule2D_undefined;
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::nurbs9>
  {
    static DRT::UTILS::GaussRule2D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return DRT::UTILS::intrule2D_undefined;
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::line2>
  {
    static DRT::UTILS::GaussRule1D GetGaussRule(int degree)
    {
      switch(degree)
      {
      case 0:
      case 1:
        return DRT::UTILS::intrule_line_1point;
        break;
      case 2:
      case 3:
        return DRT::UTILS::intrule_line_2point;
        break;
      case 4:
      case 5:
        return DRT::UTILS::intrule_line_3point;
        break;
      case 6:
      case 7:
        return DRT::UTILS::intrule_line_4point;
        break;
      case 8:
      case 9:
        return DRT::UTILS::intrule_line_5point;
        break;
      case 10:
      case 11:
        return DRT::UTILS::intrule_line_6point;
        break;
      case 12:
      case 13:
        return DRT::UTILS::intrule_line_7point;
        break;
      case 14:
      case 15:
        return DRT::UTILS::intrule_line_8point;
        break;
      case 16:
      case 17:
        return DRT::UTILS::intrule_line_9point;
        break;
      case 18:
      case 19:
        return DRT::UTILS::intrule_line_10point;
        break;
      default:
        dserror("Integration rule only until degree 19 for LINE elements defined. You used a degree of %d",degree);
        return DRT::UTILS::intrule1D_undefined;
        break;
      }
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::line3>
  {
    static DRT::UTILS::GaussRule1D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for LINE elements defined. Feel free to add the Gauss rule.");
      return DRT::UTILS::intrule1D_undefined;
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::nurbs2>
  {
    static DRT::UTILS::GaussRule1D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return DRT::UTILS::intrule1D_undefined;
    };
  };
  template<> struct DisTypeToMatGaussRule<DRT::Element::nurbs3>
  {
    static DRT::UTILS::GaussRule1D GetGaussRule(int degree)
    {
      dserror("Gauss rule not for NURBS elements defined. Feel free to add the Gauss rule.");
      return DRT::UTILS::intrule1D_undefined;
    };
  };


  //! Template Meta Programming version of switch over discretization type
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToGaussRuleForExactSol {};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::hex8>    {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::hex20>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::hex27>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::tet4>    {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_tet_5point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::tet10>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule3D_undefined;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::wedge6>  {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule3D_undefined;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::pyramid5>{static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule3D_undefined;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs8>  {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point; };
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs27> {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::quad4>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::quad8>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::quad9>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::tri3>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule2D_undefined;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::tri6>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule2D_undefined;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs4>  {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_4point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs9>  {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::line2>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_2point; };
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::line3>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_8point; };
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs2>  {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule1D_undefined; };
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs3>  {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule1D_undefined; };

  //! for each distype provide the m_k needed for stabilization computation
  template<DRT::Element::DiscretizationType DISTYPE>
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

  //! Template Meta Programming version of switch over discretization type
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToStabGaussRule {};
  template<> struct DisTypeToStabGaussRule<DRT::Element::hex8>    {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_1point; };
  template<> struct DisTypeToStabGaussRule<DRT::Element::hex20>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_1point;};
  template<> struct DisTypeToStabGaussRule<DRT::Element::hex27>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_1point;};
  template<> struct DisTypeToStabGaussRule<DRT::Element::tet4>    {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_tet_1point; };
  template<> struct DisTypeToStabGaussRule<DRT::Element::tet10>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_tet_1point; };
  template<> struct DisTypeToStabGaussRule<DRT::Element::wedge6>  {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_wedge_1point;};
  template<> struct DisTypeToStabGaussRule<DRT::Element::pyramid5>{static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_pyramid_1point;};
  template<> struct DisTypeToStabGaussRule<DRT::Element::nurbs8>  {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_1point; };
  template<> struct DisTypeToStabGaussRule<DRT::Element::nurbs27> {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_1point; };
  template<> struct DisTypeToStabGaussRule<DRT::Element::quad4>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_1point;};
  template<> struct DisTypeToStabGaussRule<DRT::Element::quad8>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_1point;};
  template<> struct DisTypeToStabGaussRule<DRT::Element::quad9>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_1point;};
  template<> struct DisTypeToStabGaussRule<DRT::Element::nurbs4>  {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_1point;};
  template<> struct DisTypeToStabGaussRule<DRT::Element::nurbs9>  {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_1point;};
  template<> struct DisTypeToStabGaussRule<DRT::Element::tri3>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_tri_1point; };
  template<> struct DisTypeToStabGaussRule<DRT::Element::tri6>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_tri_1point; };
  template<> struct DisTypeToStabGaussRule<DRT::Element::line2>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_1point; };
  template<> struct DisTypeToStabGaussRule<DRT::Element::line3>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_1point; };
  template<> struct DisTypeToStabGaussRule<DRT::Element::nurbs2>  {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_1point; };
  template<> struct DisTypeToStabGaussRule<DRT::Element::nurbs3>  {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_1point; };


  //! determine whether there are only two charged species present or not
  bool IsBinaryElectrolyte(const std::vector<double>& valence);

  //! determine indices of the two charged species in case of an binary electrolyte
  std::vector<int> GetIndicesBinaryElectrolyte(const std::vector<double>& valence);

  //! calculate resultant diffusion coefficient for stabilization of binary electrolyte systems
  double CalResDiffCoeff(
      const std::vector<double>& valence, ///< valences
      const std::vector<double>& diffus,  ///< diffusivities
      const std::vector<int>& indices  ///< indices
  );

  //! identify elements of inflow section
  bool InflowElement(const DRT::Element* ele);

  //! convert implementation type of scalar transport elements into corresponding string for output purposes
  const std::string ImplTypeToString(const INPAR::SCATRA::ImplType impltype);
} // namespace SCATRA

#endif
