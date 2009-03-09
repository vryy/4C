/*!
\file scatra_ele_impl_utils.H

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef SCATRA_ELE_IMPL_UTILS_H
#define SCATRA_ELE_IMPL_UTILS_H


#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_element.H"


namespace SCATRA
{

  //! possible options for the stabilization parameter
  enum TauType
  {
    tau_franca_valentin,
    tau_bazilevs,
    tau_zero
  };

  //! flux parameters recognized by scalar transport implementation class
  enum FluxType
  {
    noflux,
    totalflux,
    diffusiveflux
  };

  //! decide, whether second derivatives are needed  (template version)
  /*! In convection-diffusion problems, ONLY N,xx , N,yy and N,zz are needed
   *  to evaluate the laplacian operator for the residual-based stabilization.
   *  Hence, unlike to the Navier-Stokes equations, hex8, wedge6 and pyramid5 
   *  return false although they have non-zero MIXED second derivatives.*/
  template<DRT::Element::DiscretizationType DISTYPE>
  struct Use2ndDerivs {};
  template<> struct Use2ndDerivs<DRT::Element::hex8>    {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::tet4>    {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::wedge6>  {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::pyramid5>{static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::quad4>   {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::tri3>    {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::line2>   {static const bool use = false;};
  template<> struct Use2ndDerivs<DRT::Element::hex20>   {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::hex27>   {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::tet10>   {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::quad8>   {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::quad9>   {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::tri6>    {static const bool use = true; };
  template<> struct Use2ndDerivs<DRT::Element::line3>   {static const bool use = true; };


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
  template<> struct DisTypeToOptGaussRule<DRT::Element::quad4>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_4point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::quad8>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::quad9>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToOptGaussRule<DRT::Element::tri3>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_tri_3point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::tri6>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_tri_6point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::line2>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_2point; };
  template<> struct DisTypeToOptGaussRule<DRT::Element::line3>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_3point; };

  //! Template Meta Programming version of switch over discretization type
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToGaussRuleForExactSol {};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::hex8>    {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::hex20>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::hex27>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::tet4>    {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule3D_undefined;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::tet10>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule3D_undefined;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::wedge6>  {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule3D_undefined;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::pyramid5>{static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule3D_undefined;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::quad4>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::quad8>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::quad9>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_9point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::tri3>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule2D_undefined;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::tri6>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule2D_undefined;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::line2>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule1D_undefined; };
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::line3>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule1D_undefined; };


  //! for each distype provide the m_k needed for stabilization computation
  template<DRT::Element::DiscretizationType DISTYPE>
  struct DisTypeToMK {};
  template<> struct DisTypeToMK<DRT::Element::hex8>    {static const double mk = 0.333333333333333333333;};
  template<> struct DisTypeToMK<DRT::Element::tet4>    {static const double mk = 0.333333333333333333333;};
  template<> struct DisTypeToMK<DRT::Element::wedge6>  {static const double mk = 0.333333333333333333333;};
  template<> struct DisTypeToMK<DRT::Element::pyramid5>{static const double mk = 0.333333333333333333333;};
  template<> struct DisTypeToMK<DRT::Element::quad4>   {static const double mk = 0.333333333333333333333;};
  template<> struct DisTypeToMK<DRT::Element::tri3>    {static const double mk = 0.333333333333333333333;};
  template<> struct DisTypeToMK<DRT::Element::line2>   {static const double mk = 0.333333333333333333333;};

  template<> struct DisTypeToMK<DRT::Element::hex20>   {static const double mk = 0.083333333333333333333;};
  template<> struct DisTypeToMK<DRT::Element::hex27>   {static const double mk = 0.083333333333333333333;};
  template<> struct DisTypeToMK<DRT::Element::tet10>   {static const double mk = 0.083333333333333333333;};
  template<> struct DisTypeToMK<DRT::Element::quad8>   {static const double mk = 0.083333333333333333333;};
  template<> struct DisTypeToMK<DRT::Element::quad9>   {static const double mk = 0.083333333333333333333;};
  template<> struct DisTypeToMK<DRT::Element::tri6>    {static const double mk = 0.083333333333333333333;};
  template<> struct DisTypeToMK<DRT::Element::line3>   {static const double mk = 0.083333333333333333333;};


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
  template<> struct DisTypeToStabGaussRule<DRT::Element::quad4>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_1point;};
  template<> struct DisTypeToStabGaussRule<DRT::Element::quad8>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_1point;};
  template<> struct DisTypeToStabGaussRule<DRT::Element::quad9>   {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_quad_1point;};
  template<> struct DisTypeToStabGaussRule<DRT::Element::tri3>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_tri_1point; };
  template<> struct DisTypeToStabGaussRule<DRT::Element::tri6>    {static const DRT::UTILS::GaussRule2D rule =  DRT::UTILS::intrule_tri_1point; };
  template<> struct DisTypeToStabGaussRule<DRT::Element::line2>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_1point; };
  template<> struct DisTypeToStabGaussRule<DRT::Element::line3>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_1point; };


} // namespace SCATRA

#endif
#endif  // #ifdef CCADISCRET
