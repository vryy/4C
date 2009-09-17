/*!
\file thermo_ele_impl_utils.H

<pre>
Maintainer: Caroline Danowski
            danowski@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef THERMO_ELE_IMPL_UTILS_H
#define THERMO_ELE_IMPL_UTILS_H


#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_element.H"


namespace THR
{

  //! flux parameters recognized by thermo implementation class (Temper)
  enum FluxType
  {
    noflux,
    totalflux,
    diffusiveflux
  };

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
  inline double MK()
  {
    switch (DISTYPE)
    {
    case DRT::Element::tet4:
    case DRT::Element::pyramid5:
    case DRT::Element::hex8:
    case DRT::Element::wedge6:
      return 0.333333333333333333333;
      break;
    case DRT::Element::hex20:
    case DRT::Element::hex27:
    case DRT::Element::tet10:
    case DRT::Element::wedge15:
      return 0.083333333333333333333;
      break;
      // do the 2D case after the more usual 3D case
    case DRT::Element::tri3:
    case DRT::Element::quad4:
      return 0.333333333333333333333;
      break;
    case DRT::Element::tri6:
    case DRT::Element::quad8:
    case DRT::Element::quad9:
      return 0.083333333333333333333;
      break;
      // finally, do the 1D case
    case DRT::Element::line2:
      return 0.333333333333333333333;
      break;
    case DRT::Element::line3:
      return 0.083333333333333333333;
      break;
    default:
      dserror("Element shape not supported.");
    }
    return -1.0;
  }

} // namespace THR

#endif  // THERMO_ELE_IMPL_UTILS_H
#endif  // CCADISCRET
