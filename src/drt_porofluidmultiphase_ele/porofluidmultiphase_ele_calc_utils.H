/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_ele_calc_utils.H

 \brief helpful template definitions for the porofluidmultiphase element

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_CALC_UTILS_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_CALC_UTILS_H_

#include "../drt_lib/drt_element.H"
#include "../drt_fem_general/drt_utils_integration.H"

namespace POROFLUIDMULTIPHASE
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
  struct DisTypeToGaussRuleForExactSol {};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::hex8>    {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::hex20>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::hex27>   {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule_hex_27point;};
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::tet4>    {static const DRT::UTILS::GaussRule3D rule =  DRT::UTILS::intrule3D_undefined;};
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
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::line2>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule_line_2point; }; // not tested
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::line3>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule1D_undefined; };
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs2>  {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule1D_undefined; };
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::nurbs3>  {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule1D_undefined; };


} // namespace POROFLUIDMULTIPHASE



#endif /* SRC_DRT_POROFLUIDMULTIPHASE_ELE_POROFLUIDMULTIPHASE_ELE_CALC_UTILS_H_ */
