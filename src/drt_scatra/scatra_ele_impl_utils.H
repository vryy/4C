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

#ifndef SCATRA_ELE_IMPL_UTILS_H
#define SCATRA_ELE_IMPL_UTILS_H


#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_element.H"
#include "../drt_inpar/inpar_scatra.H"


namespace SCATRA
{

  //! does the given scatratype belong to the ELCH group of problems?
  bool IsElchProblem(const enum INPAR::SCATRA::ScaTraType scatratype);

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
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::line2>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule1D_undefined; };
  template<> struct DisTypeToGaussRuleForExactSol<DRT::Element::line3>   {static const DRT::UTILS::GaussRule1D rule =  DRT::UTILS::intrule1D_undefined; };
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

} // namespace SCATRA

#endif
