/*----------------------------------------------------------------------*/
/*!
 * \file fluid3_stabilization.H
 * 
 * \brief stabilization related, templated functions
 * 
 * <pre>
 * Maintainer: Axel Gerstenberger
 *             gerstenberger@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15236
 * </pre>
 */
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef XFLUID3_STABILIZATION_H
#define XFLUID3_STABILIZATION_H

#include <blitz/array.h>
#include "../drt_lib/drt_utils.H"
#include "../drt_xfem/xfem.H"


namespace XFLUID
{
  
  using namespace XFEM;
  
  //! for each distype there should be an m_k needed for stabilization computation
  template<DRT::Element::DiscretizationType DISTYPE>
  double MK()
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
      default:
        dserror("type unknown!\n");
    }   
    return -1.0;
  }
  
  
  //! get one point gauss rule to calculate tau at element center
  template<DRT::Element::DiscretizationType DISTYPE>
  DRT::UTILS::GaussRule3D getIntegrationRuleForStabilization()
  {
    switch (DISTYPE)
    {
      case DRT::Element::hex8:
      case DRT::Element::hex20:
      case DRT::Element::hex27:
        return DRT::UTILS::intrule_hex_1point;
        break;
      case DRT::Element::tet4:
      case DRT::Element::tet10:
        return DRT::UTILS::intrule_tet_1point;
        break;
      case DRT::Element::wedge6:
      case DRT::Element::wedge15:
        return DRT::UTILS::intrule_wedge_1point;
        break;
      case DRT::Element::pyramid5:
        return DRT::UTILS::intrule_pyramid_1point;
        break;
      default:
        dserror("invalid discretization type for fluid3");
        return DRT::UTILS::intrule3D_undefined;
    }
  }
  
  
  //! calculate characteristic element size for stabilization
  template<DRT::Element::DiscretizationType DISTYPE, class M1>
  double HK(
      const BlitzMat&                         evelnp,
      const double                            visc,
      const M1&                               xyze
  )
  {
    const int numnode = DRT::UTILS::_switchDisType<DISTYPE>::numNodePerElement;
    
    // use one point gauss rule to calculate tau at element center
    const DRT::UTILS::GaussRule3D integrationrule_stabili = getIntegrationRuleForStabilization<DISTYPE>();
    
    // gaussian points
    const DRT::UTILS::IntegrationPoints3D intpoints(integrationrule_stabili);
    
    // shape functions and derivs at element center
    static BlitzVec3 e;
    e(0) = intpoints.qxg[0][0];
    e(1) = intpoints.qxg[0][1];
    e(2) = intpoints.qxg[0][2];
    const double wquad = intpoints.qwgt[0];
    
    static blitz::TinyMatrix<double,3,numnode> deriv_;
    DRT::UTILS::shape_function_deriv1<DISTYPE>(e, deriv_);
    
    // get Jacobian matrix and determinant
    //const BlitzMat xjm_(blitz::sum(deriv_(i,k)*xyze(j,k),k));
    static BlitzMat3x3 xjm_;
    BLITZTINY::MMt_product<3,3,numnode>(deriv_,xyze,xjm_);
    const double det = xjm_(0,0)*xjm_(1,1)*xjm_(2,2)+
    xjm_(0,1)*xjm_(1,2)*xjm_(2,0)+
    xjm_(0,2)*xjm_(1,0)*xjm_(2,1)-
    xjm_(0,2)*xjm_(1,1)*xjm_(2,0)-
    xjm_(0,0)*xjm_(1,2)*xjm_(2,1)-
    xjm_(0,1)*xjm_(1,0)*xjm_(2,2);
    const double vol = wquad*det;
    
    // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
    const double hk = pow((6.*vol/PI),(1.0/3.0))/sqrt(3.0);
    
    return hk;
  }
  
}  // end namespace XFLUID

#endif

#endif
#endif
