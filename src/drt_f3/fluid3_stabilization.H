/*----------------------------------------------------------------------*/
/*!
\file fluid3_stabilization.H

\brief stabilization related, templated functions

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#ifndef XFLUID3_STABILIZATION_H
#define XFLUID3_STABILIZATION_H

#include <blitz/array.h>
#include "../drt_lib/drt_utils.H"


typedef blitz::Array<double,1> BlitzVec;
typedef blitz::Array<double,2> BlitzMat;


namespace XFLUID
{

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
    DRT::UTILS::GaussRule3D integrationrule_stabili = DRT::UTILS::intrule3D_undefined;
    switch (DISTYPE)
    {
    case DRT::Element::hex8:
    case DRT::Element::hex20:
    case DRT::Element::hex27:
      integrationrule_stabili = DRT::UTILS::intrule_hex_1point;
      break;
    case DRT::Element::tet4:
    case DRT::Element::tet10:
      integrationrule_stabili = DRT::UTILS::intrule_tet_1point;
      break;
    case DRT::Element::wedge6:
    case DRT::Element::wedge15:
      integrationrule_stabili = DRT::UTILS::intrule_wedge_1point;
      break;
    case DRT::Element::pyramid5:
      integrationrule_stabili = DRT::UTILS::intrule_pyramid_1point;
      break;
    default:
      dserror("invalid discretization type for fluid3");
    }
    return integrationrule_stabili;
}


//
// calculate stabilization parameter
//
template<DRT::Element::DiscretizationType DISTYPE>
double HK(
        const BlitzMat&                         evelnp,
        const double                            visc,
        const BlitzMat&                         xyze
  )
{
  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex k;    // Placeholder for the third index
  blitz::fourthIndex l;   // Placeholder for the fourth index

  // use one point gauss rule to calculate tau at element center
  const DRT::UTILS::GaussRule3D integrationrule_stabili = getIntegrationRuleForStabilization<DISTYPE>();

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(integrationrule_stabili);

  // shape functions and derivs at element center
  const double e1    = intpoints.qxg[0][0];
  const double e2    = intpoints.qxg[0][1];
  const double e3    = intpoints.qxg[0][2];
  const double wquad = intpoints.qwgt[0];

  const BlitzMat deriv_(DRT::UTILS::shape_function_3D_deriv1<DISTYPE>(e1,e2,e3));

  // get Jacobian matrix and determinant
  const blitz::Array<double,2> xjm_(blitz::sum(deriv_(i,k)*xyze(j,k),k));
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


//
// calculate stabilization parameter
//
template <DRT::Element::DiscretizationType DISTYPE>
BlitzVec CalTauStationary(
  const BlitzMat&                         evelnp,
  const double                            visc,
  const int                               fssgv,
  const BlitzMat&                         xyze
  )
{
  
  BlitzVec tau_(3);       /// tau_Mu, tau_Mp, tau_C
  
  blitz::firstIndex i;    // Placeholder for the first index
  blitz::secondIndex j;   // Placeholder for the second index
  blitz::thirdIndex k;    // Placeholder for the third index
  blitz::fourthIndex l;   // Placeholder for the fourth index

  // use one point gauss rule to calculate tau at element center
  const DRT::UTILS::GaussRule3D integrationrule_stabili = getIntegrationRuleForStabilization<DISTYPE>();

  // gaussian points
  const DRT::UTILS::IntegrationPoints3D intpoints(integrationrule_stabili);

  // shape functions and derivs at element center
  const double e1    = intpoints.qxg[0][0];
  const double e2    = intpoints.qxg[0][1];
  const double e3    = intpoints.qxg[0][2];
  const double wquad = intpoints.qwgt[0];

  const int numnode = DRT::UTILS::getNumberOfElementNodes<DISTYPE>();
  const BlitzVec funct_(DRT::UTILS::shape_function_3D<DISTYPE>(e1,e2,e3));
  const BlitzMat deriv_(DRT::UTILS::shape_function_3D_deriv1<DISTYPE>(e1,e2,e3));

  // get element type constant for tau
  const double mk = MK<DISTYPE>();

  // get velocities at element center
  const BlitzVec velint_(blitz::sum(funct_(j)*evelnp(i,j),j));

  // get Jacobian matrix and determinant
  const BlitzMat xjm(blitz::sum(deriv_(i,k)*xyze(j,k),k));
  const double det = xjm(0,0)*xjm(1,1)*xjm(2,2)+
                     xjm(0,1)*xjm(1,2)*xjm(2,0)+
                     xjm(0,2)*xjm(1,0)*xjm(2,1)-
                     xjm(0,2)*xjm(1,1)*xjm(2,0)-
                     xjm(0,0)*xjm(1,2)*xjm(2,1)-
                     xjm(0,1)*xjm(1,0)*xjm(2,2);
  const double vol = wquad*det;

  // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
  const double hk = pow((6.*vol/PI),(1.0/3.0))/sqrt(3.0);

  // inverse of jacobian
  static BlitzMat xji(3,3);
  xji(0,0) = (  xjm(1,1)*xjm(2,2) - xjm(2,1)*xjm(1,2))/det;
  xji(1,0) = (- xjm(1,0)*xjm(2,2) + xjm(2,0)*xjm(1,2))/det;
  xji(2,0) = (  xjm(1,0)*xjm(2,1) - xjm(2,0)*xjm(1,1))/det;
  xji(0,1) = (- xjm(0,1)*xjm(2,2) + xjm(2,1)*xjm(0,2))/det;
  xji(1,1) = (  xjm(0,0)*xjm(2,2) - xjm(2,0)*xjm(0,2))/det;
  xji(2,1) = (- xjm(0,0)*xjm(2,1) + xjm(2,0)*xjm(0,1))/det;
  xji(0,2) = (  xjm(0,1)*xjm(1,2) - xjm(1,1)*xjm(0,2))/det;
  xji(1,2) = (- xjm(0,0)*xjm(1,2) + xjm(1,0)*xjm(0,2))/det;
  xji(2,2) = (  xjm(0,0)*xjm(1,1) - xjm(1,0)*xjm(0,1))/det;

  // compute global derivates
  const BlitzMat derxy_(blitz::sum(xji(i,k)*deriv_(k,j),k));

  // get velocity norm
  const double vel_norm = sqrt(blitz::sum(velint_*velint_));

  // normed velocity at element centre
  static BlitzVec velino(3);
  if (vel_norm>=1e-6)
  {
    velino = velint_/vel_norm;
  }
  else
  {
    velino = 0.;
    velino(0) = 1;
  }

  // get streamlength
  const double val = blitz::sum(blitz::abs(blitz::sum(velino(j)*derxy_(j,i),j)));
  const double strle = 2.0/val;

  // calculate tau
  // stabilization parameters for stationary case

  // compute tau_Mu
  const double re_tau_mu = mk * vel_norm * strle / (2.0 * visc);   /* convective : viscous forces */
  const double xi_tau_mu = DMAX(re_tau_mu, 1.0);
  tau_(0) = (DSQR(strle)*mk)/(4.0*visc*xi_tau_mu);

  // compute tau_Mp
  const double re_tau_mp = mk * vel_norm * hk / (2.0 * visc);      /* convective : viscous forces */
  const double xi_tau_mp = DMAX(re_tau_mp,1.0);
  tau_(1) = (DSQR(hk)*mk)/(4.0*visc*xi_tau_mp);

  // compute tau_C
  const double xi_tau_c = DMIN(re_tau_mp, 1.0);
  tau_(2) = 0.5*vel_norm*hk*xi_tau_c;

  return tau_;
}
}  // end namespace XFLUID

#endif

#endif
#endif
