/*----------------------------------------------------------------------*/
/*!
 * \file fluid_ele_calc_stabilization.H
 *
 * \brief stabilization related, templated functions
 *
 * <pre>
 * Maintainer: Florian Henke
 *             henke@lnm.mw.tum.de
 *             http://www.lnm.mw.tum.de
 *             089 - 289-15265
 * </pre>
 */
/*----------------------------------------------------------------------*/

#ifndef COMBUST3_STABILIZATION_H
#define COMBUST3_STABILIZATION_H

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_local_connectivity_matrices.H"
#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_lib/drt_element.H"
#include "../linalg/linalg_serialdensevector.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_xfem/enrichment_utils.H"


namespace COMBUST
{
namespace UTILS
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
        return 0.333333333333333333333; // 1/3
        break;
      case DRT::Element::hex20:
      case DRT::Element::hex27:
      case DRT::Element::tet10:
      case DRT::Element::wedge15:
        return 0.083333333333333333333; // 1/12
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
  template <class V>
  double Streamlength(
      const V&     shp_dx,
      const V&     shp_dy,
      const V&     shp_dz,
      const LINALG::Matrix<3,1>& gpvelnp,  /// velocity at Gaussian point
      const double               vel_norm, /// velocity norm
      const int                  numparamvelx
  )
  {
    // normed velocity at Gauss point
    LINALG::Matrix<3,1> velino;
    if (vel_norm >= 1e-6)
    {
      velino.Update(1.0/vel_norm,gpvelnp);
    }
    else
    {
      //cout << "/!\\ set x-velocity for streamlength computation" << endl;
      velino.Clear();
      velino(0) = 1.0;
    }

    // get streamlength
    //const double val = sum(abs(velino(j)*derxy(j,i)));
    LINALG::SerialDenseVector velinoder(numparamvelx,true);
    for (int iparam = 0; iparam < numparamvelx; ++iparam)
    {
      velinoder(iparam) += velino(0)*shp_dx(iparam);
      velinoder(iparam) += velino(1)*shp_dy(iparam);
      velinoder(iparam) += velino(2)*shp_dz(iparam);
    }
    // remark: "const double val = velinoder.Norm1();" ;
    double val = 0.0;
    for (int iparam = 0; iparam < numparamvelx; ++iparam)
      val += fabs(velinoder(iparam));

    if (val == 0.0) dserror("infinite stream length detected!");
    const double strle = 2.0/val;

    return strle;
  }

  //! calculate characteristic element size for stabilization
  // to be used for XFEM problems with normal enrichment strategy
  template <class V>
  double Streamlength(
      const V&                   shp,
      const LINALG::Matrix<3,1>& gpvelnp,  /// velocity at Gaussian point
      const double               vel_norm, /// velocity norm
      const int                  numparamvel
  )
  {
#ifdef DEBUG
  if (numparamvel!=8) dserror("This will only work for hex8 elements");
#endif
    // normed velocity at Gauss point
    LINALG::Matrix<3,1> velino(true);
    if (vel_norm >= 1e-6)
    {
      velino.Update(1.0/vel_norm,gpvelnp);
    }
    else
    {
      //cout << "/!\\ set x-velocity for streamlength computation" << endl;
      velino.Clear();
      velino(0) = 1.0;
    }

    // get streamlength
    //const double val = sum(abs(velino(j)*derxy(j,i)));
    // remark: extend vector from 8 to 16 to bring in normal shape functions correctly
    LINALG::SerialDenseVector velinoder(numparamvel+8,true);
    for (int iparam = 0; iparam < numparamvel; ++iparam)
    {
      velinoder(iparam) += velino(0)*shp.velx.dx.s(iparam);
      velinoder(iparam) += velino(1)*shp.vely.dy.s(iparam);
      velinoder(iparam) += velino(2)*shp.velz.dz.s(iparam);

      // store values belonging to normal shape functions
      velinoder(iparam+8) += velino(0)*shp.velx.dx.n(iparam);
      velinoder(iparam+8) += velino(1)*shp.vely.dy.n(iparam);
      velinoder(iparam+8) += velino(2)*shp.velz.dz.n(iparam);
    }
    // remark: "const double val = velinoder.Norm1();" ;
    double val = 0.0;
    for (int iparam = 0; iparam < numparamvel+8; ++iparam)
      val += fabs(velinoder(iparam));

    if (val == 0.0) dserror("infinite stream length detected!");
    const double strle = 2.0/val;

    return strle;
  }


  //! calculate volume based on numerical integration
  template<DRT::Element::DiscretizationType DISTYPE, class M2>
  double VolumeViaNumIntegration(
      const M2&  xyze     ///< element nodal position array
  )
  {
    const int numnode = DRT::UTILS::DisTypeToNumNodePerEle<DISTYPE>::numNodePerElement;

    // use one point integration rule to calculate hk at element center
    const DRT::UTILS::GaussRule3D integrationrule_stabili = getIntegrationRuleForStabilization<DISTYPE>();

    // integration points
    const DRT::UTILS::IntegrationPoints3D intpoints(integrationrule_stabili);

    // shape functions and derivs at element center
    LINALG::Matrix<3,1> e;
    e(0) = intpoints.qxg[0][0];
    e(1) = intpoints.qxg[0][1];
    e(2) = intpoints.qxg[0][2];
    const double wquad = intpoints.qwgt[0];

    LINALG::Matrix<3,numnode> deriv;
    DRT::UTILS::shape_function_3D_deriv1(deriv, e(0), e(1), e(2), DISTYPE);

    // get Jacobian matrix and determinant
    // xjm_ = deriv_(i,k)*xyze(j,k);
    LINALG::Matrix<3,3> xjm;
    xjm.MultiplyNT(deriv,xyze);

    const double vol = wquad * xjm.Determinant();

    return vol;
  }

  //! calculate characteristic element size for stabilization
  template<DRT::Element::DiscretizationType DISTYPE, class M1, class M2>
  double HK(
      const M1&  evelnp,  ///< nodal velocity array
      const M2&  xyze     ///< element nodal position array
  )
  {
    const double vol = VolumeViaNumIntegration<DISTYPE>(xyze);

    // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
    const double hk = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

    return hk;
  }


  void computeStabilizationParams(
      const LINALG::Matrix<3,1>& gpvelnp,  /// velocity at Gaussian point
      const LINALG::Matrix<3,3>& xji,      /// inverse of transposed Jacobian matrix
      const bool   instationary,
      const double dynvisc,                /// dyamic viscosity
      const double dens,                   /// density
      const double vel_norm,
      const double strle,
      const double hk,
      const double mk,
      const double timefac,
      const double dt,
      const enum INPAR::FLUID::TauType tautype,
      double& tau_stab_Mu,
      double& tau_stab_Mp,
      double& tau_stab_C
      );

} // namespace UTILS
} // namespace COMBUST

#endif
