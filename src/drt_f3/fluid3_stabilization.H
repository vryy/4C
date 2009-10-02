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

#include "../drt_lib/linalg_serialdensevector.H"


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
  template<DRT::Element::DiscretizationType DISTYPE, class M1, class M2>
  double HK(
      const M1&  evelnp,  ///< nodal velocity array
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

    // get element length for tau_Mp/tau_C: volume-equival. diameter/sqrt(3)
    const double hk = pow((6.*vol/PI),(1.0/3.0))/sqrt(3.0);

    return hk;
  }

  template <class V>
  void computeStabilization(
      const V&       shp_dx,
      const V&       shp_dy,
      const V&       shp_dz,
      const LINALG::Matrix<3,1>&      gpvelnp,
      const int   numparamvelx,
      const bool instationary,
      const double visc,
      const double hk,
      const double mk,
      const double timefac,
      double& tau_stab_M,
      double& tau_stab_Mp,
      double& tau_stab_C
      )
  {
    const int nsd = 3;

    // get velocity norm
    const double vel_norm = std::sqrt(gpvelnp(0)*gpvelnp(0)+gpvelnp(1)*gpvelnp(1)+gpvelnp(2)*gpvelnp(2));

    // normed velocity at Gauss point
    LINALG::Matrix<3,1> velino;
    if (vel_norm>=1e-6)
    {
        for (int isd = 0; isd < nsd; ++isd)
            velino(isd) = gpvelnp(isd)/vel_norm;
    }
    else
    {
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

    double val = 0.0;
    for (int iparam = 0; iparam < numparamvelx; ++iparam)
        val += fabs(velinoder(iparam));

    if (val == 0) // can happen, if ALL aprrox derivatives (shp_dx, shp_dy and shp_dz) are zero as in XFEM simulations
      val = 1.0;
    const double strle = 2.0/val;

    if (std::isinf(strle))
      dserror("Infinite stream-length detected! Quitting...");

    if (instationary)
    {
        // calculate tau: stabilization parameters for stationary case

        const double visceff = visc;
        /* viscous : reactive forces */
        const double re1 = 4.0 * timefac * visceff / (mk * DSQR(strle));

        /* convective : viscous forces */
        const double re2 = mk * vel_norm * strle / (2.0 * visceff);

        const double xi1 = max(re1,1.0);
        const double xi2 = max(re2,1.0);

        tau_stab_M = DSQR(strle) / (DSQR(strle)*xi1+( 4.0 * timefac*visceff/mk)*xi2);

        // compute tau_Mp
        //    stability parameter definition according to Franca and Valentin (2000)
        //                                       and Barrenechea and Valentin (2002)

         /* viscous : reactive forces */
        const double re_viscous = 4.0 * timefac * visceff / (mk * DSQR(hk));
        /* convective : viscous forces */
        const double re_convect = mk * vel_norm * hk / (2.0 * visceff);

        const double xi_viscous = max(re_viscous,1.0);
        const double xi_convect = max(re_convect,1.0);

        /*
                        xi1,xi2 ^
                                |      /
                                |     /
                                |    /
                              1 +---+
                                |
                                |
                                |
                                +--------------> re1,re2
                                    1
        */
        tau_stab_Mp = DSQR(hk) / (DSQR(hk) * xi_viscous + ( 4.0 * timefac * visceff/mk) * xi_convect);

        /*------------------------------------------------------ compute tau_C ---*/
        /*-- stability parameter definition according to Codina (2002), CMAME 191
         *
         * Analysis of a stabilized finite element approximation of the transient
         * convection-diffusion-reaction equation using orthogonal subscales.
         * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
         *
         * */
        //tau[2] = sqrt(DSQR(visc)+DSQR(0.5*vel_norm*hk));

        // Wall Diss. 99
        /*
                            xi2 ^
                                |
                              1 |   +-----------
                                |  /
                                | /
                                |/
                                +--------------> Re2
                                    1
        */
        const double xi_tau_c = min(re2,1.0);
        tau_stab_C = vel_norm * hk * 0.5 * xi_tau_c /timefac;
    }
    else
    {
        // calculate tau: stabilization parameters for stationary case

        // compute tau_Mu
        const double re_tau_mu = mk * vel_norm * strle / (2.0 * visc);   /* convective : viscous forces */
        const double xi_tau_mu = max(re_tau_mu, 1.0);
        tau_stab_M = (DSQR(strle)*mk)/(4.0*visc*xi_tau_mu);

        // compute tau_Mp
        const double re_tau_mp = mk * vel_norm * hk / (2.0 * visc);      /* convective : viscous forces */
        const double xi_tau_mp = max(re_tau_mp,1.0);
        tau_stab_Mp = (DSQR(hk)*mk)/(4.0*visc*xi_tau_mp);

        // compute tau_C
        const double xi_tau_c = min(re_tau_mp, 1.0);
        tau_stab_C = 0.5*vel_norm*hk*xi_tau_c;
    }
  }

}  // end namespace XFLUID

#endif

#endif
#endif
