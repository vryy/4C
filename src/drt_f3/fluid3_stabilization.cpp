/*----------------------------------------------------------------------*/
/*!
\file fluid3_stabilization.cpp

\brief compute stabilization paramters

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3_stabilization.H"

void FLD::UTILS::computeStabilizationParams(
    const LINALG::Matrix<3,1>& gpvelnp,  /// velocity at Gaussian point
    const LINALG::Matrix<3,3>& xji,      /// inverse of transposed Jacobian matrix
    const bool   instationary,
    const double dynvisc,                /// dynamic viscosity
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
    )
{
  // ---------------------------------------------------------------
  // computation of stabilization parameter tau
  // ---------------------------------------------------------------
  // compute stabilization parameters for instationary case
  if (instationary)
  {
    if (tautype == INPAR::FLUID::tautype_franca_barrenechea_valentin_wall)
    {
      /*----------------------------------------------------- compute tau_Mu ---*/
      /* stability parameter definition according to
       Barrenechea, G.R. and Valentin, F.: An unusual stabilized finite
       element method for a generalized Stokes problem. Numerische
       Mathematik, Vol. 92, pp. 652-677, 2002.
       http://www.lncc.br/~valentin/publication.htm

       and:

      Franca, L.P. and Valentin, F.: On an Improved Unusual Stabilized
      Finite Element Method for the Advective-Reactive-Diffusive
      Equation. Computer Methods in Applied Mechanics and Enginnering,
      Vol. 190, pp. 1785-1800, 2000.
      http://www.lncc.br/~valentin/publication.htm                   */

      /* viscous : reactive forces */
      const double re01 = 4.0 * timefac * dynvisc / (mk * dens * DSQR(strle));

      /* convective : viscous forces */
      const double re02 = mk * dens * vel_norm * strle / (2.0 * dynvisc);

      const double xi01 = DMAX(re01,1.0);
      const double xi02 = DMAX(re02,1.0);

      tau_stab_Mu = timefac*DSQR(strle) / (DSQR(strle)*dens*xi01 + (4.0*timefac*dynvisc/mk)*xi02);

      // compute tau_Mp
      //    stability parameter definition according to Franca and Valentin (2000)
      //                                       and Barrenechea and Valentin (2002)

       /* viscous : reactive forces */
      const double re11 = 4.0 * timefac * dynvisc / (mk * dens * DSQR(hk));
      /* convective : viscous forces */
      const double re12 = mk * dens * vel_norm * hk / (2.0 * dynvisc);

      const double xi11 = DMAX(re11,1.0);
      const double xi12 = DMAX(re12,1.0);
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
      tau_stab_Mp = timefac*DSQR(hk) / (DSQR(hk)*dens*xi11+(4.0 * timefac * dynvisc/mk) * xi12);
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
      const double xi_tau_c = DMIN(re02,1.0);
      tau_stab_C = dens * vel_norm * hk * 0.5 * xi_tau_c;
    }
    else if(tautype == INPAR::FLUID::tautype_bazilevs)
    {
      /* INSTATIONARY FLOW PROBLEM, ONE-STEP-THETA, BDF2

      tau_M: Bazilevs et al.
                                                                 1.0
                   +-                                       -+ - ---
                   |                                         |   2.0
                   | 4.0    n+1       n+1          2         |
            tau  = | --- + u     * G u     + C * nu  * G : G |
               M   |   2           -          I        -   - |
                   | dt            -                   -   - |
                   +-                                       -+
      tau_C: Bazilevs et al., derived from the fine scale complement Shur
            operator of the pressure equation
                                      1.0
                      tau  = -----------------
                         C            /     \
                              tau  * | g * g |
                                 M    \-   -/
      */
      /*            +-           -+   +-           -+   +-           -+
                    |             |   |             |   |             |
                    |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
              G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
               ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                    |    i     j  |   |    i     j  |   |    i     j  |
                    +-           -+   +-           -+   +-           -+
      */
      /*            +----
                     \
            G : G =   +   G   * G
            -   -    /     ij    ij
            -   -   +----
                     i,j
      */
      /*                      +----
             n+1       n+1     \     n+1          n+1
            u     * G u     =   +   u    * G   * u
                    -          /     i     -ij    j
                    -         +----        -
                               i,j
      */
      double G;
      double normG = 0.0;
      double Gnormu = 0.0;
      const double dens_sqr = dens*dens;
      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          G = xji(nn,0)*xji(rr,0) + xji(nn,1)*xji(rr,1) + xji(nn,2)*xji(rr,2);
          normG+=G*G;
          Gnormu+=dens_sqr*gpvelnp(nn,0)*G*gpvelnp(rr,0);
        }
      }

      // definition of constant
      // (Akkerman et al. (2008) used 36.0 for quadratics, but Stefan
      //  brought 144.0 from Austin...)
      const double CI = 12.0/mk;

      /*                                                         1.0
               +-                                       -+ - ---
               |                                         |   2.0
               | 4.0    n+1       n+1          2         |
        tau  = | --- + u     * G u     + C * nu  * G : G |
           M   |   2           -          I        -   - |
               | dt            -                   -   - |
               +-                                       -+
       */
      tau_stab_Mu = 1.0/(sqrt((4.0*dens_sqr)/(dt*dt)+Gnormu+CI*dynvisc*dynvisc*normG));
      tau_stab_Mp = tau_stab_Mu;
       /*           +-     -+   +-     -+   +-     -+
                   |       |   |       |   |       |
                   |  dr   |   |  ds   |   |  dt   |
              g  = |  ---  | + |  ---  | + |  ---  |
               i   |  dx   |   |  dx   |   |  dx   |
                   |    i  |   |    i  |   |    i  |
                   +-     -+   +-     -+   +-     -+
      */
      /*           +----
                    \
           g * g =   +   g * g
           -   -    /     i   i
                   +----
                     i
      */
      double g;
      double normgsq = 0.0;
      for (int rr=0;rr<3;++rr)
      {
        g = xji(rr,0) + xji(rr,1) + xji(rr,2);
        normgsq += g*g;
      }

      /*
                              1.0
                tau  = -----------------
                   C            /     \
                        tau  * | g * g |
                           M    \-   -/
      */
      tau_stab_C = 1.0/(tau_stab_Mu*normgsq);
    }
    else dserror("unknown definition of tau\n");
  }
  // compute stabilization parameters for stationary case
  else
  {
      // compute tau_Mu
      const double re_tau_mu = mk * dens * vel_norm * strle / (2.0 * dynvisc);   /* convective : viscous forces */
      const double xi_tau_mu = DMAX(re_tau_mu, 1.0);
      tau_stab_Mu = (DSQR(strle)*mk)/(4.0*dynvisc*xi_tau_mu);
       // compute tau_Mp
      const double re_tau_mp = mk * dens * vel_norm * hk / (2.0 * dynvisc);      /* convective : viscous forces */
      const double xi_tau_mp = DMAX(re_tau_mp,1.0);
      tau_stab_Mp = (DSQR(hk)*mk)/(4.0*dynvisc*xi_tau_mp);
       // compute tau_C
      const double xi_tau_c = min(re_tau_mp, 1.0);
      tau_stab_C = dens*0.5*vel_norm*hk*xi_tau_c;
  }
}


#endif
#endif
