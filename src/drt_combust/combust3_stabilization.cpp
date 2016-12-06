/*----------------------------------------------------------------------*/
/*!
\file combust3_stabilization.cpp

\brief compute stabilization paramters

\level 2

<pre>
\maintainer Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
*/
/*----------------------------------------------------------------------*/


#include "combust3_stabilization.H"

void COMBUST::UTILS::computeStabilizationParams(
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
  //---------------------------------------------------------------------
  // preliminary definition of values which will already be computed for
  // tau_M and later be used for tau_C again by some of the subsequent
  // stabilization parameter definitions
  //---------------------------------------------------------------------
  double traceG = 0.0;
  double Gnormu = 0.0;
  double Gvisc  = 0.0;

  double re12   = 0.0;
  double c3     = 0.0;

  //---------------------------------------------------------------------
  // first step: computation of tau_M with the following options
  // (both with or without inclusion of dt-part):
  // A) definition according to Taylor et al. (1998)
  //    -> see also Gravemeier and Wall (2010) for version for
  //       variable-density flow at low Mach number
  // B) combined definition according to Franca and Valentin (2000) as
  //    well as Barrenechea and Valentin (2002)
  //    -> differentiating tau_Mu and tau_Mp for this definition
  // C) definition according to Shakib (1989) / Shakib and Hughes (1991)
  //    -> differentiating tau_Mu and tau_Mp for this definition
  // D) definition according to Codina (1998)
  //    -> differentiating tau_Mu and tau_Mp for this definition
  //---------------------------------------------------------------------
  // computation depending on which parameter definition is used
  switch (tautype)
  {
  case INPAR::FLUID::tau_taylor_hughes_zarins:
  case INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt:
  case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen:
  case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt:
  case INPAR::FLUID::tau_taylor_hughes_zarins_scaled:
  case INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt:
  {
    /*

    literature:
    1) C.A. Taylor, T.J.R. Hughes, C.K. Zarins, Finite element modeling
       of blood flow in arteries, Comput. Methods Appl. Mech. Engrg. 158
       (1998) 155-196.
    2) V. Gravemeier, W.A. Wall, An algebraic variational multiscale-
       multigrid method for large-eddy simulation of turbulent variable-
       density flow at low Mach number, J. Comput. Phys. 229 (2010)
       6047-6070.
       -> version for variable-density low-Mach-number flow as implemented
          here, which corresponds to version for incompressible flow as
          given in the previous publications when density is constant

                                                                           1
                     +-                                               -+ - -
                     |        2                                        |   2
                     | c_1*rho                                  2      |
          tau  = C * | -------   +  c_2*rho*u*G*rho*u  +  c_3*mu *G:G  |
             M       |     2                                           |
                     |   dt                                            |
                     +-                                               -+

          with the constants and covariant metric tensor defined as follows:

          C   = 1.0 (not explicitly defined here),
          c_1 = 4.0 (for version with dt), 0.0 (for version without dt),
          c_2 = 1.0 (not explicitly defined here),
          c_3 = 12.0/m_k (36.0 for linear and 144.0 for quadratic elements)

                  +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+

                  +----
                   \
          G : G =   +   G   * G
                   /     ij    ij
                  +----
                   i,j
                             +----
                             \
          rho*u*G*rho*u  =   +   rho*u * G  *rho*u
                             /        i   ij      j
                            +----
                              i,j
    */

    // definition of constants as described above
    double c1 = 4.0;
    if ((tautype == INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt or
         tautype == INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt or
         tautype == INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt) or
         not instationary)
      c1 = 0.0;
    c3 = 12.0/mk;

    // computation of various values derived from covariant metric tensor
    // (trace of covariant metric tensor required for computation of tau_C below)
    double G;
    double normG = 0.0;
    const double dens_sqr = dens*dens;
    for (int nn=0;nn<3;++nn)
    {
      traceG += xji(nn,0)*xji(nn,0) + xji(nn,1)*xji(nn,1) + xji(nn,2)*xji(nn,2);
      for (int rr=0;rr<3;++rr)
      {
        G = xji(nn,0)*xji(rr,0) + xji(nn,1)*xji(rr,1) + xji(nn,2)*xji(rr,2);
        normG+=G*G;
        Gnormu+=dens_sqr*gpvelnp(nn,0)*G*gpvelnp(rr,0);
      }
    }

    // compute viscous part
    Gvisc = c3*dynvisc*dynvisc*normG;

    // computation of stabilization parameters tau_Mu and tau_Mp
    // -> identical for the present definitions
    tau_stab_Mu = 1.0/(sqrt(c1*dens_sqr/(dt*dt) + Gnormu + Gvisc));
    tau_stab_Mp = tau_stab_Mu;
  }
  break;

  case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall:
  case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt:
  {
    // relating convective to viscous part (re02: tau_Mu, re12: tau_Mp)
    const double re02 = mk * dens * vel_norm * strle / (2.0 * dynvisc);
                 re12 = mk * dens * vel_norm * hk / (2.0 * dynvisc);

    // respective "switching" parameters
    const double xi02 = std::max(re02,1.0);
    const double xi12 = std::max(re12,1.0);

    if (tautype == INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt or
        not instationary)
    {
      tau_stab_Mu = (DSQR(strle)*mk)/(4.0*dynvisc*xi02);
      tau_stab_Mp = (DSQR(hk)*mk)/(4.0*dynvisc*xi12);
    }
    else
    {
      // various parameter computations for case with dt:
      // relating viscous to reactive part (re01: tau_Mu, re11: tau_Mp)
      const double re01 = 4.0 * timefac * dynvisc / (mk * dens * DSQR(strle));
      const double re11 = 4.0 * timefac * dynvisc / (mk * dens * DSQR(hk));

      // respective "switching" parameters
      const double xi01 = std::max(re01,1.0);
      const double xi11 = std::max(re11,1.0);

      tau_stab_Mu = timefac*DSQR(strle)/(DSQR(strle)*dens*xi01+(4.0*timefac*dynvisc/mk)*xi02);
      tau_stab_Mp = timefac*DSQR(hk)/(DSQR(hk)*dens*xi11+(4.0*timefac*dynvisc/mk)*xi12);
    }
  }
  break;

  case INPAR::FLUID::tau_shakib_hughes_codina:
  case INPAR::FLUID::tau_shakib_hughes_codina_wo_dt:
  {
    /*

    literature on franca_barrenechea_valentin:
    1) L.P. Franca, F. Valentin, On an improved unusual stabilized
       finite element method for the advective-reactive-diffusive
       equation, Comput. Methods Appl. Mech. Engrg. 190 (2000) 1785-1800.
    2) G.R. Barrenechea, F. Valentin, An unusual stabilized finite
       element method for a generalized Stokes problem, Numer. Math.
       92 (2002) 652-677.


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


    literature on shakib_codina:
    1) F. Shakib, Finite element analysis of the compressible Euler and
       Navier-Stokes equations, PhD thesis, Division of Applied Mechanics,
       Stanford University, Stanford, CA, USA, 1989.
    2) F. Shakib, T.J.R. Hughes, A new finite element formulation for
       computational fluid dynamics: IX. Fourier analysis of space-time
       Galerkin/least-squares algorithms, Comput. Methods Appl. Mech.
       Engrg. 87 (1991) 35-58.
    3) R. Codina, Stabilized finite element approximation of transient
       incompressible flows using orthogonal subscales, Comput. Methods
       Appl. Mech. Engrg. 191 (2002) 4295-4321.

       constants defined as in Shakib (1989) / Shakib and Hughes (1991),
       merely slightly different with respect to c_3:

       c_1 = 4.0 (for version with dt), 0.0 (for version without dt),
       c_2 = 4.0,
       c_3 = 4.0/(m_k*m_k) (36.0 for linear, 576.0 for quadratic ele.)

       Codina (2002) proposed present version without dt and explicit
       definition of constants.

    */

    // definition of constants as described above
    double c1 = 4.0;
    if ((tautype == INPAR::FLUID::tau_shakib_hughes_codina_wo_dt) or
        not instationary)
      c1 = 0.0;
    const double c2 = 4.0;
    c3 = 4.0/(mk*mk);
    // alternative value as proposed in Shakib (1989): c3 = 16.0/(mk*mk);

    tau_stab_Mu = 1.0/(sqrt(c1*DSQR(dens)/DSQR(dt)
                          + c2*DSQR(dens)*DSQR(vel_norm)/DSQR(strle)
                          + c3*DSQR(dynvisc)/(DSQR(strle)*DSQR(strle))));
    tau_stab_Mp = 1.0/(sqrt(c1*DSQR(dens)/DSQR(dt)
                          + c2*DSQR(dens)*DSQR(vel_norm)/DSQR(hk)
                          + c3*DSQR(dynvisc)/(DSQR(hk)*DSQR(hk))));
  }
  break;

  case INPAR::FLUID::tau_codina:
  case INPAR::FLUID::tau_codina_wo_dt:
  {
    /*

      literature:
         R. Codina, Comparison of some finite element methods for solving
         the diffusion-convection-reaction equation, Comput. Methods
         Appl. Mech. Engrg. 156 (1998) 185-210.

         constants:
         c_1 = 1.0 (for version with dt), 0.0 (for version without dt),
         c_2 = 2.0,
         c_3 = 4.0/m_k (12.0 for linear, 48.0 for quadratic elements)

         Codina (1998) proposed present version without dt.

    */

    // definition of constants as described above
    double c1 = 1.0;
    if ((tautype == INPAR::FLUID::tau_shakib_hughes_codina_wo_dt) or
        not instationary)
      c1 = 0.0;
    const double c2 = 2.0;
    c3 = 4.0/mk;

    tau_stab_Mu = 1.0/(sqrt(c1*dens/dt
                          + c2*dens*vel_norm/strle
                          + c3*dynvisc/DSQR(strle)));
    tau_stab_Mp = 1.0/(sqrt(c1*dens/dt
                          + c2*dens*vel_norm/hk
                          + c3*dynvisc/DSQR(hk)));
  }
  break;

  default: dserror("unknown definition for tau_M\n %i  ", tautype);
  }  // end switch (tautype)


  //---------------------------------------------------------------------
  // second step: computation of tau_C with the following options:
  // A) definition according to Whiting (1999)/Whiting and Jansen (2001)
  // B) definition according to Taylor et al. (1998)
  // C) definition according to Wall (1999)
  // D) definition according to Codina (2002)
  //---------------------------------------------------------------------
  // computation depending on which parameter definition is used
  switch (tautype)
  {
  case INPAR::FLUID::tau_taylor_hughes_zarins:
  case INPAR::FLUID::tau_taylor_hughes_zarins_wo_dt:
  {
    /*

    literature:
       C.A. Taylor, T.J.R. Hughes, C.K. Zarins, Finite element modeling
       of blood flow in arteries, Comput. Methods Appl. Mech. Engrg. 158
       (1998) 155-196.

                                              1/2
                           (c_2*rho*u*G*rho*u)
                    tau  = -------------------
                       C       trace (G)


       -> see respective definitions for computation of tau_M above

    */

    tau_stab_C = sqrt(Gnormu)/traceG;
  }
  break;

  case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen:
  case INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen_wo_dt:
  {
    /*

    literature:
    1) C.H. Whiting, Stabilized finite element methods for fluid dynamics
       using a hierarchical basis, PhD thesis, Rensselaer Polytechnic
       Institute, Troy, NY, USA, 1999.
    2) C.H. Whiting, K.E. Jansen, A stabilized finite element method for
       the incompressible Navier-Stokes equations using a hierarchical
       basis, Int. J. Numer. Meth. Fluids 35 (2001) 93-116.

                                  1.0
                    tau  = ------------------
                       C    tau  * trace (G)
                               M

       -> see respective definitions for computation of tau_M above

    */

    tau_stab_C = 1.0/(tau_stab_Mu*traceG);
  }
  break;

  case INPAR::FLUID::tau_taylor_hughes_zarins_scaled:
  case INPAR::FLUID::tau_taylor_hughes_zarins_scaled_wo_dt:
  {
    /*

      Caution: This is an experimental version of a stabilization
               parameter definition which scales the definition
               for tau_C by Taylor et al. (1998) in a similar
               way as proposed below by Franca and Frey (1992)
               and Wall (1999) by appropriately defining an
               element Reynolds number based on the covariant
               metric tensor.

                  /                        1/2    \
                  |  /                    \       |                       1/2
                  | |  c_2*rho*u*G*rho*u  |       |    (c_2*rho*u*G*rho*u)
      tau  =  MIN | | ------------------- | | 1.0 | *  -------------------
         C        | |          2          |       |         trace (G)
                  | \    c_3*mu *G:G      /       |
                  \                               /
                    |                     |
                    -----------------------
                    element Reynolds number
                      based on covariant
                        metric tensor

       -> see respective definitions for computation of tau_M above

    */

    // element Reynolds number based on covariant metric tensor
    const double reG = sqrt(Gnormu/Gvisc);

    // "switching" parameter
    const double xi_tau_c = fmin(reG,1.0);

    tau_stab_C = xi_tau_c*sqrt(Gnormu)/traceG;
  }
  break;

  case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall:
  case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt:
  {
    /*

    literature:
    1) L.P. Franca, S.L. Frey, Stabilized finite element methods:
       II. The incompressible Navier-Stokes equations, Comput. Methods
       Appl. Mech. Engrg. 99 (1992) 209-293.
    2) W.A. Wall, Fluid-Struktur-Interaktion mit stabilisierten Finiten
       Elementen, Dissertation, Universitaet Stuttgart, 1999.

                 xi_tau_c ^
                          |
                        1 |   +-----------
                          |  /
                          | /
                          |/
                          +--------------> re12
                              1

       -> see respective definitions for computation of tau_M above

    */

    // "switching" parameter
    const double xi_tau_c = fmin(re12,1.0);

    tau_stab_C = 0.5 * dens * vel_norm * hk * xi_tau_c;
  }
  break;

  case INPAR::FLUID::tau_shakib_hughes_codina:
  case INPAR::FLUID::tau_shakib_hughes_codina_wo_dt:
  case INPAR::FLUID::tau_codina:
  case INPAR::FLUID::tau_codina_wo_dt:
  {
    /*

    literature:
       R. Codina, Stabilized finite element approximations of transient
       incompressible flows using orthogonal subscales, Comput. Methods
       Appl. Mech. Engrg. 191 (2002) 4295-4321.

       -> see respective definitions for computation of tau_M above

    */

    tau_stab_C = DSQR(hk)/(sqrt(c3)*tau_stab_Mp);
  }
  break;

  default: dserror("unknown definition for tau_C\n %i  ", tautype);
  }  // end switch (tautype)
}


void COMBUST::UTILS::computeStabilizationParamsEdgeBased(
    const double dynvisc,                /// dynamic viscosity
    const double dens,                   /// density
    const double vel_norm,
    const double normal_vel,             /// projection of velocity in normal direction of face
    const double hk,
    const enum INPAR::FLUID::EOS_TauType tautype,
    double& tau_conv,
    double& tau_div,
    double& tau_p,
    const bool instationary,
    const double timefac,
    const bool add_ghost_penalties,
    double& tau_ghost_first_order,
    double& tau_ghost_second_order,
    double& tau_pres_second_order
    )
{
  // non-dimensional factors
  double gamma_u   = 0.0; //scaling factor for streamline stabilization
  double gamma_div = 0.0; //scaling factor for divergence stabilization
  double gamma_p   = 0.0; //scaling factor for pressure stabilization

  //--------------------------------------------------------------------------------------------------------------
  //                       edge-oriented stabilization (EOS), continuous interior penalty (CIP)
  //--------------------------------------------------------------------------------------------------------------

  switch(tautype)
  {
    case INPAR::FLUID::EOS_tau_burman_fernandez_hansbo:
    case INPAR::FLUID::EOS_tau_burman_fernandez_hansbo_wo_dt:
    {
      // E.Burman, M.A.Fernandez and P.Hansbo 2006
      // "Continuous interior penalty method for Oseen's equations"
      //
      // velocity:
      //
      //             gamma_u * xi * h_K^2 * || u * n ||_0,inf_,K
      //
      //
      // divergence:
      //
      //             gamma_div * h_K^2 * || u ||_0,inf_,K * xi
      //
      // pressure:
      //
      //                         h_K^2
      //  gamma_p * xi * ---------------------
      //                   || u ||_0,inf_,K
      //
      //                    || u ||_0,inf_,K   *   h_K
      //  with    Re_K =  --------------------------------   and  xi = min(1,Re_K)
      //                                nu

      // values obtained by analyzing various test cases, i.e., Beltrami flow, lid-driven cavity flow, turbulent channel flow, ...
      // for further details see upcoming paper
      // -> for further details see upcoming paper
      // the following values worked best for beltrami flow, i.e., yield lowest errors
      // gamma_p = 3.0 / 100.0;
      // gamma_u = 5.0 / 100.0;
      // gamma_div = 5.0 / 10.0;
      // however, the divergence jump term comes along with numerical dissipation analogously to the
      // grad-div term of residual-based stabilizations
      // as this dissipation increases with increasing gamma_div, the previous choices are not appropriate
      // to compute turbulent flows
      // therefore, we suggest to use the same values as given in Burman 2007
      // note: we integrate each face once; Burman twice -> factor two for gamma compared to Burman
      gamma_p = 0.02;
      gamma_u = 0.02;
      gamma_div = 0.05 * gamma_u;

      // original values
      // gamma_p = 2.0 / 8.0;
      // gamma_u = 2.0 / 8.0;
      // gamma_div = 2.0 / 8.0;

      // element Reynold's number Re_K
      double Re_K = vel_norm*hk*dens/ (dynvisc * 6.0);
      double xi = std::min(1.0,Re_K);
      double xip = std::max(1.0,Re_K);

      //-----------------------------------------------
      // streamline
      tau_conv = dens * gamma_u * xi * hk*hk * normal_vel;

      //-----------------------------------------------
      // pressure
      //    if (max_vel_L2_norm > 1.0e-9)
      //      tau_p = gamma_p * xi/max_vel_L2_norm * hk*hk / density;
      //    else
      //      tau_p = gamma_p * hk * hk*hk/ kinvisc / density;
      // this does the same, but avoids division by velocity
      // note: this expression is closely related to the definition of
      //       tau_Mp according to Franca_Barrenechea_Valentin_Frey_Wall
      if (tautype == INPAR::FLUID::EOS_tau_burman_fernandez_hansbo_wo_dt or (not instationary))
        tau_p = gamma_p * hk * hk * hk / (dynvisc * xip);
      else
      {
        // respective "switching" parameter for transient / reactive part
        const double Re_R = 12.0 * timefac * dynvisc/ (dens * hk*hk);
        const double xir = std::max(1.0,Re_R);

        tau_p = gamma_p * hk * hk * hk / (dynvisc * xip + dens * hk * hk * xir / (timefac * 12.0));
      }

      //-----------------------------------------------
      // divergence
      tau_div= dens * gamma_div * xi * vel_norm * hk*hk;
    }
    break;
    case INPAR::FLUID::EOS_tau_burman_hansbo_dangelo_zunino:
    case INPAR::FLUID::EOS_tau_burman_hansbo_dangelo_zunino_wo_dt:
    {
      // this definition is derived form the following papers
      // E. Burman, P. Hansbo
      // "Edge stabilization for the generalized Stokes problem: A continuous interior penalty method"
      // Comput. Methods Appl. Mech. Engrg. 2006
      // C. D'Angelo, P. Zunino
      // "Numerical approximation with Nitsche's coupling of transient Stokes'/Darcy's flow problems applied to hemodynamics"
      // Applied Numerical Mathematics 2012

      // each face has to be evaluated only once
      gamma_p = 0.005;
      gamma_u = 0.05;
      gamma_div = 0.05 * gamma_u;

      //-----------------------------------------------
      // streamline
      tau_conv = dens * normal_vel * gamma_u * hk * hk; // Braack et al. 2006

      //-----------------------------------------------
      // divergence
      tau_div= dens * gamma_div * vel_norm * hk * hk; // Braack et al. 2006

      //-----------------------------------------------
      // pressure
      if ((tautype == INPAR::FLUID::EOS_tau_burman_hansbo_dangelo_zunino_wo_dt) or (not instationary))
        tau_p = gamma_p * hk * hk / (dynvisc / hk + dens * vel_norm / 6.0);
      else
        tau_p = gamma_p * hk * hk / (dens* hk / (timefac * 12.0) + dynvisc / hk + dens * vel_norm / 6.0);
    }
    break;
    case INPAR::FLUID::EOS_tau_burman:
    case INPAR::FLUID::EOS_tau_burman_fernandez:
    case INPAR::FLUID::EOS_tau_braack_burman_john_lube:
    case INPAR::FLUID::EOS_tau_Taylor_Hughes_Zarins_Whiting_Jansen_Codina_scaling:
    case INPAR::FLUID::EOS_tau_braack_burman_john_lube_wo_divjump:
    case INPAR::FLUID::EOS_tau_franca_barrenechea_valentin_wall:
    {
      dserror("edge-based tau-def not implemented yet");
    }
    break;
    default: dserror("unknown definition for tau\n %i  ", tautype); break;
  }


  //--------------------------------------------------------------------------------------------------------------
  //                                               ghost penalty
  //--------------------------------------------------------------------------------------------------------------

  if (add_ghost_penalties)
  {
    double gamma_ghost_penalty = 0.05;

    // visc ghost penalty + potential transient part
    tau_ghost_first_order = gamma_ghost_penalty*dynvisc*hk;
    if (instationary and
        ((tautype == INPAR::FLUID::EOS_tau_burman_hansbo_dangelo_zunino) or
          (tautype == INPAR::FLUID::EOS_tau_burman_fernandez_hansbo)))
      tau_ghost_first_order += (gamma_ghost_penalty * dens * hk*hk*hk/timefac);
                          // convective contribution contained in tau_conv
                          // + gamma_ghost_penalty * dens * normal_vel * hk*hk;

    tau_ghost_second_order = tau_ghost_first_order * hk*hk;
    tau_pres_second_order = tau_p * hk*hk;

    // set tau_conv
    // compared to, e.g, form of Burman, Fernandez and Hansbo, no xi-scaling is used here
    gamma_u = 0.05;
    tau_conv = gamma_u * dens * normal_vel * hk*hk;

  }
  return;
}

