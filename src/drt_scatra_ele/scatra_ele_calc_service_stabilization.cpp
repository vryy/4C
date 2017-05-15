/*----------------------------------------------------------------------*/
/*!
\file scatra_ele_calc_service_stabilization.cpp

\brief Internal implementation of ScaTra element

\level 1

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>
*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_calc.H"

#include "scatra_ele_parameter_std.H"
#include "scatra_ele_parameter_timint.H"

#include "scatra_ele_calc_utils.H"

#include "../drt_lib/standardtypes_cpp.H"  // for EPS13 and so on
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"


/*----------------------------------------------------------------------*
 |  calculate stabilization parameter  (private)              gjb 06/08 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcTau(
  double&                         tau,         //!< the stabilisation parameters (one per transported scalar)
  const double                    diffus,      //!< diffusivity or viscosity
  const double                    reacoeff,    //!< reaction coefficient
  const double                    densnp,      //!< density at t_(n+1)
  const LINALG::Matrix<nsd_,1>&   convelint,   //!< convective velocity at integration point
  const double                    vol          //!< element volume
  )
{
  //----------------------------------------------------------------------
  // computation of stabilization parameters depending on definition used
  //----------------------------------------------------------------------
  switch (scatrapara_->TauDef())
  {
  case INPAR::SCATRA::tau_taylor_hughes_zarins:
  case INPAR::SCATRA::tau_taylor_hughes_zarins_wo_dt:
  {
    CalcTauTaylorHughesZarins(tau,diffus,reacoeff,densnp,convelint);
    break;
  }
  case INPAR::SCATRA::tau_franca_valentin:
  case INPAR::SCATRA::tau_franca_valentin_wo_dt:
  {
    CalcTauFrancaValentin(tau,diffus,reacoeff,densnp,convelint,vol);
    break;
  }
  case INPAR::SCATRA::tau_shakib_hughes_codina:
  case INPAR::SCATRA::tau_shakib_hughes_codina_wo_dt:
  {
    CalcTauFrancaShakibCodina(tau,diffus,reacoeff,densnp,convelint,vol);
    break;
  }
  case INPAR::SCATRA::tau_codina:
  case INPAR::SCATRA::tau_codina_wo_dt:
  {
    CalcTauCodina(tau,diffus,reacoeff,densnp,convelint,vol);
    break;
  }
  case INPAR::SCATRA::tau_franca_madureira_valentin:
  case INPAR::SCATRA::tau_franca_madureira_valentin_wo_dt:
  {
    CalcTauFrancaMadureiraValentin(tau,diffus,reacoeff,densnp,vol);
    break;
  }
  case INPAR::SCATRA::tau_exact_1d:
  {
    CalcTau1DExact(tau,diffus,reacoeff,densnp,convelint,vol);
    break;
  }
  case INPAR::SCATRA::tau_zero:
  {
    // set tau's to zero (-> no stabilization effect)
    tau = 0.0;
    break;
  }
  default:
  {
    dserror("unknown definition for stabilization parameter tau\n");
    break;
  }
  }

  return;
} //ScaTraEleCalc::CalcTau

/*----------------------------------------------------------------------*
 |  calculation of tau according to Taylor, Hughes and Zarins  vg 01/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcTauTaylorHughesZarins(
  double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
  const double                 diffus,    //!< diffusivity or viscosity
  const double                 reacoeff,  //!< reaction coefficient
  const double                 densnp,    //!< density at t_(n+1)
  const LINALG::Matrix<nsd_,1>& convelint  //!< convective velocity at integration point
  )
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
  -> version for variable-density scalar transport equation as
  implemented here, which corresponds to constant-density
  version as given in the previous publication when density
  is constant

                                                                  1
            +-                                               -+ - -
            |        2                                        |   2
            | c_1*rho                                  2      |
  tau = C * | -------   +  c_2*rho*u*G*rho*u  +  c_3*mu *G:G  |
            |     2                                           |
            |   dt                                            |
            +-                                               -+

  with the constants and covariant metric tensor defined as follows:

  C   = 1.0 (not explicitly defined here),
  c_1 = 4.0,
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

  // get element-type constant for tau
  const double mk = SCATRA::MK<distype>();

  // effective velocity at element center:
  // (weighted) convective velocity + individual migration velocity
  LINALG::Matrix<nsd_,1> veleff(convelint);

  // total reaction coefficient sigma_tot: sum of "artificial" reaction
  // due to time factor and reaction coefficient (reaction coefficient
  // ensured to be zero in GetMaterialParams for non-reactive material)
  double sigma_tot = reacoeff;
  if (scatrapara_->TauDef() == INPAR::SCATRA::tau_taylor_hughes_zarins) sigma_tot += 1.0/scatraparatimint_->Dt();

  // computation of various values derived from covariant metric tensor
  double G;
  double normG(0.0);
  double Gnormu(0.0);
  const double dens_sqr = densnp*densnp;
  for (unsigned nn=0;nn<nsd_;++nn)
  {
    for (unsigned rr=0;rr<nsd_;++rr)
    {
      G = xij_(nn,0)*xij_(rr,0);
      for (unsigned tt=1;tt<nsd_;tt++)
      {
        G += xij_(nn,tt)*xij_(rr,tt);
      }
      normG+=G*G;
      Gnormu+=dens_sqr*veleff(nn,0)*G*veleff(rr,0);
    }
  }

  // definition of constants as described above
  const double c1 = 4.0;
  const double c3 = 12.0/mk;

  // compute diffusive part
  const double Gdiff = c3*diffus*diffus*normG;

  // computation of stabilization parameter tau
  tau = 1.0/(sqrt(c1*dens_sqr*DSQR(sigma_tot) + Gnormu + Gdiff));

  return;
}


/*----------------------------------------------------------------------*
 |  calculation of tau according to Franca and Valentin        vg 01/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcTauFrancaValentin(
  double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
  const double                 diffus,    //!< diffusivity or viscosity
  const double                 reacoeff,  //!< reaction coefficient
  const double                 densnp,    //!< density at t_(n+1)
  const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
  const double                 vol        //!< element volume
  )
{
  /*

  literature:
  L.P. Franca, F. Valentin, On an improved unusual stabilized
  finite element method for the advective-reactive-diffusive
  equation, Comput. Methods Appl. Mech. Engrg. 190 (2000) 1785-1800.


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

  // get element-type constant for tau
  const double mk = SCATRA::MK<distype>();

  // get Euclidean norm of (weighted) velocity at element center
  double vel_norm;
  vel_norm = convelint.Norm2();

  // total reaction coefficient sigma_tot: sum of "artificial" reaction
  // due to time factor and reaction coefficient (reaction coefficient
  // ensured to be zero in GetMaterialParams for non-reactive material)
  double sigma_tot = reacoeff;
  if (scatrapara_->TauDef() == INPAR::SCATRA::tau_franca_valentin)
    sigma_tot += 1.0/scatraparatimint_->TimeFac();

  // calculate characteristic element length
  const double h = CalcCharEleLength(vol,vel_norm, convelint);

  // various parameter computations:
  // relating convective to viscous part
  const double epe = mk * densnp * vel_norm * h;
  // relating viscous to reactive part
  double epe1 = 0.0;
  if (scatrapara_->TauDef() == INPAR::SCATRA::tau_franca_valentin or reacoeff != 0.0)
    epe1 = 2.0*diffus/(mk*densnp*sigma_tot*DSQR(h));

  // respective "switching" parameters
  const double xi  = std::max(epe,1.0*diffus);
  const double xi1 = std::max(epe1,1.0);

  tau = DSQR(h)/(DSQR(h)*densnp*sigma_tot*xi1 + 2.0*xi/mk);

  return;
}

/*----------------------------------------------------------------------*
 |  calculation of tau according to Franca, Shakib and Codina  vg 01/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcTauFrancaShakibCodina(
  double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
  const double                 diffus,    //!< diffusivity or viscosity
  const double                 reacoeff,  //!< reaction coefficient
  const double                 densnp,    //!< density at t_(n+1)
  const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
  const double                 vol        //!< element volume
  )
{
  /*
  literature:
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

  All those proposed definitions were for non-reactive incompressible
  flow; they are adapted to potentially reactive scalar transport
  equations with potential density variations here.

  constants defined as in Shakib (1989) / Shakib and Hughes (1991),
  merely slightly different with respect to c_3:

  c_1 = 4.0,
  c_2 = 4.0,
  c_3 = 4.0/(m_k*m_k) (36.0 for linear, 576.0 for quadratic ele.)

  Codina (2002) proposed present version without dt and explicit
  definition of constants
  (condition for constants as defined here: c_2 <= sqrt(c_3)).

  */

  // get element-type constant for tau
  const double mk = SCATRA::MK<distype>();

  // get Euclidean norm of velocity
  const double vel_norm = convelint.Norm2();

  // total reaction coefficient sigma_tot: sum of "artificial" reaction
  // due to time factor and reaction coefficient (reaction coefficient
  // ensured to be zero in GetMaterialParams for non-reactive material)
  double sigma_tot = reacoeff;
  if (scatrapara_->TauDef() == INPAR::SCATRA::tau_shakib_hughes_codina) sigma_tot += 1.0/scatraparatimint_->Dt();

  // calculate characteristic element length
  const double h = CalcCharEleLength(vol,vel_norm,convelint);

  // definition of constants as described above
  const double c1 = 4.0;
  const double c2 = 4.0;
  const double c3 = 4.0/(mk*mk);
  // alternative value as proposed in Shakib (1989): c3 = 16.0/(mk*mk);

  tau = 1.0/(sqrt(c1*DSQR(densnp)*DSQR(sigma_tot)
                      + c2*DSQR(densnp)*DSQR(vel_norm)/DSQR(h)
                      + c3*DSQR(diffus)/(DSQR(h)*DSQR(h))));

  return;
}


/*----------------------------------------------------------------------*
 |  calculation of tau according to Codina                     vg 01/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcTauCodina(
  double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
  const double                 diffus,    //!< diffusivity or viscosity
  const double                 reacoeff,  //!< reaction coefficient
  const double                 densnp,    //!< density at t_(n+1)
  const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
  const double                 vol        //!< element volume
  )
{
  /*
  literature:
  R. Codina, Comparison of some finite element methods for solving
  the diffusion-convection-reaction equation, Comput. Methods
  Appl. Mech. Engrg. 156 (1998) 185-210.

  constants:
  c_1 = 1.0,
  c_2 = 2.0,
  c_3 = 4.0/m_k (12.0 for linear, 48.0 for quadratic elements)

  Codina (1998) proposed present version without dt.

  */

  // get element-type constant for tau
  const double mk = SCATRA::MK<distype>();

  // get Euclidean norm of velocity
  const double vel_norm = convelint.Norm2();

  // total reaction coefficient sigma_tot: sum of "artificial" reaction
  // due to time factor and reaction coefficient (reaction coefficient
  // ensured to be zero in GetMaterialParams for non-reactive material)
  double sigma_tot = reacoeff;
  if (scatrapara_->TauDef() == INPAR::SCATRA::tau_codina) sigma_tot += 1.0/scatraparatimint_->Dt();

  // calculate characteristic element length
  const double h = CalcCharEleLength(vol,vel_norm,convelint);

  // definition of constants as described above
  const double c1 = 1.0;
  const double c2 = 2.0;
  const double c3 = 4.0/mk;

  tau = 1.0/(c1*densnp*sigma_tot
                 + c2*densnp*vel_norm/h
                 + c3*diffus/(h*h));

  return;
}

/*---------------------------------------------------------------------------*
 |  calculation of tau according to Franca, Madureira and Valentin  vg 01/11 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcTauFrancaMadureiraValentin(
  double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
  const double                 diffus,    //!< diffusivity or viscosity
  const double                 reacoeff,  //!< reaction coefficient
  const double                 densnp,    //!< density at t_(n+1)
  const double                 vol        //!< element volume
  )
{
  /*
  This stabilization parameter is only intended to be used for
  reactive-diffusive problems such as structure-based scalar
  transport problems in case of potentially dominating reaction.

  literature:
  L.P. Franca, A.L. Madureira, F. Valentin, Towards multiscale
  functions: enriching finite element spaces with local but not
  bubble-like functions, Comput. Methods Appl. Mech. Engrg. 194
  (2005) 3006-3021.

  */

  // get element-type constant for tau
  const double mk = SCATRA::MK<distype>();

  // total reaction coefficient sigma_tot: sum of "artificial" reaction
  // due to time factor and reaction coefficient (reaction coefficient
  // ensured to be zero in GetMaterialParams for non-reactive material)
  double sigma_tot = reacoeff;
  if (scatrapara_->TauDef() == INPAR::SCATRA::tau_franca_madureira_valentin)
    sigma_tot += 1.0/scatraparatimint_->TimeFac();

  // calculate characteristic element length
  // -> currently: cubic/square root of element volume/area or
  //    element length (3-/2-/1-D)
  // cast dimension to a double variable -> pow()
  const double dim = (double) nsd_;
  const double h = std::pow(vol,1/dim);

  // parameter relating reactive to diffusive part
  double epe = 0.0;
  if (scatrapara_->TauDef() == INPAR::SCATRA::tau_franca_madureira_valentin or reacoeff != 0.0)
    epe = 2.0*diffus/(mk*densnp*sigma_tot*DSQR(h));

  // respective "switching" parameter
  const double xi = std::max(epe,1.0);

  // constant c_u as suggested in Badia and Codina (2010), method A
  // is set to be 1.0 here as in Franca et al. (2005)
  // alternative: 4.0 as suggested in Badia and Codina (2010) for
  // Darcy flow
  const double c_u = 1.0;

  if (scatrapara_->TauDef() == INPAR::SCATRA::tau_franca_madureira_valentin or reacoeff != 0.0)
    tau = DSQR(h)/(c_u*DSQR(h)*densnp*sigma_tot*xi + (2.0*diffus/mk));

  return;
}


/*----------------------------------------------------------------------*
 |  exact calculation of tau for 1D                            vg 01/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcTau1DExact(
  double&                      tau,       //!< the stabilisation parameters (one per transported scalar)
  const double                 diffus,    //!< diffusivity or viscosity
  const double                 reacoeff,  //!< reaction coefficient
  const double                 densnp,    //!< density at t_(n+1)
  const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
  const double                 vol        //!< element volume
  )
{
  // get number of dimensions (convert from int to double)
  const double dim = (double) nsd_;

  // get characteristic element length
  double h = std::pow(vol,(1.0/dim)); // equals streamlength in 1D

  // get Euclidean norm of (weighted) velocity at element center
  double vel_norm(0.0);
  vel_norm = convelint.Norm2();

  if (diffus < EPS14) dserror("Invalid diffusion coefficent");
  double epe = 0.5 * densnp * vel_norm * h / diffus;

  const double pp = exp(epe);
  const double pm = exp(-epe);
  double xi = 0.0;
  if (epe >= 700.0)
    tau = 0.5*h/vel_norm;
  else if (epe < 700.0 and epe > EPS15)
  {
    xi = (((pp+pm)/(pp-pm))-(1.0/epe)); // xi = coth(epe) - 1/epe
    // compute optimal stabilization parameter
    tau = 0.5*h*xi/vel_norm;

  }
  else tau = 0.0;

  return;
}


/*----------------------------------------------------------------------*
 |  calculation of characteristic element length               vg 01/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
double DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcCharEleLength(
  const double                  vol,        //!< element volume
  const double                  vel_norm,   //!< norm of velocity
  const LINALG::Matrix<nsd_,1>& convelint   //!< convective velocity at integration point
  )
{
  // define and initialize streamlength
  double h = 0.0;

  //---------------------------------------------------------------------
  // select from various definitions for characteristic element length
  //---------------------------------------------------------------------
  switch (scatrapara_->CharEleLength())
  {
    // a) streamlength due to Tezduyar et al. (1992) -> default
    // normed velocity vector
    case INPAR::SCATRA::streamlength:
    {
      LINALG::Matrix<nsd_,1> velino(true);
      if (vel_norm>=1e-6) velino.Update(1.0/vel_norm,convelint);
      else
      {
        velino.Clear();
        velino(0,0) = 1.0;
      }

      // get streamlength using the normed velocity at element centre
      LINALG::Matrix<nen_,1> tmp;
      tmp.MultiplyTN(derxy_,velino);
      const double val = tmp.Norm1();
      h = 2.0/val; // h=streamlength
    }
    break;

    // b) volume-equivalent diameter (warning: 3-D formula!)
    case INPAR::SCATRA::volume_equivalent_diameter:
    {
      h = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);
    }
    break;

    // c) cubic/square root of element volume/area or element length (3-/2-/1-D)
    case INPAR::SCATRA::root_of_volume:
    {
      // cast dimension to a double varibale -> pow()
      const double dim = double (nsd_ele_);
      h = std::pow(vol,1.0/dim);
    }
    break;

    default: dserror("unknown characteristic element length\n");
    break;
  } //switch (charelelength_)

  return h;
}


/*----------------------------------------------------------------------*
 |  calculate artificial diffusivity                           vg 10/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcArtificialDiff(
  const double                  vol,        //!< element volume
  const int                     k,          //!< id of current scalar
  const double                  densnp,     //!< density at t_(n+1)
  const LINALG::Matrix<nsd_,1>& convelint,  //!< convective velocity at integration point
  const LINALG::Matrix<nsd_,1>& gradphi,    //!< scalar gradient
  const double                  conv_phi,   //!< convective contribution
  const double                  scatrares,  //!< residual of convection-diffusion-reaction eq
  const double                  tau         //!< the stabilisation parameter
  )
{
  // get number of dimensions
  const double dim = (double) nsd_;

  // get characteristic element length as cubic root of element volume
  // (2D: square root of element area, 1D: element length)
  const double h = std::pow(vol,(1.0/dim));
//  const double h = CalcCharEleLength(vol,convelint.Norm2(),convelint); //std::pow(vol,(1.0/dim));

  // artificial diffusivity
  double artdiff = 0.0;

  // classical linear artificial all-scale subgrid diffusivity
  if (scatrapara_->ASSGDType() == INPAR::SCATRA::assgd_artificial)
  {
    // get element-type constant
    const double mk = SCATRA::MK<distype>();

    // velocity norm
    const double vel_norm = convelint.Norm2();

    // parameter relating convective and diffusive forces + respective switch
    double epe = 0.0;
    double xi = 1.0;
    if (diffmanager_->GetIsotropicDiff(k)>1.0e-8)
    {
      epe = 0.5 * mk * densnp * vel_norm * h / diffmanager_->GetIsotropicDiff(k);
      xi = std::min(epe,1.0);
    }

    // compute subgrid diffusivity
    artdiff = xi*0.5*densnp*vel_norm*h;
  }
  else if (scatrapara_->ASSGDType() == INPAR::SCATRA::assgd_lin_reinit)
  {
    artdiff = 0.005*h;
  }
  else if (scatrapara_->ASSGDType() == INPAR::SCATRA::assgd_codina)
  {
    double alpha = std::max(0.0,(0.7-2.0*diffmanager_->GetIsotropicDiff(k)/convelint.Norm2()/h));

    // gradient norm
    const double grad_norm = gradphi.Norm2();
    if (grad_norm > EPS8)
      artdiff = 0.5 * alpha * h * std::abs(scatrares) / grad_norm;
  }
  else if (scatrapara_->ASSGDType() == INPAR::SCATRA::assgd_yzbeta)
  {
   // phiref is the tuning parameter for this form of artificial diffusion
    const double phiref = 0.01;

    // gradient norm
    const double grad_norm = gradphi.Norm2();

    if (phiref > EPS12 and grad_norm > EPS12)
    {
      // normalized gradient of phi
      LINALG::Matrix<nsd_,1> normalized_gradphi(true);
      normalized_gradphi.Update(1.0,gradphi,0.0);
      normalized_gradphi.Scale(1.0/grad_norm);

      // compute reference length
      double h_sum(0.0);
      for (unsigned inode=0; inode<nen_; inode++)
      {
        double val(0.0);
        for (unsigned idim=0; idim<nsd_; idim++)
          val += normalized_gradphi(idim,0)*derxy_(idim,inode);

        h_sum += std::abs(val);
      }
      double h_dc(0.0);
      if (h_sum > EPS12)
        h_dc = 2.0/h_sum;
      else
        h_dc = h;

      // compute intermediate quantity
      double kappa_inter(0.0);
      for (unsigned idim=0; idim<nsd_; idim++)
      {
        double val = gradphi(idim,0)/phiref;
        kappa_inter += (val*val);
      }
//      if (kappa_inter < EPS10) dserror("Too low value");

      // smoothness parameter beta: 1 (smoother layers) or 2 (sharper layers)
      // note for 1.0, this form is equivalent to the Codina form above, except
      // for a different definition of the reference length
      const double beta = 2.0;

      // finally compute artificial diffusion
      artdiff = std::abs(scatrares/phiref) * pow(kappa_inter,beta/2.0-1.0) * pow(h_dc/2.0,beta);
    }
    else
      artdiff = 0.0;
  }
  else
  {
    // gradient norm
    const double grad_norm = gradphi.Norm2();

    if (grad_norm > EPS10)
    {
      // for the present definitions, sigma and a specific term (either
      // residual or convective term) are different
      double sigma = 0.0;
      double specific_term = 0.0;
      switch (scatrapara_->ASSGDType())
      {
        case INPAR::SCATRA::assgd_hughes:
        {
          if (eid_ == 0)
          {
            std::cout << "WARNING: Nonlinear isotropic artificial diffusion according to Hughes et al. (1986)\n";
            std::cout << "         is implemented based on the exact tau for 1D stationary problems!" << std::endl;
          }
          // remark on this warning:
          // 1. Here, tau is calculated based on the exact formula for 1-D stationary problems. This is inconsistent
          //    if the problem is not 1-D and/or other definitions than the ecaxt tau are chosen in the input file.
          //    Consistently, one has to use the same definition here.
          // 2. Instead of using sigma = tau_bhbar, Hughes et al. suggested to use sigma = tau_bhbar - tau to not double
          //    the SUPG stabilization. This is another inconsitent aspect on this implementation. To have the right tau
          //    here (i.e, not the one the last gauss point or even last step), one has to calaculate tau first. Then,
          //    sigma and, hence, the addition diffusion is computed based on this tau. Next, tau is recomputed with the
          //    diffusivity repaced by the original (physical) diffusivity plus the estimated artificial diffusivity. This
          //    is pobably not a good choice, because, first, tau is considered in the estimation of the artificial diffusion
          //    and then this artificial diffusion is incorporated into tau. This would reduce the effect. Perhaps, one
          //    should either consider tau in sigma or the artificial diffusion in tau. When changing this aspect, be aware
          //    that tau has to be computed after the subgrid-scale velocity has been calcuated since, for this calculation
          //    tau is overwritten by its value in the fluid field. Note that similar considerations may also hold for
          //    the methods by do Carmo and Almeida.

          // get norm of velocity vector b_h^par
          const double vel_norm_bhpar = abs(conv_phi/grad_norm);

          // compute stabilization parameter based on b_h^par
          // (so far, only exact formula for stationary 1-D implemented)
          // element Peclet number relating convective and diffusive forces
          double epe = 0.5 * vel_norm_bhpar * h / diffmanager_->GetIsotropicDiff(k);
          const double pp = exp(epe);
          const double pm = exp(-epe);
          double xi = 0.0;
          double tau_bhpar = 0.0;
          if (epe >= 700.0) tau_bhpar = 0.5*h/vel_norm_bhpar;
          else if (epe < 700.0 and epe > EPS15)
          {
            xi = (((pp+pm)/(pp-pm))-(1.0/epe)); // xi = coth(epe) - 1/epe
            // compute optimal stabilization parameter
            tau_bhpar = 0.5*h*xi/vel_norm_bhpar;
          }

          // compute sigma
          sigma = std::max(0.0,tau_bhpar-tau);

          // set specific term to convective term
          specific_term = conv_phi;
        }
        break;
        case INPAR::SCATRA::assgd_tezduyar:
        case INPAR::SCATRA::assgd_tezduyar_wo_phizero:
        {
          // velocity norm
          const double vel_norm = convelint.Norm2();

          // calculate stream length
          // according to John and Knobloch stream length in direction of b_h^par should be used
          //const double h_stream = CalcCharEleLength(vol,vel_norm);

          // get norm of velocity vector b_h^par
          const double vel_norm_bhpar = abs(conv_phi/grad_norm);

          // compute stabilization parameter based on b_h^par
          // (so far, only exact formula for stationary 1-D implemented)

          // compute sigma (version 1 according to John and Knobloch (2007))
          if (scatrapara_->ASSGDType() == INPAR::SCATRA::assgd_tezduyar_wo_phizero)
          {
            if (vel_norm > EPS10)
              sigma = (h/vel_norm)*(1.0-(vel_norm_bhpar/vel_norm));
          }
          else
          {
            // compute sigma (version 2 according to John and Knobloch (2007))
            // setting scaling phi_0=1.0 as in John and Knobloch (2007)
            const double phi0 = 1.0;
            if (vel_norm > EPS10)
              sigma = (h*h*grad_norm/(vel_norm*phi0))*(1.0-(vel_norm_bhpar/vel_norm));
          }

          // set specific term to convective term
          specific_term = conv_phi;
        }
        break;
        case INPAR::SCATRA::assgd_docarmo:
        case INPAR::SCATRA::assgd_almeida:
        {
          // velocity norm
          const double vel_norm = convelint.Norm2();

          // get norm of velocity vector z_h
          const double vel_norm_zh = abs(scatrares/grad_norm);

          // parameter zeta differentiating approaches by doCarmo and Galeao (1991)
          // and Almeida and Silva (1997)
          double zeta = 0.0;
          if (scatrapara_->ASSGDType() == INPAR::SCATRA::assgd_docarmo)
            zeta = 1.0;
          else
          {
            if (abs(scatrares) > EPS10)
              zeta = std::max(1.0,(conv_phi/scatrares));
          }

          // compute sigma
          if (vel_norm_zh > EPS10)
            sigma = tau*std::max(0.0,(vel_norm/vel_norm_zh)-zeta);

          // set specific term to residual
          specific_term = scatrares;
        }
        break;
        default: dserror("unknown type of all-scale subgrid diffusivity\n"); break;
      } //switch (whichassgd)

      // computation of subgrid diffusivity
      artdiff = sigma*scatrares*specific_term/(grad_norm*grad_norm);
      if (artdiff < 0.0)
      {
//        std::cout << "WARNING: isotropic artificial diffusion sgdiff < 0.0\n";
//        std::cout << "         -> set sgdiff to abs(sgdiff)!" << std::endl;
        artdiff = abs(sigma*scatrares*specific_term/(grad_norm*grad_norm));
      }
    }
    else artdiff = 0.0;
  }

//  if (artdiff>EPS8)
//    std::cout<<__FILE__<<__LINE__<<"\t artdiff=\t"<<artdiff<<std::endl;
  diffmanager_->SetIsotropicSubGridDiff(artdiff,k);

  return;
} //ScaTraEleCalc::CalcSubgrDiff


/*-------------------------------------------------------------------------------*
 |  calculation of strong residual for stabilization                             |
 | (depending on respective stationary or time-integration scheme)      vg 10/11 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcStrongResidual(
  const int                   k,          //!< index of current scalar
  double&                     scatrares,  //!< residual of convection-diffusion-reaction eq
  const double                densam,     //!< density at t_(n+am)
  const double                densnp,     //!< density at t_(n+1)
  const double                rea_phi,    //!< reactive contribution
  const double                rhsint,     //!< rhs at gauss point
  const double                tau         //!< the stabilisation parameter
  )
{
  // scalar at t_(n+1)
  const double   phinp = scatravarmanager_->Phinp(k);
  // history of time integration
  const double   hist = scatravarmanager_->Hist(k);
  // convective contribution
  const double   conv_phi = scatravarmanager_->ConvPhi(k);

  // diffusive part used in stabilization terms
  double diff_phi(0.0);
  LINALG::Matrix<nen_,1> diff(true);

  // diffusive term using current scalar value for higher-order elements
  // Note: has to be recomputed here every time, since the diffusion coefficient may have changed since the last call
  if (use2ndderiv_)
  {
    // diffusive part:  diffus * ( N,xx  +  N,yy +  N,zz )
    GetLaplacianStrongForm(diff);
    diff.Scale(diffmanager_->GetIsotropicDiff(k));
    diff_phi = diff.Dot(ephinp_[k]);
  }

  if (scatraparatimint_->IsGenAlpha())
  {
    // time derivative stored on history variable
    scatrares  = densam*hist + densnp*conv_phi
                     - diff_phi + rea_phi - rhsint;
  }
  else
  {
    // stationary residual
    scatrares = densnp*conv_phi - diff_phi + rea_phi - rhsint;

    if (not scatraparatimint_->IsStationary())
    {
      scatrares *= scatraparatimint_->TimeFac()/scatraparatimint_->Dt();
      scatrares += densnp*(phinp - hist)/scatraparatimint_->Dt();
    }
  }

  return;
} //ScaTraEleCalc::CalcStrongResidual


/*----------------------------------------------------------------------*
 |  calculate subgrid-scale velocity                           vg 10/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcSubgrVelocity(
  const DRT::Element*  ele,        //!< the element we are dealing with
  LINALG::Matrix<nsd_,1>& sgvelint, //!< subgrid velocity at integration point
  const double   densam,     //!< density at t_(n+am)
  const double   densnp,     //!< density at t_(n+1)
  const double   visc,       //!< fluid viscosity
  const LINALG::Matrix<nsd_,1>& convelint, //!< convective velocity at integration point
  const double   tau                       //!< the stabilisation parameter
  )
{
  // definitions
  LINALG::Matrix<nsd_,1> acc;
  LINALG::Matrix<nsd_,nsd_> vderxy;
  LINALG::Matrix<nsd_,1> conv;
  LINALG::Matrix<nsd_,1> gradp;
  LINALG::Matrix<nsd_,1> epsilonvel;
  LINALG::Matrix<nsd_,1> bodyforce;
  LINALG::Matrix<nsd_,1> pressuregrad;
  LINALG::Matrix<nsd_,nen_> nodebodyforce;
  LINALG::Matrix<nsd_,nen_> nodepressuregrad;

  // get acceleration or momentum history data
  acc.Multiply(eaccnp_,funct_);

  // get velocity derivatives
  vderxy.MultiplyNT(evelnp_,derxy_);

  // compute convective fluid term
  conv.Multiply(vderxy,convelint);

  // get pressure gradient
  gradp.Multiply(derxy_,eprenp_);

  //--------------------------------------------------------------------
  // get nodal values of fluid body force
  //--------------------------------------------------------------------
  std::vector<DRT::Condition*> myfluidneumcond;

  // check whether all nodes have a unique Fluid Neumann condition
  switch(nsd_)
  {
  case 3:
    DRT::UTILS::FindElementConditions(ele, "FluidVolumeNeumann", myfluidneumcond);
    break;
  case 2:
    DRT::UTILS::FindElementConditions(ele, "FluidSurfaceNeumann", myfluidneumcond);
    break;
  case 1:
    DRT::UTILS::FindElementConditions(ele, "FluidLineNeumann", myfluidneumcond);
    break;
  default:
    dserror("Illegal number of space dimensions: %d",nsd_); break;
  }

  if (myfluidneumcond.size()>1)
    dserror("more than one Fluid Neumann condition on one node");

  if (myfluidneumcond.size()==1)
  {
    const std::string* condtype = myfluidneumcond[0]->Get<std::string>("type");

    // get values and switches from the condition
    const std::vector<int>*    onoff = myfluidneumcond[0]->Get<std::vector<int> >   ("onoff");
    const std::vector<double>* val   = myfluidneumcond[0]->Get<std::vector<double> >("val"  );
    const std::vector<int>*    funct = myfluidneumcond[0]->Get<std::vector<int> >   ("funct");

    // factor given by spatial function
    double functfac = 1.0;
    int functnum = -1;

    // set this condition to the body-force array
    for (unsigned isd = 0; isd < nsd_; isd++)
    {
      // get factor given by spatial function
      if (funct) functnum = (*funct)[isd];
      else       functnum = -1;

      double num = (*onoff)[isd]*(*val)[isd];

      for (unsigned jnode = 0; jnode < nen_; ++jnode )
      {
        if (functnum > 0)
        {
          // time factor for the intermediate step
          // (negative time value indicates error)
          if (scatraparatimint_->Time() >= 0.0)
          {
            // evaluate function at the position of the current node
            // ------------------------------------------------------
            // comment: this introduces an additional error compared to an
            // evaluation at the integration point. However, we need a node
            // based element bodyforce vector for prescribed pressure gradients
            // in some fancy turbulance stuff.
            functfac = DRT::Problem::Instance()->Funct(functnum-1).Evaluate(isd,
                                                                            (ele->Nodes()[jnode])->X(),
                                                                             scatraparatimint_->Time());
          }
          else dserror("Negative time value in body force calculation: time = %f",scatraparatimint_->Time());

        }
        else functfac = 1.0;

        // compute body force
        if (*condtype == "neum_dead" or *condtype == "neum_live") nodebodyforce(isd,jnode) = num*functfac;
        else nodebodyforce.Clear();

        // compute prescribed pressure gradient
        if (*condtype == "neum_pgrad") nodepressuregrad(isd,jnode) = num*functfac;
        else nodepressuregrad.Clear();
      }
    }
  }
  else
  {
    nodebodyforce.Clear();
    nodepressuregrad.Clear();
  }

  // get fluid body force
  bodyforce.Multiply(nodebodyforce,funct_);
  // or prescribed pressure gradient
  pressuregrad.Multiply(nodepressuregrad,funct_);

  // get viscous term
  if (use2ndderiv_)
    /*--- viscous term: div(epsilon(u)) --------------------------------*/
    /*   /                                                \
         |  2 N_x,xx + N_x,yy + N_y,xy + N_x,zz + N_z,xz  |
         1 |                                                |
         - |  N_y,xx + N_x,yx + 2 N_y,yy + N_z,yz + N_y,zz  |
         2 |                                                |
         |  N_z,xx + N_x,zx + N_y,zy + N_z,yy + 2 N_z,zz  |
         \                                                /

         with N_x .. x-line of N
         N_y .. y-line of N                                             */

    /*--- subtraction for low-Mach-number flow: div((1/3)*(div u)*I) */
    /*   /                            \
         |  N_x,xx + N_y,yx + N_z,zx  |
         1 |                            |
         -  - |  N_x,xy + N_y,yy + N_z,zy  |
         3 |                            |
         |  N_x,xz + N_y,yz + N_z,zz  |
         \                            /

         with N_x .. x-line of N
         N_y .. y-line of N                                             */
    CalcSubgrVelocityVisc(epsilonvel);
  else
    epsilonvel.Clear();

  //--------------------------------------------------------------------
  // calculation of subgrid-scale velocity based on momentum residual
  // and stabilization parameter
  // (different for generalized-alpha and other time-integration schemes)
  //--------------------------------------------------------------------
  if (scatraparatimint_->IsGenAlpha())
  {
    for (unsigned rr=0;rr<nsd_;++rr)
    {
      sgvelint(rr) = -tau*(densam*acc(rr)+densnp*conv(rr)
                                +gradp(rr)-2*visc*epsilonvel(rr)
                                -densnp*bodyforce(rr)-pressuregrad(rr));
    }
  }
  else
  {
    for (unsigned rr=0;rr<nsd_;++rr)
    {
      sgvelint(rr) = -tau*(densnp*convelint(rr)+scatraparatimint_->TimeFac()*(densnp*conv(rr)
                                                                   +gradp(rr)-2*visc*epsilonvel(rr)
                                                                   -densnp*bodyforce(rr)-pressuregrad(rr))
                                                                   -densnp*acc(rr))/scatraparatimint_->Dt();
    }
  }

  return;
} //ScaTraEleCalc::CalcSubgrVelocity


/*---------------------------------------------------------------*
 | calculate viscous part of subgrid-scale velocity   fang 02/15 |
 *---------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype,int probdim>
void DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcSubgrVelocityVisc(
    LINALG::Matrix<nsd_,1>&   epsilonvel
    )
{
  if(nsd_ == 3)
  {
    for (unsigned i=0; i<nen_; ++i)
    {
      double sum = (derxy2_(0,i)+derxy2_(1,i)+derxy2_(2,i));

      epsilonvel(0) += (sum*evelnp_(0,i))/2.0;
      epsilonvel(1) += (sum*evelnp_(1,i))/2.0;
      epsilonvel(2) += (sum*evelnp_(2,i))/2.0;
    }
  }

  else if(nsd_ == 2)
  {
    for (unsigned i=0; i<nen_; ++i)
    {
      double sum = (derxy2_(0,i)+derxy2_(1,i));

      epsilonvel(0) += (sum*evelnp_(0,i))/2.0;
      epsilonvel(1) += (sum*evelnp_(1,i))/2.0;
    }
  }

  else
   dserror("Epsilon(u) is not implemented for the 1D case!");

  return;
} // DRT::ELEMENTS::ScaTraEleCalc<distype,probdim>::CalcSubgrVelocityVisc


// template classes

#include "scatra_ele_calc_fwd.hpp"
