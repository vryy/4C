/*!----------------------------------------------------------------------
\file fluid_ele_calc_immersed.cpp

\brief calc class for immersed problems

\level 3

<pre>
\maintainer  Andreas Rauch
             rauch@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289--15240
</pre>
*----------------------------------------------------------------------*/

#include "fluid_ele_calc_immersed.H"

#include "fluid_ele_tds.H"
#include "fluid_ele_parameter_std.H"
#include "../drt_fluid_ele/fluid_ele_immersed_base.H"
#include "../drt_lib/drt_globalproblem.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcImmersed<distype>*
DRT::ELEMENTS::FluidEleCalcImmersed<distype>::Instance(bool create)
{
  static FluidEleCalcImmersed<distype>* instance;
  if (create)
  {
    if (instance == NULL)
    {
      instance = new FluidEleCalcImmersed<distype>();
    }
  }
  else
  {
    if (instance != NULL) delete instance;
    instance = NULL;
  }
  return instance;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcImmersed<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance(false);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidEleCalcImmersed<distype>::FluidEleCalcImmersed()
    : DRT::ELEMENTS::FluidEleCalc<distype>::FluidEleCalc(), immersedele_(NULL), gp_iquad_(0)
{
  my::fldpara_ = DRT::ELEMENTS::FluidEleParameterStd::Instance();
}

/*----------------------------------------------------------------------*
 * Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalcImmersed<distype>::Evaluate(DRT::ELEMENTS::Fluid* ele,
    DRT::Discretization& discretization, const std::vector<int>& lm, Teuchos::ParameterList& params,
    Teuchos::RCP<MAT::Material>& mat, Epetra_SerialDenseMatrix& elemat1_epetra,
    Epetra_SerialDenseMatrix& elemat2_epetra, Epetra_SerialDenseVector& elevec1_epetra,
    Epetra_SerialDenseVector& elevec2_epetra, Epetra_SerialDenseVector& elevec3_epetra,
    bool offdiag)
{
  // get integration rule for fluid elements cut by structural boundary
  int num_gp_fluid_bound =
      DRT::Problem::Instance()->ImmersedMethodParams().get<int>("NUM_GP_FLUID_BOUND");
  int degree_gp_fluid_bound = 3;
  if (num_gp_fluid_bound == 8)
    degree_gp_fluid_bound = 3;
  else if (num_gp_fluid_bound == 64)
    degree_gp_fluid_bound = 7;
  else if (num_gp_fluid_bound == 125)
    degree_gp_fluid_bound = 9;
  else if (num_gp_fluid_bound == 343)
    degree_gp_fluid_bound = 13;
  else if (num_gp_fluid_bound == 729)
    degree_gp_fluid_bound = 17;
  else if (num_gp_fluid_bound == 1000)
    degree_gp_fluid_bound = 19;
  else
    dserror(
        "Invalid value for parameter NUM_GP_FLUID_BOUND (valid parameters are 8, 64, 125, 343, 729 "
        "and 1000).");

  // initialize integration rules
  const DRT::UTILS::GaussIntegration intpoints_fluid_bound(distype, degree_gp_fluid_bound);
  const DRT::UTILS::GaussIntegration intpoints_std(distype);

  // store current element
  immersedele_ = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(ele);

  // use different integration rule for fluid elements that are cut by the structural boundary
  if (immersedele_->IsBoundaryImmersed())
  {
    return my::Evaluate(ele, discretization, lm, params, mat, elemat1_epetra, elemat2_epetra,
        elevec1_epetra, elevec2_epetra, elevec3_epetra, intpoints_fluid_bound, offdiag);
  }
  else
  {
    return my::Evaluate(ele, discretization, lm, params, mat, elemat1_epetra, elemat2_epetra,
        elevec1_epetra, elevec2_epetra, elevec3_epetra, intpoints_std, offdiag);
  }
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcImmersed<distype>::ComputeSubgridScaleVelocity(
    const LINALG::Matrix<my::nsd_, my::nen_>& eaccam, double& fac1, double& fac2, double& fac3,
    double& facMtau, int iquad, double* saccn, double* sveln, double* svelnp)
{
  // set number of current gp
  gp_iquad_ = iquad;
  // compute convective conservative term from previous iteration u_old(u_old*nabla)
  LINALG::Matrix<my::nsd_, 1> conv_old_cons(true);
  conv_old_cons.Update(my::vdiv_, my::convvelint_, 0.0);

  //----------------------------------------------------------------------
  // compute residual of momentum equation
  // -> different for generalized-alpha and other time-integration schemes
  //----------------------------------------------------------------------
  if (my::fldparatimint_->IsGenalpha())
  {
    // rhs of momentum equation: density*bodyforce at n+alpha_F
    if (my::fldpara_->PhysicalType() == INPAR::FLUID::boussinesq)
    {
      // safety check
      if (my::fldparatimint_->AlphaF() != 1.0 or my::fldparatimint_->Gamma() != 1.0)
        dserror(
            "Boussinesq approximation in combination with generalized-alpha time integration "
            "has only been tested for BDF2-equivalent time integration parameters! "
            "Feel free to remove this error at your own risk!");

      my::rhsmom_.Update(my::deltadens_, my::bodyforce_, 0.0);
    }
    else
      my::rhsmom_.Update(my::densaf_, my::bodyforce_, 0.0);

    // add pressure gradient prescribed as body force (caution: not density weighted)
    my::rhsmom_.Update(1.0, my::generalbodyforce_, 1.0);

    // get acceleration at time n+alpha_M at integration point
    my::accint_.Multiply(eaccam, my::funct_);

    // evaluate momentum residual once for all stabilization right hand sides
    for (int rr = 0; rr < my::nsd_; ++rr)
    {
      if (immersedele_->HasProjectedDirichlet())
        my::momres_old_(rr) = my::densam_ * my::accint_(rr) +
                              my::densaf_ * (my::conv_old_(rr) + conv_old_cons(rr)) +
                              my::gradp_(rr) - 2 * my::visceff_ * my::visc_old_(rr) +
                              my::reacoeff_ * my::velint_(rr) - my::rhsmom_(rr);
      else
        my::momres_old_(rr) = my::densam_ * my::accint_(rr) + my::densaf_ * my::conv_old_(rr) +
                              my::gradp_(rr) - 2 * my::visceff_ * my::visc_old_(rr) +
                              my::reacoeff_ * my::velint_(rr) - my::rhsmom_(rr);
    }

    // add consistency terms for MFS if applicable
    my::MultfracSubGridScalesConsistentResidual();
  }
  else
  {
    if (not my::fldparatimint_->IsStationary())
    {
      // rhs of instationary momentum equation:
      // density*theta*bodyforce at n+1 + density*(histmom/dt)
      // in the case of a Boussinesq approximation: f = rho_0*[(rho - rho_0)/rho_0]*g = (rho -
      // rho_0)*g else:                                      f = rho * g Changed density from densn_
      // to densaf_. Makes the OST consistent with the gen-alpha.
      if (my::fldpara_->PhysicalType() == INPAR::FLUID::boussinesq)
        my::rhsmom_.Update((my::densaf_ / my::fldparatimint_->Dt() / my::fldparatimint_->Theta()),
            my::histmom_, my::deltadens_, my::bodyforce_);
      else
        my::rhsmom_.Update((my::densaf_ / my::fldparatimint_->Dt() / my::fldparatimint_->Theta()),
            my::histmom_, my::densaf_, my::bodyforce_);

      // add pressure gradient prescribed as body force (caution: not density weighted)
      my::rhsmom_.Update(1.0, my::generalbodyforce_, 1.0);

      // compute instationary momentum residual:
      // momres_old = u_(n+1)/dt + theta ( ... ) - histmom_/dt - theta*bodyforce_
      for (int rr = 0; rr < my::nsd_; ++rr)
      {
        if (immersedele_->HasProjectedDirichlet())
          my::momres_old_(rr) = ((my::densaf_ * my::velint_(rr) / my::fldparatimint_->Dt() +
                                     my::fldparatimint_->Theta() *
                                         (my::densaf_ * (my::conv_old_(rr) + conv_old_cons(rr)) +
                                             my::gradp_(rr) - 2 * my::visceff_ * my::visc_old_(rr) +
                                             my::reacoeff_ * my::velint_(rr))) /
                                    my::fldparatimint_->Theta()) -
                                my::rhsmom_(rr);
        else
          my::momres_old_(rr) =
              ((my::densaf_ * my::velint_(rr) / my::fldparatimint_->Dt() +
                   my::fldparatimint_->Theta() * (my::densaf_ * my::conv_old_(rr) + my::gradp_(rr) -
                                                     2 * my::visceff_ * my::visc_old_(rr) +
                                                     my::reacoeff_ * my::velint_(rr))) /
                  my::fldparatimint_->Theta()) -
              my::rhsmom_(rr);
      }
    }
    else
    {
      // rhs of stationary momentum equation: density*bodyforce
      // in the case of a Boussinesq approximation: f = rho_0*[(rho - rho_0)/rho_0]*g = (rho -
      // rho_0)*g else:                                      f = rho * g and pressure gradient
      // prescribed as body force (not density weighted)
      if (my::fldpara_->PhysicalType() == INPAR::FLUID::boussinesq)
        my::rhsmom_.Update(my::deltadens_, my::bodyforce_, 1.0, my::generalbodyforce_);
      else
        my::rhsmom_.Update(my::densaf_, my::bodyforce_, 1.0, my::generalbodyforce_);

      // compute stationary momentum residual:
      for (int rr = 0; rr < my::nsd_; ++rr)
      {
        if (immersedele_->HasProjectedDirichlet())
          my::momres_old_(rr) = my::densaf_ * (my::conv_old_(rr) + conv_old_cons(rr)) +
                                my::gradp_(rr) - 2 * my::visceff_ * my::visc_old_(rr) +
                                my::reacoeff_ * my::velint_(rr) - my::rhsmom_(rr);
        else
          my::momres_old_(rr) = my::densaf_ * my::conv_old_(rr) + my::gradp_(rr) -
                                2 * my::visceff_ * my::visc_old_(rr) +
                                my::reacoeff_ * my::velint_(rr) - my::rhsmom_(rr);
      }

      // add consistency terms for MFS if applicable
      my::MultfracSubGridScalesConsistentResidual();
    }
  }

  //----------------------------------------------------------------------
  // compute subgrid-scale velocity
  //----------------------------------------------------------------------
  // 1) quasi-static subgrid scales
  // Definition of subgrid-scale velocity is not consistent for the SUPG term and Franca, Valentin,
  // ... Definition of subgrid velocity used by Hughes
  if (my::fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
  {
    my::sgvelint_.Update(-my::tau_(1), my::momres_old_, 0.0);
  }
  // 2) time-dependent subgrid scales
  else
  {
    // some checking
    if (my::fldparatimint_->IsStationary())
      dserror("there is no time dependent subgrid scale closure for stationary problems\n");
    if (saccn == NULL or sveln == NULL or svelnp == NULL) dserror("no subscale array provided");

    // parameter definitions
    double alphaF = my::fldparatimint_->AlphaF();
    double alphaM = my::fldparatimint_->AlphaM();
    double gamma = my::fldparatimint_->Gamma();
    double dt = my::fldparatimint_->Dt();

    /*
                                            1.0
       facMtau =  -------------------------------------------------------
                     n+aM                      n+aF
                  rho     * alphaM * tauM + rho     * alphaF * gamma * dt
    */
    facMtau =
        1.0 / (my::densam_ * alphaM * my::tau_(1) + my::densaf_ * my::fldparatimint_->Afgdt());

    /*
       factor for old subgrid velocities:

                 n+aM                      n+aF
       fac1 = rho     * alphaM * tauM + rho     * gamma * dt * (alphaF-1)
    */
    fac1 =
        (my::densam_ * alphaM * my::tau_(1) + my::densaf_ * gamma * dt * (alphaF - 1.0)) * facMtau;
    /*
      factor for old subgrid accelerations

                 n+aM
       fac2 = rho     * tauM * dt * (alphaM-gamma)
    */
    fac2 = (my::densam_ * dt * my::tau_(1) * (alphaM - gamma)) * facMtau;
    /*
      factor for residual in current subgrid velocities:

       fac3 = gamma * dt * tauM
    */
    fac3 = (gamma * dt * my::tau_(1)) * facMtau;

    // warning: time-dependent subgrid closure requires generalized-alpha time
    // integration
    if (!my::fldparatimint_->IsGenalpha())
    {
      dserror("the time-dependent subgrid closure requires a genalpha time integration\n");
    }

    /*         +-                                       -+
        ~n+1   |        ~n           ~ n            n+1  |
        u    = | fac1 * u  + fac2 * acc  -fac3 * res     |
         (i)   |                                    (i)  |
               +-                                       -+
    */

    /* compute the intermediate value of subscale velocity

            ~n+af            ~n+1                   ~n
            u     = alphaF * u     + (1.0-alphaF) * u
             (i)              (i)

    */

    static LINALG::Matrix<1, my::nsd_> sgvelintaf(true);
    sgvelintaf.Clear();
    for (int rr = 0; rr < my::nsd_; ++rr)
    {
      my::tds_->UpdateSvelnpInOneDirection(fac1, fac2, fac3, my::momres_old_(rr),
          my::fldparatimint_->AlphaF(), rr, iquad,
          my::sgvelint_(
              rr),  // sgvelint_ is set to sgvelintnp, but is then overwritten below anyway!
          sgvelintaf(rr));

      int pos = rr + my::nsd_ * iquad;

      /*
       *  ~n+1           ~n           ~ n            n+1
       *  u    =  fac1 * u  + fac2 * acc  -fac3 * res
       *   (i)
       *
       */

      svelnp[pos] = fac1 * sveln[pos] + fac2 * saccn[pos] - fac3 * my::momres_old_(rr);

      /* compute the intermediate value of subscale velocity
       *
       *          ~n+af            ~n+1                   ~n
       *          u     = alphaF * u     + (1.0-alphaF) * u
       *           (i)              (i)
       *
       */
      my::sgvelint_(rr) = alphaF * svelnp[pos] + (1.0 - alphaF) * sveln[pos];
    }
  }  // end time dependent subgrid scale closure

  //----------------------------------------------------------------------
  // include computed subgrid-scale velocity in convective term
  // -> only required for cross- and Reynolds-stress terms
  //----------------------------------------------------------------------
  if (my::fldpara_->Cross() != INPAR::FLUID::cross_stress_stab_none or
      my::fldpara_->Reynolds() != INPAR::FLUID::reynolds_stress_stab_none or
      my::fldpara_->ContiCross() != INPAR::FLUID::cross_stress_stab_none or
      my::fldpara_->ContiReynolds() != INPAR::FLUID::reynolds_stress_stab_none)
    my::sgconv_c_.MultiplyTN(my::derxy_, my::sgvelint_);
  else
    my::sgconv_c_.Clear();

  return;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcImmersed<distype>::LinGalMomResU(
    LINALG::Matrix<my::nsd_ * my::nsd_, my::nen_>& lin_resM_Du, const double& timefacfac)
{
  /*
      instationary                                 conservative          cross-stress, part 1
       +-----+                             +-------------------------+   +-------------------+
       |     |                             |                         |   |                   |

                 /       n+1       \           /               n+1  \     /      ~n+1       \
       rho*Du + |   rho*u   o nabla | Du  + Du |   rho*nabla o u     | + |   rho*u   o nabla | Du +
                 \      (i)        /           \               (i)  /     \       (i)        /

                 /                \  n+1     n+1  /                 \
              + |   rho*Du o nabla | u    + u    | rho*nabla o Du    |  + sigma*Du
                 \                /   (i)    (i)  \                 /
                |                        | |                         |    |       |
                +------------------------+ +-------------------------+    +-------+
                        Newton                Newton (conservative)       reaction
  */


  // convective form of momentum residual
  my::LinGalMomResU(lin_resM_Du, timefacfac);

  // add conservative terms for covered fluid integration points
  const double timefacfac_densaf = timefacfac * my::densaf_;
  int idim_nsd_p_idim[my::nsd_];
  for (int idim = 0; idim < my::nsd_; ++idim)
  {
    idim_nsd_p_idim[idim] = idim * my::nsd_ + idim;
  }


  // convection, convective part (conservative addition)
  if (immersedele_->HasProjectedDirichlet())
  {
    for (int ui = 0; ui < my::nen_; ++ui)
    {
      const double v = timefacfac_densaf * my::vdiv_;

      for (int idim = 0; idim < my::nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim], ui) += my::funct_(ui) * v;
      }
    }

    //  convection, reactive part (conservative addition) (only for Newton)
    if (my::fldpara_->IsNewton())
    {
      for (int ui = 0; ui < my::nen_; ++ui)
      {
        for (int idim = 0; idim < my::nsd_; ++idim)
        {
          const double temp = timefacfac_densaf * my::velint_(idim);
          const int idim_nsd = idim * my::nsd_;

          for (int jdim = 0; jdim < my::nsd_; ++jdim)
          {
            lin_resM_Du(idim_nsd + jdim, ui) += temp * my::derxy_(jdim, ui);
          }
        }
      }
    }
  }

  return;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcImmersed<distype>::InertiaConvectionReactionGalPart(
    LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_>& estif_u,
    LINALG::Matrix<my::nsd_, my::nen_>& velforce,
    LINALG::Matrix<my::nsd_ * my::nsd_, my::nen_>& lin_resM_Du,
    LINALG::Matrix<my::nsd_, 1>& resM_Du, const double& rhsfac)
{
  my::InertiaConvectionReactionGalPart(estif_u, velforce, lin_resM_Du, resM_Du, rhsfac);

  if (immersedele_->HasProjectedDirichlet())
  {
    for (int idim = 0; idim < my::nsd_; ++idim)
    {
      /* convection (conservative addition) on right-hand side */
      double v = -rhsfac * my::densaf_ * my::velint_(idim) * my::vdiv_;

      if (my::fldpara_->PhysicalType() == INPAR::FLUID::loma)
        v += rhsfac * my::velint_(idim) * my::densaf_ * my::scaconvfacaf_ * my::conv_scaaf_;
      else if (my::fldpara_->PhysicalType() == INPAR::FLUID::varying_density)
        v -= rhsfac * my::velint_(idim) * my::conv_scaaf_;

      for (int vi = 0; vi < my::nen_; ++vi) velforce(idim, vi) += v * my::funct_(vi);
    }
  }

  return;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcImmersed<distype>::ContinuityGalPart(
    LINALG::Matrix<my::nen_, my::nen_ * my::nsd_>& estif_q_u, LINALG::Matrix<my::nen_, 1>& preforce,
    const double& timefacfac, const double& timefacfacpre, const double& rhsfac)
{
  for (int vi = 0; vi < my::nen_; ++vi)
  {
    const double v = timefacfacpre * my::funct_(vi);
    for (int ui = 0; ui < my::nen_; ++ui)
    {
      const int fui = my::nsd_ * ui;

      for (int idim = 0; idim < my::nsd_; ++idim)
      {
        /* continuity term */
        /*
             /                \
            |                  |
            | nabla o Du  , q  |
            |                  |
             \                /
        */
        estif_q_u(vi, fui + idim) += v * my::derxy_(idim, ui);
      }
    }
  }  // end for(idim)

  // continuity term on right-hand side
  if (not immersedele_->IsImmersed()) preforce.Update(-rhsfac * my::vdiv_, my::funct_, 1.0);

  if (immersedele_->IsBoundaryImmersed())
    preforce.Update(rhsfac * immersedele_->ProjectedIntPointDivergence(gp_iquad_), my::funct_, 1.0);


  return;
}

template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalcImmersed<distype>::ConservativeFormulation(
    LINALG::Matrix<my::nen_ * my::nsd_, my::nen_ * my::nsd_>& estif_u,
    LINALG::Matrix<my::nsd_, my::nen_>& velforce, const double& timefacfac, const double& rhsfac)
{
  if (not(immersedele_->HasProjectedDirichlet()))
    my::ConservativeFormulation(estif_u, velforce, timefacfac, rhsfac);
  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// Ursula is responsible for this comment!
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::wedge15>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalcImmersed<DRT::Element::nurbs27>;
