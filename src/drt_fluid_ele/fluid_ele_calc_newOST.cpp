/*----------------------------------------------------------------------*/
/*!

\brief main file containing routines for calculation of fluid element
       new one-step theta time integration variant

\maintainer Martin Kronbichler

\level 3

*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc.H"
#include "fluid_ele_parameter.H"
#include "fluid_ele_parameter_std.H"
#include "fluid_ele_parameter_timint.H"
#include "fluid_ele.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "fluid_ele_tds.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_fem_general/drt_utils_nurbs_shapefunctions.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"

#include "../drt_geometry/position_array.H"

#include "../drt_inpar/inpar_turbulence.H"
#include "../drt_inpar/inpar_topopt.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/standardtypes_cpp.H"

#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/herschelbulkley.H"
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/modpowerlaw.H"
#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/permeablefluid.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/tempdepwater.H"
#include "../drt_mat/yoghurt.H"
#include "../drt_mat/matlist.H"

#include "../drt_nurbs_discret/drt_nurbs_utils.H"


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)   g.bau 03/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::SysmatOSTNew(
    const LINALG::Matrix<nsd_, nen_>& ebofoaf, const LINALG::Matrix<nsd_, nen_>& eprescpgaf,
    const LINALG::Matrix<nsd_, nen_>& ebofon, const LINALG::Matrix<nsd_, nen_>& eprescpgn,
    const LINALG::Matrix<nsd_, nen_>& evelaf, const LINALG::Matrix<nsd_, nen_>& evelam,
    const LINALG::Matrix<nsd_, nen_>& eveln, const LINALG::Matrix<nsd_, nen_>& evelnp,
    const LINALG::Matrix<nsd_, nen_>& fsevelaf, const LINALG::Matrix<nen_, 1>& fsescaaf,
    const LINALG::Matrix<nsd_, nen_>& evel_hat,
    const LINALG::Matrix<nsd_ * nsd_, nen_>& ereynoldsstress_hat,
    const LINALG::Matrix<nen_, 1>& epreaf, const LINALG::Matrix<nen_, 1>& epream,
    const LINALG::Matrix<nen_, 1>& epren, const LINALG::Matrix<nen_, 1>& eprenp,
    const LINALG::Matrix<nsd_, nen_>& eaccam, const LINALG::Matrix<nen_, 1>& escaaf,
    const LINALG::Matrix<nen_, 1>& escaam, const LINALG::Matrix<nen_, 1>& escadtam,
    const LINALG::Matrix<nen_, 1>& escabofoaf, const LINALG::Matrix<nen_, 1>& escabofon,
    const LINALG::Matrix<nsd_, nen_>& emhist, const LINALG::Matrix<nsd_, nen_>& edispnp,
    const LINALG::Matrix<nsd_, nen_>& egridv, const LINALG::Matrix<nsd_, nen_>& egridvn,
    LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& estif,
    LINALG::Matrix<(nsd_ + 1) * nen_, (nsd_ + 1) * nen_>& emesh,
    LINALG::Matrix<(nsd_ + 1) * nen_, 1>& eforce, const LINALG::Matrix<nen_, 1>& eporo,
    const LINALG::Matrix<nsd_, 2 * nen_>& egradphi, const LINALG::Matrix<nen_, 2 * 1>& ecurvature,
    const double thermpressaf, const double thermpressam, const double thermpressdtaf,
    const double thermpressdtam, Teuchos::RCP<const MAT::Material> material, double& Cs_delta_sq,
    double& Ci_delta_sq, double& Cv, bool isale, double* saccn, double* sveln, double* svelnp,
    const DRT::UTILS::GaussIntegration& intpoints)
{
  //------------------------------------------------------------------------
  //  preliminary definitions and evaluations
  //------------------------------------------------------------------------
  // definition of matrices
  estif_u_.Clear();
  estif_p_v_.Clear();
  estif_q_u_.Clear();
  ppmat_.Clear();

  // definition of vectors
  preforce_.Clear();
  velforce_.Clear();

  // definition of velocity-based momentum residual vectors
  lin_resM_Du_.Clear();
  resM_Du_.Clear();

  // if polynomial pressure projection: reset variables
  if (fldpara_->PPP())
  {
    D_ = 0;
    E_.Clear();
  }

  // evaluate shape functions and derivatives at element center
  EvalShapeFuncAndDerivsAtEleCenter();

  // set element area or volume
  const double vol = fac_;

  //------------------------------------------------------------------------
  // potential evaluation of material parameters, subgrid viscosity
  // and/or stabilization parameters at element center
  //------------------------------------------------------------------------
  // get material parameters at element center
  if (not fldpara_->MatGp() or not fldpara_->TauGp())
  {
    GetMaterialParams(material, evelaf, epreaf, epream, escaaf, escaam, escabofoaf, thermpressaf,
        thermpressam, thermpressdtaf, thermpressdtam, vol);

    // calculate all-scale or fine-scale subgrid viscosity at element center
    visceff_ = visc_;

    if (fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or
        fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky or
        fldpara_->TurbModAction() == INPAR::FLUID::vreman or
        fldpara_->TurbModAction() == INPAR::FLUID::dynamic_vreman)
    {
      if (fldpara_->TurbModAction() == INPAR::FLUID::dynamic_vreman)
        Cs_delta_sq = Cv;  // use the declaration of Cs_delta_sq for the dynamic Vreman constant
      CalcSubgrVisc(evelaf, vol, Cs_delta_sq, Ci_delta_sq);
      // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
      visceff_ += sgvisc_;
    }
    else if (fldpara_->Fssgv() != INPAR::FLUID::no_fssgv)
      CalcFineScaleSubgrVisc(evelaf, fsevelaf, vol);
  }

  // potential evaluation of multifractal subgrid-scales at element center
  // coefficient B of fine-scale velocity
  static LINALG::Matrix<nsd_, 1> B_mfs(true);
  B_mfs.Clear();

  // coefficient D of fine-scale scalar (loma only)
  double D_mfs = 0.0;
  if (fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales)
  {
    if (not fldpara_->BGp())
    {
      // make sure to get material parameters at element center
      if (fldpara_->MatGp())
        // GetMaterialParams(material,evelaf,epreaf,epream,escaaf,escaam,thermpressaf,thermpressam,thermpressdtam,vol);
        GetMaterialParams(material, evelaf, epreaf, epream, escaaf, escaam, escabofoaf,
            thermpressaf, thermpressam, thermpressdtaf, thermpressdtam, vol);

      // provide necessary velocities and gradients at element center
      velint_.Multiply(evelaf, funct_);
      fsvelint_.Multiply(fsevelaf, funct_);
      vderxy_.MultiplyNT(evelaf, derxy_);
      // calculate parameters of multifractal subgrid-scales and, finally,
      // calculate coefficient for multifractal modeling of subgrid velocity
      // if loma, calculate coefficient for multifractal modeling of subgrid scalar
      PrepareMultifractalSubgrScales(B_mfs, D_mfs, evelaf, fsevelaf, vol);
      // clear all velocities and gradients
      velint_.Clear();
      fsvelint_.Clear();
      vderxy_.Clear();
    }
  }


  // calculate stabilization parameter at element center
  if (not fldpara_->TauGp() and fldpara_->StabType() == INPAR::FLUID::stabtype_residualbased)
  {
    // get convective velocity at element center
    // for evaluation of stabilization parameter
    velint_.Multiply(evelaf, funct_);

    // get the grid velocity in case of ALE
    if (isale)
    {
      gridvelint_.Multiply(egridv, funct_);
      gridvelintn_.Multiply(egridvn, funct_);
    }

    // get convective velocity at integration point
    SetConvectiveVelint(isale);


    if (fldpara_->Tds() == INPAR::FLUID::subscales_time_dependent)
    {
      // get velocity derivatives at integration point
      // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
      vderxy_.MultiplyNT(evelaf, derxy_);  // required for time-dependent subscales

      // compute velnp at integration point (required for time-dependent subscales)
      static LINALG::Matrix<nsd_, 1> velintnp(true);
      velintnp.Multiply(evelnp, funct_);
      vel_normnp_ = velintnp.Norm2();
    }

    // calculate stabilization parameters at element center
    CalcStabParameter(vol);
  }

  // get Gaussian integration points
  // const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);
  // const DRT::UTILS::IntPointsAndWeights<nsd_>
  // intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  // for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)

  for (DRT::UTILS::GaussIntegration::const_iterator iquad = intpoints.begin();
       iquad != intpoints.end(); ++iquad)
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(), iquad.Weight());

    //----------------------------------------------------------------------
    //  evaluation of various values at integration point:
    //  1) velocity (including derivatives, fine-scale and grid velocity)
    //  2) pressure (including derivatives)
    //  3) body-force vector
    //  4) "history" vector for momentum equation
    //----------------------------------------------------------------------
    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    velint_.Multiply(evelaf, funct_);

    // get velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    vderxy_.MultiplyNT(evelaf, derxy_);

    // get fine-scale velocity and its derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    if (fldpara_->Fssgv() != INPAR::FLUID::no_fssgv)
    {
      fsvderxy_.MultiplyNT(fsevelaf, derxy_);
    }
    else
    {
      fsvderxy_.Clear();
    }
    if (fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales)
    {
      fsvelint_.Multiply(fsevelaf, funct_);
      fsvderxy_.MultiplyNT(fsevelaf, derxy_);
    }
    else
    {
      fsvelint_.Clear();
    }

    // get the grid velocity in case of ALE
    if (isale)
    {
      gridvelint_.Multiply(egridv, funct_);
      gridvelintn_.Multiply(egridvn, funct_);
    }

    // get convective velocity at integration point
    SetConvectiveVelint(isale);


    // get pressure at integration point
    // (value at n+alpha_F for generalized-alpha scheme,
    //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
    double press = 0.0;
    if (fldparatimint_->IsGenalphaNP())
      press = funct_.Dot(eprenp);
    else
      press = funct_.Dot(epreaf);

    // get pressure gradient at integration point
    // (value at n+alpha_F for generalized-alpha scheme,
    //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
    if (fldparatimint_->IsGenalphaNP())
      gradp_.Multiply(derxy_, eprenp);
    else
      gradp_.Multiply(derxy_, epreaf);

    // get bodyforce at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    bodyforce_.Multiply(ebofoaf, funct_);
    bodyforcen_.Multiply(ebofon, funct_);

    // get prescribed pressure gradient acting as body force
    // (required for turbulent channel flow)
    generalbodyforce_.Multiply(eprescpgaf, funct_);
    generalbodyforcen_.Multiply(eprescpgn, funct_);

    // get momentum history data at integration point
    // (only required for one-step-theta and BDF2 time-integration schemes)
    histmom_.Multiply(emhist, funct_);


    //----------------------------------------------------------------------
    // potential evaluation of material parameters, subgrid viscosity
    // and/or stabilization parameters at integration point
    //----------------------------------------------------------------------
    // get material parameters at integration point

    if (fldpara_->MatGp())
    {
      GetMaterialParams(material, evelaf, epreaf, epream, escaaf, escaam, escabofoaf, thermpressaf,
          thermpressam, thermpressdtaf, thermpressdtam, vol);

      // calculate all-scale or fine-scale subgrid viscosity at integration point
      visceff_ = visc_;

      if (fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or
          fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky or
          fldpara_->TurbModAction() == INPAR::FLUID::vreman)
      {
        CalcSubgrVisc(evelaf, vol, Cs_delta_sq, Ci_delta_sq);
        // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
        visceff_ += sgvisc_;
      }
      else if (fldpara_->Fssgv() != INPAR::FLUID::no_fssgv)
        CalcFineScaleSubgrVisc(evelaf, fsevelaf, vol);
    }

    // get reaction coefficient due to porosity for topology optimization
    // !do this only at gauss point since this is nonlinear!
    if (fldpara_->ReactionTopopt()) GetPorosityAtGP(eporo);

    // calculate stabilization parameter at integration point
    if (fldpara_->TauGp() and fldpara_->StabType() == INPAR::FLUID::stabtype_residualbased)
      CalcStabParameter(vol);

    // potential evaluation of coefficient of multifractal subgrid-scales at integration point
    if (fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales)
    {
      if (fldpara_->BGp())
      {
        // make sure to get material parameters at gauss point
        if (not fldpara_->MatGp())
        {
          // GetMaterialParams(material,evelaf,epreaf,epream,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam);
          // would overwrite materials at the element center, hence BGp() should always be combined
          // with MatGp()
          dserror(
              "evaluation of B and D at gauss-point should always be combined with evaluation "
              "material at gauss-point!");
        }

        // calculate parameters of multifractal subgrid-scales
        PrepareMultifractalSubgrScales(B_mfs, D_mfs, evelaf, fsevelaf, vol);
      }

      // calculate fine-scale velocity, its derivative and divergence for multifractal subgrid-scale
      // modeling
      for (int idim = 0; idim < nsd_; idim++)
        mffsvelint_(idim, 0) = fsvelint_(idim, 0) * B_mfs(idim, 0);

      for (int idim = 0; idim < nsd_; idim++)
      {
        for (int jdim = 0; jdim < nsd_; jdim++)
          mffsvderxy_(idim, jdim) = fsvderxy_(idim, jdim) * B_mfs(idim, 0);
      }

      mffsvdiv_ = mffsvderxy_(0, 0) + mffsvderxy_(1, 1) + mffsvderxy_(2, 2);

      // only required for variable-density flow at low Mach number
      if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
      {
        mfssgscaint_ = D_mfs * funct_.Dot(fsescaaf);
        grad_fsscaaf_.Multiply(derxy_, fsescaaf);
        for (int dim = 0; dim < nsd_; dim++) grad_fsscaaf_(dim, 0) *= D_mfs;
      }
      else
      {
        mfssgscaint_ = 0.0;
        grad_fsscaaf_.Clear();
      }

      if (isale) dserror("Multifractal subgrid-scales with ale not supported");
    }
    else
    {
      mffsvelint_.Clear();
      mffsvderxy_.Clear();
      mffsvdiv_ = 0.0;
    }

    // TODO: Need gradphi and curvature at time n. (Can be implemented for interface values at
    // timestep t^(n+1)) Adds surface tension force to the Gausspoint.
    // Note: has to be called after GetMaterialParams(), otherwise gamma_ is uninitialized!!
    if (fldpara_->GetIncludeSurfaceTension())
      AddSurfaceTensionForce(escaaf, escaam, egradphi, ecurvature);

    //----------------------------------------------------------------------
    //  evaluation of various partial operators at integration point
    //  1) convective term from previous iteration and convective operator
    //  2) viscous term from previous iteration and viscous operator
    //  3) divergence of velocity from previous iteration
    //----------------------------------------------------------------------

    // compute convective term from previous iteration and convective operator
    conv_old_.Multiply(vderxy_, convvelint_);
    conv_c_.MultiplyTN(derxy_, convvelint_);

    // compute viscous term from previous iteration and viscous operator
    if (is_higher_order_ele_)
      CalcDivEps(evelaf, eveln);
    else
    {
      visc_old_.Clear();
      visc_oldn_.Clear();
      viscs2_.Clear();
    }

    // compute divergence of velocity from previous iteration
    vdiv_ = 0.0;
    vdivn_ = 0.0;
    vderxyn_.MultiplyNT(eveln, derxy_);
    if (not fldparatimint_->IsGenalphaNP())
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        vdiv_ += vderxy_(idim, idim);
        vdivn_ += vderxyn_(idim, idim);
      }
    }
    else
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        // get vdiv at time n+1 for np_genalpha,
        static LINALG::Matrix<nsd_, nsd_> vderxy(true);
        vderxy.MultiplyNT(evelnp, derxy_);
        vdiv_ += vderxy(idim, idim);
      }
    }

    // New One Step Theta variables (for old time step):
    velintn_.Multiply(eveln, funct_);
    // get convective velocity at integration point
    SetConvectiveVelintN(isale);
    conv_oldn_.Multiply(vderxyn_, convvelintn_);
    const double pressn = funct_.Dot(epren);
    gradpn_.Multiply(derxy_, epren);


    //-----------------------------------------------------------------------
    //       |          timefac         |  timefacpre     |    timefacrhs   |
    // ----------------------------------------------------------------------
    // OST   |                        dt*theta                              |
    //-----------------------------------------------------------------------
    // BDF2  |                        2/3 * dt                              |
    //-----------------------------------------------------------------------
    // Af GA |          alphaF*gamma*dt/alphaM            | gamma*dt/alphaM |
    //----------------------------------------------------------------------
    // NP GA | alphaF*gamma*dt/alphaM   | gamma*dt/alphaM | gamma*dt/alphaM |
    //-----------------------------------------------------------------------
    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    //----------------------------------------------------------------------
    const double timefacfac = fldparatimint_->TimeFac() * fac_;
    const double timefacfacpre = fldparatimint_->TimeFacPre() * fac_;
    const double rhsfac = fldparatimint_->TimeFacRhs() * fac_;

    // For One Step Theta rhsfac is: \Delta t (1 - \theta)
    const double rhsfacn = (1 - fldparatimint_->Theta()) * fldparatimint_->Dt() * fac_;
    // For One Step Theta,
    // the density multiplied with the instationary term has to be evaluated at time = ( n + \theta
    // ).
    dens_theta_ = fldparatimint_->Theta() * densaf_ + (1 - fldparatimint_->Theta()) * densn_;

    //----------------------------------------------------------------------
    // computation of various subgrid-scale values and residuals
    //----------------------------------------------------------------------
    // compute residual of momentum equation and subgrid-scale velocity
    // -> residual of momentum equation different for generalized-alpha
    //    and other time-integration schemes
    double fac1 = 0.0;
    double fac2 = 0.0;
    double fac3 = 0.0;
    double facMtau = 0.0;
    ComputeSubgridScaleVelocityOSTNew(
        eaccam, fac1, fac2, fac3, facMtau, *iquad, saccn, sveln, svelnp);

    // compute residual of continuity equation
    // residual contains velocity divergence only for incompressible flow
    conres_old_ = vdiv_;
    conres_oldn_ = vdivn_;

    // following computations only required for variable-density flow at low Mach number
    if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
    {
      // compute additional Galerkin terms on right-hand side of continuity equation
      // -> different for generalized-alpha and other time-integration schemes
      ComputeGalRHSContEq(eveln, escaaf, escaam, escadtam, isale);

      // add to residual of continuity equation
      conres_old_ -= rhscon_;

      if (fldpara_->UpdateMat() or fldpara_->ContiSUPG() or
          fldpara_->ContiCross() != INPAR::FLUID::cross_stress_stab_none or
          fldpara_->ContiReynolds() != INPAR::FLUID::reynolds_stress_stab_none or
          fldpara_->MultiFracLomaConti())
      {
        // compute subgrid-scale part of scalar
        // -> different for generalized-alpha and other time-integration schemes
        ComputeSubgridScaleScalar(escaaf, escaam);

        // update material parameters including subgrid-scale part of scalar
        if (fldpara_->UpdateMat())
        {
          // since we update the viscosity in the next step, a potential subgrid-scale velocity
          // would be overwritten
          if (fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or
              fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky or
              fldpara_->TurbModAction() == INPAR::FLUID::vreman)
            dserror("No material update in combination with smagorinsky model!");

          if (fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales)
            UpdateMaterialParams(material, evelaf, epreaf, epream, escaaf, escaam, thermpressaf,
                thermpressam, mfssgscaint_);
          else
            UpdateMaterialParams(material, evelaf, epreaf, epream, escaaf, escaam, thermpressaf,
                thermpressam, sgscaint_);
          visceff_ = visc_;
          if (fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or
              fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky or
              fldpara_->TurbModAction() == INPAR::FLUID::vreman)
            visceff_ += sgvisc_;
        }

        // right-hand side of continuity equation based on updated material parameters
        // and including all stabilization terms
        // -> different for generalized-alpha and other time-integration schemes
        RecomputeGalAndComputeCrossRHSContEq();
      }
    }
    else if (fldpara_->PhysicalType() == INPAR::FLUID::artcomp)
    {
      // compute additional Galerkin terms on right-hand side of continuity equation
      // -> different for generalized-alpha and other time-integration schemes
      ComputeGalRHSContEqArtComp(epreaf, epren, escadtam);

      // add to residual of continuity equation
      conres_old_ -= rhscon_;
    }

    // set velocity-based momentum residual vectors to zero
    lin_resM_Du_.Clear();
    resM_Du_.Clear();

    // compute first version of velocity-based momentum residual containing
    // inertia term, convection term (convective and reactive part),
    // reaction term and cross-stress term
    LinGalMomResUOSTNew(lin_resM_Du_, timefacfac);

    // potentially rescale first version of velocity-based momentum residual
    if (fldpara_->Tds() == INPAR::FLUID::subscales_time_dependent &&
        fldpara_->Transient() == INPAR::FLUID::inertia_stab_keep)
    {
      LinGalMomResU_subscales(estif_p_v_, lin_resM_Du_, resM_Du_, timefacfac, facMtau);
    }

    //----------------------------------------------------------------------
    // computation of standard Galerkin and stabilization contributions to
    // element matrix and right-hand-side vector
    //----------------------------------------------------------------------
    // 1) standard Galerkin inertia, convection and reaction terms
    //    (convective and reactive part for convection term)
    //    as well as first part of cross-stress term on left-hand side
    InertiaConvectionReactionGalPart(estif_u_, velforce_, lin_resM_Du_, resM_Du_, rhsfac, rhsfacn);

    // 2) standard Galerkin viscous term
    //    (including viscous stress computation,
    //     excluding viscous part for low-Mach-number flow)
    static LINALG::Matrix<nsd_, nsd_> viscstress(true);
    viscstress.Clear();

    ViscousGalPart(estif_u_, velforce_, viscstress, timefacfac, rhsfac, rhsfacn);

    // 3) stabilization of continuity equation,
    //    standard Galerkin viscous part for low-Mach-number flow and
    //    right-hand-side part of standard Galerkin viscous term
    if (fldpara_->CStab() or fldpara_->PhysicalType() == INPAR::FLUID::loma)
      ContStab(estif_u_, velforce_, fldparatimint_->TimeFac(), timefacfac, timefacfacpre, rhsfac,
          rhsfacn);

    // 4) standard Galerkin pressure term
    PressureGalPart(
        estif_p_v_, velforce_, timefacfac, timefacfacpre, rhsfac, rhsfacn, press, pressn);

    // 5) standard Galerkin continuity term
    ContinuityGalPart(estif_q_u_, preforce_, timefacfac, timefacfacpre, rhsfac, rhsfacn);

    // 6) standard Galerkin bodyforce term on right-hand side
    BodyForceRhsTerm(velforce_, rhsfac, rhsfacn);

    // 7) additional standard Galerkin terms due to conservative formulation
    //    New One Step Theta not implemented for this as of yet!
    if (fldpara_->IsConservative())
    {
      ConservativeFormulation(estif_u_, velforce_, timefacfac, rhsfac);
    }

    // 8) additional standard Galerkin terms for low-Mach-number flow and
    //    artificial compressibility (only right-hand side in latter case)
    //    New One Step Theta not implemented for LOMA as of yet.
    if (fldpara_->PhysicalType() == INPAR::FLUID::loma or
        fldpara_->PhysicalType() == INPAR::FLUID::artcomp)
    {
      LomaGalPart(estif_q_u_, preforce_, timefacfac, rhsfac);
    }

    // 9) additional standard Galerkin term for temporal derivative of pressure
    //    in case of artificial compressibility (only left-hand side)
    if (fldpara_->PhysicalType() == INPAR::FLUID::artcomp and not fldparatimint_->IsStationary())
      ArtCompPressureInertiaGalPartandContStab(estif_p_v_, ppmat_);

    //----------------------------------------------------------------------
    // compute second version of velocity-based momentum residual containing
    // inertia term, convection term (convective and reactive part) and
    // viscous term
    //----------------------------------------------------------------------
    StabLinGalMomResU(lin_resM_Du_, timefacfac);

    // 10) PSPG term
    if (fldpara_->PSPG())
    {
      PSPGOSTNew(estif_q_u_, ppmat_, preforce_, lin_resM_Du_, fac3, timefacfac, timefacfacpre,
          rhsfac, *iquad);
    }

    // 11) SUPG term as well as first part of Reynolds-stress term on
    //     left-hand side and Reynolds-stress term on right-hand side
    if (fldpara_->SUPG())
    {
      SUPGOSTNew(estif_u_, estif_p_v_, velforce_, preforce_, lin_resM_Du_, fac3, timefacfac,
          timefacfacpre, rhsfac);
    }

    // 12) reactive stabilization term
    if (fldpara_->RStab() != INPAR::FLUID::reactive_stab_none)
    {
      ReacStab(
          estif_u_, estif_p_v_, velforce_, lin_resM_Du_, timefacfac, timefacfacpre, rhsfac, fac3);
    }

    // 13) viscous stabilization term
    if (is_higher_order_ele_ and (fldpara_->VStab() != INPAR::FLUID::viscous_stab_none))
    {
      ViscStab(
          estif_u_, estif_p_v_, velforce_, lin_resM_Du_, timefacfac, timefacfacpre, rhsfac, fac3);
    }

    // if ConvDivStab for XFEM
    //    {
    //      ConvDivStab(estif_u,
    //           velforce,
    //           timefacfac,
    //           rhsfac);
    //    }


    // 14) cross-stress term: second part on left-hand side (only for Newton
    //     iteration) as well as cross-stress term on right-hand side
    if (fldpara_->Cross() != INPAR::FLUID::cross_stress_stab_none)
    {
      CrossStressStab(
          estif_u_, estif_p_v_, velforce_, lin_resM_Du_, timefacfac, timefacfacpre, rhsfac, fac3);
    }

    // 15) Reynolds-stress term: second part on left-hand side
    //     (only for Newton iteration)
    if (fldpara_->Reynolds() == INPAR::FLUID::reynolds_stress_stab and fldpara_->IsNewton())
    {
      ReynoldsStressStab(estif_u_, estif_p_v_, lin_resM_Du_, timefacfac, timefacfacpre, fac3);
    }

    // 16) fine-scale subgrid-viscosity term
    //     (contribution only to right-hand-side vector)
    if (fldpara_->Fssgv() != INPAR::FLUID::no_fssgv)
    {
      const double fssgviscfac = fssgvisc_ * rhsfac;

      FineScaleSubGridViscosityTerm(velforce_, fssgviscfac);
    }

    // 17) subgrid-stress term (multifractal subgrid scales)
    if (fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales)
    {
      MultfracSubGridScalesCross(estif_u_, velforce_, timefacfac, rhsfac);

      MultfracSubGridScalesReynolds(estif_u_, velforce_, timefacfac, rhsfac);
    }

    // 18) polynomial pressure projection term (Dohrmann, Bochev IJNME 2004)
    //     (parameter-free inf-sub-stabilization, e.g. used instead of PSPG)
    if (fldpara_->PPP())
    {
      PressureProjection(ppmat_);
    }

    // linearization wrt mesh motion
    if (emesh.IsInitialized())
    {
      if (nsd_ == 3)
        LinMeshMotion_3D(emesh, evelaf, press, fldparatimint_->TimeFac(), timefacfac);
      else if (nsd_ == 2)
        LinMeshMotion_2D(emesh, evelaf, press, fldparatimint_->TimeFac(), timefacfac);
      else
        dserror("Linearization of the mesh motion is not available in 1D");
    }
  }
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  // if polynomial pressure projection: finalize matrices and rhs
  if (fldpara_->PPP())
  {
    if (fldparatimint_->IsGenalphaNP())
      PressureProjectionFinalize(ppmat_, preforce_, eprenp);
    else
      PressureProjectionFinalize(ppmat_, preforce_, epreaf);
  }

  //------------------------------------------------------------------------
  //  add contributions to element matrix and right-hand-side vector
  //------------------------------------------------------------------------
  // add pressure part to right-hand-side vector
  for (int vi = 0; vi < nen_; ++vi)
  {
    eforce(numdofpernode_ * vi + nsd_) += preforce_(vi);
  }

  // add velocity part to right-hand-side vector
  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      eforce(numdofpernode_ * vi + idim) += velforce_(idim, vi);
    }
  }

  // add pressure-pressure part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int fuippp = numdofpernode_ * ui + nsd_;

    for (int vi = 0; vi < nen_; ++vi)
    {
      const int numdof_vi_p_nsd = numdofpernode_ * vi + nsd_;

      estif(numdof_vi_p_nsd, fuippp) += ppmat_(vi, ui);
    }
  }

  // add velocity-velocity part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int numdof_ui = numdofpernode_ * ui;
    const int nsd_ui = nsd_ * ui;

    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      const int numdof_ui_jdim = numdof_ui + jdim;
      const int nsd_ui_jdim = nsd_ui + jdim;

      for (int vi = 0; vi < nen_; ++vi)
      {
        const int numdof_vi = numdofpernode_ * vi;
        const int nsd_vi = nsd_ * vi;

        for (int idim = 0; idim < nsd_; ++idim)
        {
          estif(numdof_vi + idim, numdof_ui_jdim) += estif_u_(nsd_vi + idim, nsd_ui_jdim);
        }
      }
    }
  }

  // add velocity-pressure part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int numdof_ui_nsd = numdofpernode_ * ui + nsd_;

    for (int vi = 0; vi < nen_; ++vi)
    {
      const int nsd_vi = nsd_ * vi;
      const int numdof_vi = numdofpernode_ * vi;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        estif(numdof_vi + idim, numdof_ui_nsd) += estif_p_v_(nsd_vi + idim, ui);
      }
    }
  }

  // add pressure-velocity part to matrix
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int numdof_ui = numdofpernode_ * ui;
    const int nsd_ui = nsd_ * ui;

    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      const int numdof_ui_jdim = numdof_ui + jdim;
      const int nsd_ui_jdim = nsd_ui + jdim;

      for (int vi = 0; vi < nen_; ++vi)
        estif(numdofpernode_ * vi + nsd_, numdof_ui_jdim) += estif_q_u_(vi, nsd_ui_jdim);
    }
  }

  return;
}

template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::CalcDivEps(
    const LINALG::Matrix<nsd_, nen_>& evelaf, const LINALG::Matrix<nsd_, nen_>& eveln)
{
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

  // set visc_old to zero
  visc_old_.Clear();
  visc_oldn_.Clear();

  double prefac;
  if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
  // if(loma_)
  {
    prefac = 1.0 / 3.0;
    derxy2_.Scale(prefac);
  }
  else
    prefac = 1.0;
  // reconstruction of second derivative via projection or superconvergent patch recovery
  if (fldpara_->IsReconstructDer())
  {
    if (is_higher_order_ele_ == false) dserror("this doesn't make sense");

    // global second derivatives of evalaf (projected velgrad)
    // not symmetric!
    LINALG::Matrix<nsd_ * nsd_, nsd_> evelgradderxy;
    evelgradderxy.MultiplyNT(evelafgrad_, derxy_);
    //*VELNP*
    if (nsd_ == 3)
    {
      // assemble div epsilon(evelaf)
      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int nsd_idim = idim * nsd_;
        // select diagonal entries (u,xx + u,yy + u,zz)
        double sum1 = (evelgradderxy(nsd_idim, 0) + evelgradderxy(nsd_idim + 1, 1) +
                          evelgradderxy(nsd_idim + 2, 2)) /
                      prefac;
        // interpolate mixed terms
        double sum2;
        if (idim == 0)
        {
          // uy,xy + uz,xz
          sum2 = 0.5 * (evelgradderxy(3, 1) + evelgradderxy(4, 0) + evelgradderxy(6, 2) +
                           evelgradderxy(8, 0));
        }
        else if (idim == 1)
        {
          // ux,xy + uz,yz
          sum2 = 0.5 * (evelgradderxy(1, 0) + evelgradderxy(0, 1) + evelgradderxy(7, 2) +
                           evelgradderxy(8, 1));
        }
        else
        {
          // ux,xz + uy,yz
          sum2 = 0.5 * (evelgradderxy(2, 0) + evelgradderxy(0, 2) + evelgradderxy(5, 1) +
                           evelgradderxy(4, 2));
        }
        // assemble each row of div epsilon(evelaf)
        visc_old_(idim) = 0.5 * (sum1 + evelgradderxy(nsd_idim + idim, idim) + sum2);
      }
    }
    else if (nsd_ == 2)
    {
      // assemble div epsilon(evelaf)
      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int nsd_idim = idim * nsd_;
        // select diagonal entries (u,xx + u,yy)
        double sum1 = (evelgradderxy(nsd_idim, 0) + evelgradderxy(nsd_idim + 1, 1)) / prefac;
        // interpolate mixed terms
        double sum2;
        if (idim == 0)
        {
          // uy,xy
          sum2 = 0.5 * (evelgradderxy(2, 1) + evelgradderxy(3, 0));
        }
        else
        {
          // ux,xy
          sum2 = 0.5 * (evelgradderxy(1, 0) + evelgradderxy(0, 1));
        }
        // assemble each row of div epsilon(evelaf)
        visc_old_(idim) = 0.5 * (sum1 + evelgradderxy(nsd_idim + idim, idim) + sum2);
      }
    }
    else
      dserror("Epsilon(N) is not implemented for the 1D case");

    //*VELN*
    // LINALG::Matrix<nsd_*nsd_,nsd_> evelngradderxy;
    evelgradderxy.Clear();
    evelgradderxy.MultiplyNT(evelngrad_, derxy_);
    if (nsd_ == 3)
    {
      // assemble div epsilon(evelaf)
      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int nsd_idim = idim * nsd_;
        // select diagonal entries (u,xx + u,yy + u,zz)
        double sum1 = (evelgradderxy(nsd_idim, 0) + evelgradderxy(nsd_idim + 1, 1) +
                          evelgradderxy(nsd_idim + 2, 2)) /
                      prefac;
        // interpolate mixed terms
        double sum2;
        if (idim == 0)
        {
          // uy,xy + uz,xz
          sum2 = 0.5 * (evelgradderxy(3, 1) + evelgradderxy(4, 0) + evelgradderxy(6, 2) +
                           evelgradderxy(8, 0));
        }
        else if (idim == 1)
        {
          // ux,xy + uz,yz
          sum2 = 0.5 * (evelgradderxy(1, 0) + evelgradderxy(0, 1) + evelgradderxy(7, 2) +
                           evelgradderxy(8, 1));
        }
        else
        {
          // ux,xz + uy,yz
          sum2 = 0.5 * (evelgradderxy(2, 0) + evelgradderxy(0, 2) + evelgradderxy(5, 1) +
                           evelgradderxy(4, 2));
        }
        // assemble each row of div epsilon(evelaf)
        visc_oldn_(idim) = 0.5 * (sum1 + evelgradderxy(nsd_idim + idim, idim) + sum2);
      }
    }
    else if (nsd_ == 2)
    {
      // assemble div epsilon(evelaf)
      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int nsd_idim = idim * nsd_;
        // select diagonal entries (u,xx + u,yy)
        double sum1 = (evelgradderxy(nsd_idim, 0) + evelgradderxy(nsd_idim + 1, 1)) / prefac;
        // interpolate mixed terms
        double sum2;
        if (idim == 0)
        {
          // uy,xy
          sum2 = 0.5 * (evelgradderxy(2, 1) + evelgradderxy(3, 0));
        }
        else
        {
          // ux,xy
          sum2 = 0.5 * (evelgradderxy(1, 0) + evelgradderxy(0, 1));
        }
        // assemble each row of div epsilon(evelaf)
        visc_oldn_(idim) = 0.5 * (sum1 + evelgradderxy(nsd_idim + idim, idim) + sum2);
      }
    }
    else
      dserror("Epsilon(N) is not implemented for the 1D case");
  }
  else
  {
    if (nsd_ == 3)
    {
      for (int inode = 0; inode < nen_; ++inode)
      {
        double sum = (derxy2_(0, inode) + derxy2_(1, inode) + derxy2_(2, inode)) / prefac;
        viscs2_(0, inode) = 0.5 * (sum + derxy2_(0, inode));
        viscs2_(1, inode) = 0.5 * derxy2_(3, inode);
        viscs2_(2, inode) = 0.5 * derxy2_(4, inode);
        viscs2_(3, inode) = 0.5 * derxy2_(3, inode);
        viscs2_(4, inode) = 0.5 * (sum + derxy2_(1, inode));
        viscs2_(5, inode) = 0.5 * derxy2_(5, inode);
        viscs2_(6, inode) = 0.5 * derxy2_(4, inode);
        viscs2_(7, inode) = 0.5 * derxy2_(5, inode);
        viscs2_(8, inode) = 0.5 * (sum + derxy2_(2, inode));

        for (int idim = 0; idim < nsd_; ++idim)
        {
          const int nsd_idim = idim * nsd_;
          for (int jdim = 0; jdim < nsd_; ++jdim)
          {
            visc_old_(idim) += viscs2_(nsd_idim + jdim, inode) * evelaf(jdim, inode);
            visc_oldn_(idim) += viscs2_(nsd_idim + jdim, inode) * eveln(jdim, inode);
          }
        }
      }
    }
    else if (nsd_ == 2)
    {
      for (int inode = 0; inode < nen_; ++inode)
      {
        double sum = (derxy2_(0, inode) + derxy2_(1, inode)) / prefac;
        viscs2_(0, inode) = 0.5 * (sum + derxy2_(0, inode));
        viscs2_(1, inode) = 0.5 * derxy2_(2, inode);
        viscs2_(2, inode) = 0.5 * derxy2_(2, inode);
        viscs2_(3, inode) = 0.5 * (sum + derxy2_(1, inode));

        for (int idim = 0; idim < nsd_; ++idim)
        {
          const int nsd_idim = idim * nsd_;
          for (int jdim = 0; jdim < nsd_; ++jdim)
          {
            visc_old_(idim) += viscs2_(nsd_idim + jdim, inode) * evelaf(jdim, inode);
            visc_oldn_(idim) += viscs2_(nsd_idim + jdim, inode) * eveln(jdim, inode);
          }
        }
      }
    }
    else
      dserror("Epsilon(N) is not implemented for the 1D case");
  }

  return;
}


template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::ComputeSubgridScaleVelocityOSTNew(
    const LINALG::Matrix<nsd_, nen_>& eaccam, double& fac1, double& fac2, double& fac3,
    double& facMtau, int iquad, double* saccn, double* sveln, double* svelnp)
{
  //----------------------------------------------------------------------
  // compute residual of momentum equation
  // -> different for generalized-alpha and other time-integration schemes
  //----------------------------------------------------------------------
  if (fldparatimint_->IsGenalpha())
  {
    if (fldpara_->PhysicalType() == INPAR::FLUID::boussinesq)
      dserror(
          "The combination of generalized-alpha time integration and a Boussinesq approximation "
          "has not been implemented yet!");

    // rhs of momentum equation: density*bodyforce at n+alpha_F
    rhsmom_.Update(densaf_, bodyforce_, 0.0);
    // and pressure gradient prescribed as body force
    // caution: not density weighted
    rhsmom_.Update(1.0, generalbodyforce_, 1.0);

    // get acceleration at time n+alpha_M at integration point
    accint_.Multiply(eaccam, funct_);

    // evaluate momentum residual once for all stabilization right hand sides
    for (int rr = 0; rr < nsd_; ++rr)
    {
      momres_old_(rr) = densam_ * accint_(rr) + densaf_ * conv_old_(rr) + gradp_(rr) -
                        2 * visceff_ * visc_old_(rr) + reacoeff_ * velint_(rr) -
                        densaf_ * bodyforce_(rr) - generalbodyforce_(rr);
    }

    // add consistency terms for MFS if applicable
    MultfracSubGridScalesConsistentResidual();
  }
  else
  {
    if (not fldparatimint_->IsStationary())
    {
      // rhs of instationary momentum equation:
      // density*theta*bodyforce at n+1 + density*(histmom/dt)
      // in the case of a Boussinesq approximation: f = rho_0*[(rho - rho_0)/rho_0]*g = (rho -
      // rho_0)*g else:                                      f = rho * g
      if (fldpara_->PhysicalType() == INPAR::FLUID::boussinesq)
      {
        //      Made old OST impl equivalent to gen-alpha (alpha_f=alpha_m=1) (multiplied with
        //      \rho_(n+1))
        rhsmom_.Update((densaf_ / fldparatimint_->Dt() / fldparatimint_->Theta()), histmom_,
            deltadens_, bodyforce_);
        // and pressure gradient prescribed as body force
        // caution: not density weighted
        rhsmom_.Update(1.0, generalbodyforce_, 1.0);
      }
      else
      {
        //      Made old OST impl equivalent to gen-alpha (alpha_f=alpha_m=1) (multiplied with
        //      \rho_(n+1))
        rhsmom_.Update((densaf_ / fldparatimint_->Dt() / fldparatimint_->Theta()), histmom_,
            densaf_, bodyforce_);

        // and pressure gradient prescribed as body force
        // caution: not density weighted
        rhsmom_.Update(1.0, generalbodyforce_, 1.0);
      }

      // compute instationary momentum residual:
      // momres_old = u_(n+1)/dt + theta ( ... ) +(1-theta) ( ... ) - theta*bodyforce_
      if (fldparatimint_->IsNewOSTImplementation())
      {
        const double quotfac =
            (1.0 - fldparatimint_->Theta()) * fldparatimint_->Dt() / fldparatimint_->TimeFacRhs();
        for (int rr = 0; rr < nsd_; ++rr)
        {
          momres_old_(rr) = dens_theta_ * (velint_(rr) - velintn_(rr)) /
                            (fldparatimint_->Dt() * fldparatimint_->Theta());
          momres_old_(rr) -= (densaf_ * bodyforce_(rr) + densn_ * quotfac * bodyforcen_(rr));
          momres_old_(rr) +=
              reacoeff_ * (velint_(rr) + quotfac * velintn_(rr));  // TODO: Time dependant reacoef.
          momres_old_(rr) += (densaf_ * conv_old_(rr) + densn_ * quotfac * conv_oldn_(rr));
          momres_old_(rr) -= 2.0 * (visceff_ * visc_old_(rr) + viscn_ * quotfac * visc_oldn_(rr));
          momres_old_(rr) -= generalbodyforce_(rr) + quotfac * generalbodyforcen_(rr);
          if (not fldparatimint_->IsFullImplPressureAndCont())
          {
            momres_old_(rr) += gradp_(rr) + quotfac * gradpn_(rr);
          }
          else
          {
            momres_old_(rr) +=
                gradp_(rr) / fldparatimint_->Theta();  // Gradient of p with no pre-factor.
          }

          // Left for implementation in ALE, possibly....
          //        //Furthermore, should rhsmom_ be calculated like this? Not with galerkin terms
          //        and integration???! rhsmom_(rr)=- fldparatimint_->Dt() *fldparatimint_->Theta()*
          //        (- densaf_*velintn_(rr)/(fldparatimint_->Dt()*fldparatimint_->Theta())
          //                                          - densn_*quotfac*bodyforcen_(rr)
          //                                          + reacoeff_*quotfac*velintn_(rr) +
          //                                          densaf_*quotfac*conv_oldn_(rr)
          //                                          - 2*viscn_*quotfac*visc_oldn_(rr)  +
          //                                          quotfac*generalbodyforcen_(rr)
          //                                          + quotfac*gradpn_(rr)); //Check what gradpn_
          //                                          to use! (This is needed for ALE
          //                                          implementations?)
        }
      }
      else
      {
        // momres_old = u_(n+1)/dt + theta ( ... ) - histmom_/dt - theta*bodyforce_
        for (int rr = 0; rr < nsd_; ++rr)
        {
          momres_old_(rr) = ((densaf_ * velint_(rr) / fldparatimint_->Dt() +
                                 fldparatimint_->Theta() *
                                     (densaf_ * conv_old_(rr) + gradp_(rr) -
                                         2 * visceff_ * visc_old_(rr) + reacoeff_ * velint_(rr))) /
                                fldparatimint_->Theta()) -
                            rhsmom_(rr);
        }
      }  // end IsNewOSTImplementation
    }
    else
    {
      // rhs of stationary momentum equation: density*bodyforce
      // in the case of a Boussinesq approximation: f = rho_0*[(rho - rho_0)/rho_0]*g = (rho -
      // rho_0)*g else:                                      f = rho * g and pressure gradient
      // prescribed as body force (not density weighted)
      if (fldpara_->PhysicalType() == INPAR::FLUID::boussinesq)
        rhsmom_.Update(deltadens_, bodyforce_, 1.0, generalbodyforce_);
      else
        rhsmom_.Update(densaf_, bodyforce_, 1.0, generalbodyforce_);

      // compute stationary momentum residual:
      for (int rr = 0; rr < nsd_; ++rr)
      {
        momres_old_(rr) = -rhsmom_(rr);
        momres_old_(rr) += gradp_(rr);
        momres_old_(rr) += reacoeff_ * velint_(rr);
        momres_old_(rr) += densaf_ * conv_old_(rr);
        momres_old_(rr) += -2 * visceff_ * visc_old_(rr);
      }

      // add consistency terms for MFS if applicable
      MultfracSubGridScalesConsistentResidual();
    }
  }

  //----------------------------------------------------------------------
  // compute subgrid-scale velocity
  //----------------------------------------------------------------------
  // 1) quasi-static subgrid scales
  // Definition of subgrid-scale velocity is not consistent for the SUPG term and Franca, Valentin,
  // ... Definition of subgrid velocity used by Hughes
  if (fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
  {
    sgvelint_.Update(-tau_(1), momres_old_, 0.0);
  }
  // 2) time-dependent subgrid scales
  else
  {
    // some checking
    if (fldparatimint_->IsStationary())
      dserror("there is no time dependent subgrid scale closure for stationary problems\n");
    if (saccn == NULL or sveln == NULL or svelnp == NULL) dserror("no subscale array provided");

    // parameter definitions
    double alphaF = fldparatimint_->AlphaF();
    double alphaM = fldparatimint_->AlphaM();
    double gamma = fldparatimint_->Gamma();
    double dt = fldparatimint_->Dt();

    /*
                                            1.0
       facMtau =  -------------------------------------------------------
                     n+aM                      n+aF
                  rho     * alphaM * tauM + rho     * alphaF * gamma * dt
    */
    facMtau = 1.0 / (densam_ * alphaM * tau_(1) + densaf_ * fldparatimint_->Afgdt());

    /*
       factor for old subgrid velocities:

                 n+aM                      n+aF
       fac1 = rho     * alphaM * tauM + rho     * gamma * dt * (alphaF-1)
    */
    fac1 = (densam_ * alphaM * tau_(1) + densaf_ * gamma * dt * (alphaF - 1.0)) * facMtau;
    /*
      factor for old subgrid accelerations

                 n+aM
       fac2 = rho     * tauM * dt * (alphaM-gamma)
    */
    fac2 = (densam_ * dt * tau_(1) * (alphaM - gamma)) * facMtau;
    /*
      factor for residual in current subgrid velocities:

       fac3 = gamma * dt * tauM
    */
    fac3 = (gamma * dt * tau_(1)) * facMtau;

    // warning: time-dependent subgrid closure requires generalized-alpha time
    // integration
    if (!fldparatimint_->IsGenalpha())
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

    static LINALG::Matrix<1, nsd_> sgvelintaf(true);
    sgvelintaf.Clear();
    for (int rr = 0; rr < nsd_; ++rr)
    {
      tds_->UpdateSvelnpInOneDirection(fac1, fac2, fac3, momres_old_(rr), fldparatimint_->AlphaF(),
          rr, iquad,
          sgvelint_(rr),  // sgvelint_ is set to sgvelintnp, but is then overwritten below anyway!
          sgvelintaf(rr));

      int pos = rr + nsd_ * iquad;

      /*
       *  ~n+1           ~n           ~ n            n+1
       *  u    =  fac1 * u  + fac2 * acc  -fac3 * res
       *   (i)
       *
       */

      svelnp[pos] = fac1 * sveln[pos] + fac2 * saccn[pos] - fac3 * momres_old_(rr);

      /* compute the intermediate value of subscale velocity
       *
       *          ~n+af            ~n+1                   ~n
       *          u     = alphaF * u     + (1.0-alphaF) * u
       *           (i)              (i)
       *
       */
      sgvelint_(rr) = alphaF * svelnp[pos] + (1.0 - alphaF) * sveln[pos];
    }
  }  // end time dependent subgrid scale closure

  //----------------------------------------------------------------------
  // include computed subgrid-scale velocity in convective term
  // -> only required for cross- and Reynolds-stress terms
  //----------------------------------------------------------------------
  if (fldpara_->Cross() != INPAR::FLUID::cross_stress_stab_none or
      fldpara_->Reynolds() != INPAR::FLUID::reynolds_stress_stab_none or
      fldpara_->ContiCross() != INPAR::FLUID::cross_stress_stab_none or
      fldpara_->ContiReynolds() != INPAR::FLUID::reynolds_stress_stab_none)
    sgconv_c_.MultiplyTN(derxy_, sgvelint_);
  else
    sgconv_c_.Clear();
}


template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::LinGalMomResUOSTNew(
    LINALG::Matrix<nsd_ * nsd_, nen_>& lin_resM_Du, const double& timefacfac)
{
  /*
      instationary                          cross-stress, part 1
       +-----+                             +-------------------+
       |     |                             |                   |

                 /       n+1       \        /      ~n+1       \
       rho*Du + |   rho*u   o nabla | Du + |   rho*u   o nabla | Du +
                 \      (i)        /        \      (i)        /

                 /                \  n+1
              + |   rho*Du o nabla | u      +  sigma*Du
                 \                /   (i)
                |                        |     |       |
                +------------------------+     +-------+
                        Newton                  reaction
  */

  int idim_nsd_p_idim[nsd_];

  for (int idim = 0; idim < nsd_; ++idim)
  {
    idim_nsd_p_idim[idim] = idim * nsd_ + idim;
  }

  if (fldparatimint_->IsStationary() == false)
  {
    //    double fac_densam= (fldparatimint_->IsOneStepTheta()) ? fac_*dens_theta_ : fac_*densam_;
    double fac_densam;
    if (fldparatimint_->IsNewOSTImplementation())
      fac_densam = fldparatimint_->IsOneStepTheta() ? fac_ * dens_theta_ : fac_ * densam_;
    else
      fac_densam = fac_ * densam_;
    // End of IsNewOSTImplementation()

    for (int ui = 0; ui < nen_; ++ui)
    {
      const double v = fac_densam * funct_(ui);

      for (int idim = 0; idim < nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim], ui) += v;
      }
    }
  }

  const double timefacfac_densaf = timefacfac * densaf_;

  // convection, reactive
  for (int ui = 0; ui < nen_; ++ui)
  {
    const double v = timefacfac_densaf * conv_c_(ui);

    for (int idim = 0; idim < nsd_; ++idim)
    {
      lin_resM_Du(idim_nsd_p_idim[idim], ui) += v;
    }
  }

  // convection, convective (only for Newton)
  if (fldpara_->IsNewton())
  {
    for (int ui = 0; ui < nen_; ++ui)
    {
      const double temp = timefacfac_densaf * funct_(ui);

      for (int idim = 0; idim < nsd_; ++idim)
      {
        const int idim_nsd = idim * nsd_;

        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          lin_resM_Du(idim_nsd + jdim, ui) += temp * vderxy_(idim, jdim);
        }
      }
    }
  }

  if (fldpara_->Reaction())
  {
    const double fac_reac = timefacfac * reacoeff_;

    for (int ui = 0; ui < nen_; ++ui)
    {
      const double v = fac_reac * funct_(ui);

      for (int idim = 0; idim < nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim], ui) += v;
      }
    }
  }

  if (fldpara_->Cross() == INPAR::FLUID::cross_stress_stab)
  {
    // const double rhsresfac_densaf=rhsresfac*densaf_;
    for (int ui = 0; ui < nen_; ++ui)
    {
      // const double v=rhsresfac_densaf*sgconv_c_(ui);
      const double v = timefacfac_densaf * sgconv_c_(ui);

      for (int idim = 0; idim < nsd_; ++idim)
      {
        lin_resM_Du(idim_nsd_p_idim[idim], ui) += v;
      }
    }
  }

  return;
}

template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::InertiaConvectionReactionGalPart(
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u, LINALG::Matrix<nsd_, nen_>& velforce,
    LINALG::Matrix<nsd_ * nsd_, nen_>& lin_resM_Du, LINALG::Matrix<nsd_, 1>& resM_Du,
    const double& rhsfac, const double& rhsfacn)
{
  /* inertia (contribution to mass matrix) if not is_stationary */
  /*
            /              \
           |                |
           |    rho*Du , v  |
           |                |
            \              /
  */
  /* convection, convective part (convective form) */
  /*
            /                             \
           |  /       n+1       \          |
           | |   rho*u   o nabla | Du , v  |
           |  \      (i)        /          |
            \                             /
  */
  /*  convection, reactive part (convective form)
            /                               \
           |  /                \   n+1       |
           | |  rho*Du o nabla  | u     , v  |
           |  \                /   (i)       |
            \                               /
  */
  /*  reaction */
  /*
            /                \
           |                  |
           |    sigma*Du , v  |
           |                  |
            \                /
  */
  if ((fldpara_->IsNewton() or
          (is_higher_order_ele_ and fldpara_->Tds() == INPAR::FLUID::subscales_time_dependent)))
  {
    for (int ui = 0; ui < nen_; ++ui)
    {
      for (int vi = 0; vi < nen_; ++vi)
      {
        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          for (int idim = 0; idim < nsd_; ++idim)
          {
            estif_u(nsd_ * vi + idim, nsd_ * ui + jdim) +=
                funct_(vi) * lin_resM_Du(idim * nsd_ + jdim, ui);
          }  // end for (idim)
        }    // end for (jdim)
      }      // end for (vi)
    }        // end for (ui)
  }
  else
  {
    for (int ui = 0; ui < nen_; ++ui)
    {
      for (int vi = 0; vi < nen_; ++vi)
      {
        for (int idim = 0; idim < nsd_; ++idim)
        {
          estif_u(nsd_ * vi + idim, nsd_ * ui + idim) +=
              funct_(vi) * lin_resM_Du(idim * nsd_ + idim, ui);
        }  // end for (idim)
      }    // vi
    }      // ui
  }

  // inertia terms on the right hand side for instationary fluids
  if (not fldparatimint_->IsStationary())
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      if (fldparatimint_->IsGenalpha())
      {
        if (fldpara_->Tds() == INPAR::FLUID::subscales_time_dependent &&
            fldpara_->Transient() == INPAR::FLUID::inertia_stab_keep)
        {
          ;  // do nothing here! Whole term already set in LinGalMomResU_subscales()
        }
        else
          resM_Du(idim) += rhsfac * densam_ * accint_(idim);
      }
      else
      {
        if (fldparatimint_->IsNewOSTImplementation())
        {
          // this approximates \int_{t_n}^{t_{n+1}} \rho du/dt
          // It could be implemented differently like (through integration by parts):
          // \rho_{n+1}u_{n+1}-\rho_{n}u_{n} -  \int_{t_n}^{t_{n+1}} d \rho/dt u
          // But leaves the second integral to be approximated.
          resM_Du(idim) += fac_ * dens_theta_ * (velint_(idim) - velintn_(idim));
        }
        else
        {
          resM_Du(idim) += fac_ * densaf_ * velint_(idim);
        }  // end IsNewOSTImplementation
      }
    }
  }  // end if (not stationary)

  // convective terms of rhs
  for (int idim = 0; idim < nsd_; ++idim)
  {
    if (fldpara_->Tds() == INPAR::FLUID::subscales_time_dependent &&
        fldpara_->Transient() == INPAR::FLUID::inertia_stab_keep)
    {
      ;  // do nothing here! Whole term already set in LinGalMomResU_subscales()
    }
    else
    {
      resM_Du(idim) += rhsfac * densaf_ * conv_old_(idim);
      if (fldparatimint_->IsNewOSTImplementation())
      {
        if (fldparatimint_->IsOneStepTheta())
        {
          resM_Du(idim) += rhsfacn * densn_ * conv_oldn_(idim);
        }
      }  // end IsNewOSTImplementation
    }
  }  // end for(idim)

  if (fldpara_->Reaction())
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      resM_Du(idim) += rhsfac * reacoeff_ * velint_(idim);
      if (fldparatimint_->IsNewOSTImplementation())
      {
        if (fldparatimint_->IsOneStepTheta())
        {
          resM_Du(idim) += rhsfacn * reacoeff_ * velintn_(idim);
        }
      }  // end IsNewOSTImplementation
    }
  }  // end if (reaction_)

  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      velforce(idim, vi) -= resM_Du(idim) * funct_(vi);
    }
  }
  return;
}


template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::ViscousGalPart(
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u, LINALG::Matrix<nsd_, nen_>& velforce,
    LINALG::Matrix<nsd_, nsd_>& viscstress, const double& timefacfac, const double& rhsfac,
    const double& rhsfacn)
{
  const double visceff_timefacfac = visceff_ * timefacfac;

  /* viscosity term */
  /*
                   /                        \
                  |       /  \         / \   |
            2 mu  |  eps | Du | , eps | v |  |
                  |       \  /         \ /   |
                   \                        /
  */

  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      const double temp = visceff_timefacfac * derxy_(jdim, vi);

      for (int ui = 0; ui < nen_; ++ui)
      {
        for (int idim = 0; idim < nsd_; ++idim)
        {
          estif_u(nsd_ * vi + idim, nsd_ * ui + jdim) += temp * derxy_(idim, ui);
        }
      }
    }
  }

  static LINALG::Matrix<nen_, nen_> tmp_dyad;
  tmp_dyad.MultiplyTN(derxy_, derxy_);
  tmp_dyad.Scale(visceff_timefacfac);

  for (int ui = 0; ui < nen_; ++ui)
  {
    for (int vi = 0; vi < nen_; ++vi)
    {
      const double tmp_val = tmp_dyad(vi, ui);

      for (int idim = 0; idim < nsd_; ++idim)
      {
        estif_u(nsd_ * vi + idim, nsd_ * ui + idim) += tmp_val;
      }  // end for (idim)
    }    // ui
  }      // vi

  static LINALG::Matrix<nsd_, nsd_> viscstressn;

  const double v = visceff_ * rhsfac;
  const double vn = viscn_ * rhsfacn;

  for (int jdim = 0; jdim < nsd_; ++jdim)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      viscstress(idim, jdim) = v * (vderxy_(jdim, idim) + vderxy_(idim, jdim));
      viscstressn(idim, jdim) = vn * (vderxyn_(jdim, idim) + vderxyn_(idim, jdim));
    }
  }


  static LINALG::Matrix<nsd_, nen_> tmp;


  if (fldparatimint_->IsNewOSTImplementation())
  {
    if (fldparatimint_->IsOneStepTheta())
    {
      static LINALG::Matrix<nsd_, nsd_> viscstress_added;

      viscstress_added.Update(1.0, viscstress, 1.0, viscstressn, 0.0);
      tmp.Multiply(viscstress_added, derxy_);
    }
  }
  else
    tmp.Multiply(viscstress, derxy_);

  velforce.Update(-1.0, tmp, 1.0);

  return;
}

template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::ContStab(
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u, LINALG::Matrix<nsd_, nen_>& velforce,
    const double& timefac, const double& timefacfac, const double& timefacfacpre,
    const double& rhsfac, const double& rhsfacn)
{
  // In the case no continuity stabilization and no LOMA:
  // the factors 'conti_stab_and_vol_visc_fac' and 'conti_stab_and_vol_visc_rhs' are zero
  // therefore there is no contribution to the element stiffness matrix and
  // the viscous stress tensor is NOT altered!!
  //
  // ONLY
  // the rhs contribution of the viscous term is added!!

  double conti_stab_and_vol_visc_fac = 0.0;
  double conti_stab_and_vol_visc_rhs = 0.0;

  if (fldpara_->CStab())
  {
    if (not fldparatimint_->IsFullImplPressureAndCont())
    {
      conti_stab_and_vol_visc_fac += timefacfacpre * tau_(2);
      conti_stab_and_vol_visc_rhs -= rhsfac * tau_(2) * conres_old_;
      if (fldparatimint_->IsNewOSTImplementation())
      {
        if (not fldparatimint_->IsImplPressure())
        {
          if (fldparatimint_->IsOneStepTheta())
          {
            conti_stab_and_vol_visc_rhs -= rhsfacn * tau_(2) * conres_oldn_;
          }
        }
      }  // end IsNewOSTImplementation
    }
    else
    {
      // Full impl pressure weighted with \Delta t only.
      conti_stab_and_vol_visc_fac += fac_ * fldparatimint_->Dt() * tau_(2);
      conti_stab_and_vol_visc_rhs -= fac_ * fldparatimint_->Dt() * tau_(2) * conres_old_;
    }
  }
  if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
  {
    conti_stab_and_vol_visc_fac -= (2.0 / 3.0) * visceff_ * timefacfac;
    conti_stab_and_vol_visc_rhs += (2.0 / 3.0) * rhsfac * visceff_ * vdiv_;
    // additional term q_sq_ for dynamics Smagorisnky only
    if (fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky or
        fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or
        fldpara_->TurbModAction() == INPAR::FLUID::vreman)
      conti_stab_and_vol_visc_rhs += (2.0 / 3.0) * rhsfac * q_sq_;
  }

  /* continuity stabilisation on left-hand side */
  /*
              /                        \
             |                          |
        tauC | nabla o Du  , nabla o v  |
             |                          |
              \                        /
  */
  /* viscosity term - subtraction for low-Mach-number flow */
  /*
             /                             \             /                        \
            |  1                      / \   |     2 mu  |                          |
     - 2 mu |  - (nabla o u) I , eps | v |  | = - ----- | nabla o Du  , nabla o v  |
            |  3                      \ /   |       3   |                          |
             \                             /             \                        /
  */
  for (int ui = 0; ui < nen_; ++ui)
  {
    const int fui = nsd_ * ui;

    for (int idim = 0; idim < nsd_; ++idim)
    {
      const int fui_p_idim = fui + idim;
      const double v0 = conti_stab_and_vol_visc_fac * derxy_(idim, ui);
      for (int vi = 0; vi < nen_; ++vi)
      {
        const int fvi = nsd_ * vi;

        for (int jdim = 0; jdim < nsd_; ++jdim)
        {
          estif_u(fvi + jdim, fui_p_idim) += v0 * derxy_(jdim, vi);
        }
      }
    }  // end for(idim)
  }

  // computation of right-hand-side viscosity term
  for (int vi = 0; vi < nen_; ++vi)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      /* viscosity term on right-hand side */
      velforce(idim, vi) += conti_stab_and_vol_visc_rhs * derxy_(idim, vi);
    }
  }

  return;
}


template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::PressureGalPart(
    LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v, LINALG::Matrix<nsd_, nen_>& velforce,
    const double& timefacfac, const double& timefacfacpre, const double& rhsfac,
    const double& rhsfacn, const double& press, const double& pressn)
{
  for (int ui = 0; ui < nen_; ++ui)
  {
    double v;
    if (not fldparatimint_->IsFullImplPressureAndCont())
    {
      v = -timefacfacpre * funct_(ui);
    }
    else
    {
      v = -fldparatimint_->Dt() * fac_ * funct_(ui);
    }

    for (int vi = 0; vi < nen_; ++vi)
    {
      const int fvi = nsd_ * vi;
      /* pressure term */
      /*
           /                \
          |                  |
          |  Dp , nabla o v  |
          |                  |
           \                /
      */
      for (int idim = 0; idim < nsd_; ++idim)
      {
        estif_p_v(fvi + idim, ui) += v * derxy_(idim, vi);
      }
    }
  }

  // pressure term on right-hand side
  if (not fldparatimint_->IsFullImplPressureAndCont())
  {
    velforce.Update(press * rhsfac, derxy_, 1.0);
    if (fldparatimint_->IsNewOSTImplementation())
    {
      if (fldparatimint_->IsOneStepTheta())
      {
        velforce.Update(pressn * rhsfacn, derxy_, 1.0);
      }
    }  // end IsNewOSTImplementation
  }
  else
  {
    //     Full impl pressure weighted with \Delta t.
    velforce.Update(fac_ * fldparatimint_->Dt() * press, derxy_, 1.0);
  }
  return;
}



template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::ContinuityGalPart(
    LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u, LINALG::Matrix<nen_, 1>& preforce,
    const double& timefacfac, const double& timefacfacpre, const double& rhsfac,
    const double& rhsfacn)
{
  for (int vi = 0; vi < nen_; ++vi)
  {
    double v;
    if (not fldparatimint_->IsFullImplPressureAndCont())
    {
      v = timefacfacpre * funct_(vi);
    }
    else
    {
      v = fac_ * fldparatimint_->Dt() * funct_(vi);
    }

    for (int ui = 0; ui < nen_; ++ui)
    {
      const int fui = nsd_ * ui;

      for (int idim = 0; idim < nsd_; ++idim)
      {
        /* continuity term */
        /*
             /                \
            |                  |
            | nabla o Du  , q  |
            |                  |
             \                /
        */
        estif_q_u(vi, fui + idim) += v * derxy_(idim, ui);
      }
    }
  }  // end for(idim)
  if (not fldparatimint_->IsFullImplPressureAndCont())
  {
    preforce.Update(-rhsfac * vdiv_, funct_, 1.0);
    if (fldparatimint_->IsNewOSTImplementation())
    {
      if (not fldparatimint_->IsImplPressure())
      {
        if (fldparatimint_->IsOneStepTheta())
        {
          preforce.Update(-rhsfacn * vdivn_, funct_, 1.0);
        }
      }
    }  // end IsNewOSTImplementation
  }
  else
  {
    //     Full impl pressure weighted with \Delta t.
    preforce.Update(-fac_ * fldparatimint_->Dt() * vdiv_, funct_, 1.0);
  }
  return;
}

template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::BodyForceRhsTerm(
    LINALG::Matrix<nsd_, nen_>& velforce, const double& rhsfac, const double rhsfacn)
{
  for (int idim = 0; idim < nsd_; ++idim)
  {
    double scaled_rhsmom;
    if (fldparatimint_->IsGenalpha())
      scaled_rhsmom = rhsfac * rhsmom_(idim);
    else
    {
      if (fldparatimint_->IsNewOSTImplementation())
      {
        scaled_rhsmom = rhsfac * (densaf_ * bodyforce_(idim) + generalbodyforce_(idim));
        if (fldparatimint_->IsOneStepTheta())
        {
          scaled_rhsmom += rhsfacn * (densn_ * bodyforcen_(idim) + generalbodyforcen_(idim));
        }
      }
      else
      {
        scaled_rhsmom = rhsfac * rhsmom_(idim);
      }  // end IsNewOSTImplementation
    }

    for (int vi = 0; vi < nen_; ++vi)
    {
      velforce(idim, vi) += scaled_rhsmom * funct_(vi);
    }
  }  // end for(idim)

  return;
}

template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::PSPGOSTNew(
    LINALG::Matrix<nen_, nen_ * nsd_>& estif_q_u, LINALG::Matrix<nen_, nen_>& ppmat,
    LINALG::Matrix<nen_, 1>& preforce, LINALG::Matrix<nsd_ * nsd_, nen_>& lin_resM_Du,
    const double& fac3, const double& timefacfac, const double& timefacfacpre, const double& rhsfac,
    const int iquad)
{
  // conservative, stabilization terms are neglected (Hughes)

  /* pressure stabilisation:                                            */
  /*
              /                 \
             |  ~n+af            |
           - |  u     , nabla q  |
             |                   |
              \                 /
  */

  double scal_grad_q = 0.0;

  if (fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
  {
    scal_grad_q = tau_(1);
  }
  else  // time-dependent subgrid-scales
  {
    scal_grad_q = fac3;
  }

  /* pressure stabilisation: inertia if not stationary*/
  /*
            /                  \
           |                    |
           |  rho*Du , nabla q  |
           |                    |
            \                  /
  */
  /* pressure stabilisation: convection, convective part */
  /*
            /                                   \
           |  /       n+1       \                |
           | |   rho*u   o nabla | Du , nabla q  |
           |  \      (i)        /                |
            \                                   /
  */
  /* pressure stabilisation: convection, reactive part if Newton */
  /*
            /                                   \
           |  /                \   n+1           |
           | |   rho*Du o nabla | u     , grad q |
           |  \                /   (i)           |
            \                                   /
  */
  /* pressure stabilisation: reaction if included */
  /*
            /                     \
           |                      |
           |  sigma*Du , nabla q  |
           |                      |
            \                    /
  */
  /* pressure stabilisation: viscosity (-L_visc_u) */
  /*
            /                              \
           |               /  \             |
       mu  |  nabla o eps | Du | , nabla q  |
           |               \  /             |
            \                              /
  */

  if (is_higher_order_ele_ || fldpara_->IsNewton())
  {
    for (int ui = 0; ui < nen_; ++ui)
    {
      for (int vi = 0; vi < nen_; ++vi)
      {
        for (int idim = 0; idim < nsd_; ++idim)
        {
          const double temp_vi_idim = derxy_(idim, vi) * scal_grad_q;
          for (int jdim = 0; jdim < nsd_; ++jdim)
          {
            estif_q_u(vi, nsd_ * ui + jdim) += lin_resM_Du(nsd_ * idim + jdim, ui) * temp_vi_idim;
          }  // jdim
        }    // idim
      }      // vi
    }        // ui
  }          // end if (is_higher_order_ele_) or (newton_)
  else
  {
    for (int vi = 0; vi < nen_; ++vi)
    {
      for (int ui = 0; ui < nen_; ++ui)
      {
        for (int idim = 0; idim < nsd_; ++idim)
        {
          estif_q_u(vi, nsd_ * ui + idim) +=
              lin_resM_Du(nsd_ * idim + idim, ui) * derxy_(idim, vi) * scal_grad_q;
        }  // vi
      }    // ui
    }      // idim
  }        // end if not (is_higher_order_ele_) nor (newton_)


  for (int ui = 0; ui < nen_; ++ui)
  {
    /* pressure stabilisation: pressure( L_pres_p) */
    /*
         /                    \
        |                      |
        |  nabla Dp , nabla q  |
        |                      |
         \                    /
    */
    for (int vi = 0; vi < nen_; ++vi)
    {
      double sum = 0.;
      for (int idim = 0; idim < nsd_; ++idim) sum += derxy_(idim, ui) * derxy_(idim, vi);

      if (not fldparatimint_->IsFullImplPressureAndCont())
        ppmat(vi, ui) += timefacfacpre * scal_grad_q * sum;
      else
      {
        // Weighted with \Delta t
        ppmat(vi, ui) += fac_ * fldparatimint_->Dt() * scal_grad_q * sum;
      }
    }  // vi
  }    // ui

  for (int idim = 0; idim < nsd_; ++idim)
  {
    double sgvel = 0.0;
    if (fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
    {
      sgvel = sgvelint_(idim);
    }
    else  // time-dependent subgrid-scales, Np_Genal_Alpha!
    {
      sgvel = (tds_->Svelnp())(idim, iquad);
    }
    const double temp = rhsfac * sgvel;

    for (int vi = 0; vi < nen_; ++vi)
    {
      // pressure stabilisation
      preforce(vi) += temp * derxy_(idim, vi);
    }
  }  // end for(idim)

  return;
}

template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::SUPGOSTNew(
    LINALG::Matrix<nen_ * nsd_, nen_ * nsd_>& estif_u, LINALG::Matrix<nen_ * nsd_, nen_>& estif_p_v,
    LINALG::Matrix<nsd_, nen_>& velforce, LINALG::Matrix<nen_, 1>& preforce,
    LINALG::Matrix<nsd_ * nsd_, nen_>& lin_resM_Du, const double& fac3, const double& timefacfac,
    const double& timefacfacpre, const double& rhsfac)
{
  /*
                    /                                \
                   |  ~n+af    /     n+af       \     |
                 - |  u     , | rho*u    o nabla | v  |
                   |           \     (i)        /     |
                    \                                /
   */

  static LINALG::Matrix<nsd_, 1> temp;

  double supgfac;
  if (fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
    supgfac = densaf_ * tau_(0);
  else
    supgfac = densaf_ * fldparatimint_->AlphaF() * fac3;

  static LINALG::Matrix<nen_, 1> supg_test;
  for (int vi = 0; vi < nen_; ++vi)
  {
    supg_test(vi) = supgfac * conv_c_(vi);
  }

  if (fldpara_->Reynolds() == INPAR::FLUID::reynolds_stress_stab)
  {
    for (int vi = 0; vi < nen_; ++vi)
    {
      supg_test(vi) += supgfac * sgconv_c_(vi);
    }
  }

  /* supg stabilisation: inertia if not stationary */
  /*
         /                                \
        |            /     n+1       \     |
        |  rho*Du , | rho*u   o nabla | v  |
        |            \     (i)       /     |
         \                                /
  */
  /* supg stabilisation: convective part ( L_conv_u) , convective term */
  /*
         /                                                     \
        |    /       n+1        \        /      n+1       \     |
        |   |   rho*u    o nabla | Du , | rho*u    o nabla | v  |
        |    \       (i)        /        \      (i)       /     |
         \                                                     /
  */
  /* supg stabilisation: convective part ( L_conv_u) , reactive term if Newton */
  /*
         /                                                     \
        |    /       n+1        \        /     n+1        \     |
        |   |   rho*u    o nabla | Du , | rho*u    o nabla | v  |
        |    \       (i)        /        \     (i)        /     |
         \                                                     /
  */
  /* supg stabilisation: reaction if included */
  /*
         /                                  \
        |              /     n+1       \     |
        |  sigma*Du , | rho*u   o nabla | v  |
        |              \     (i)       /     |
         \                                  /
  */
  /* supg stabilisation: viscous part  (-L_visc_u) if is_higher_order_ele_ */
  /*
         /                                              \
        |               /  \    /       n+1        \     |
        |  nabla o eps | Du |, |   rho*u    o nabla | v  |
        |               \  /    \       (i)        /     |
         \                                              /
  */

  /* supg stabilisation: inertia, linearisation of testfunction if Newton */
  /*
              /                                       \
             |         n+1       /              \      |
             |    rho*u      ,  | rho*Du o nabla | v   |
             |         (i)       \              /      |
              \                                       /
  */
  /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u)
   * if Newton */
  /*
              /                                                       \
             |    /       n+1        \   n+1     /              \      |
             |   |   rho*u    o nabla | u    ,  | rho*Du o nabla | v   |
             |    \       (i)        /   (i)     \              /      |
              \                                                       /
  */
  /* supg stabilisation: reaction, linearisation of testfunction if Newton */
  /*
              /                                         \
             |           n+1       /              \      |
             |    sigma*u      ,  | rho*Du o nabla | v   |
             |           (i)       \              /      |
              \                                         /
  */
  /* supg stabilisation: pressure part, linearisation of test function if Newton ( L_pres_p) */
  /*
             /                                     \
            |         n+1    /                \     |
            |  nabla p    , |   rho*Du o nabla | v  |
            |         (i)    \                /     |
             \                                     /
  */
  /* supg stabilisation: viscous part, linearisation of test function if Newton (-L_visc_u) */
  /*
             /                                               \
            |               / n+1 \    /               \      |
            |  nabla o eps | u     |, |  rho*Du o nabla | v   |
            |               \ (i) /    \               /      |
             \                                               /
  */
  /* supg stabilisation: bodyforce part, linearisation of test function if Newton */
  /*
             /                                      \
            |                  /               \     |
            |  rho*rhsint   , |  rho*Du o nabla | v  |
            |                  \               /     |
             \                                      /
  */
  if (fldpara_->IsNewton())
  {
    if (fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
    {
      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        temp(jdim) = timefacfac * supgfac * momres_old_(jdim);
      }
    }
    else
    {
      for (int jdim = 0; jdim < nsd_; ++jdim)
      {
        temp(jdim) = -timefacfac * densaf_ * sgvelint_(jdim);
      }
    }
    for (int vi = 0; vi < nen_; ++vi)
    {
      for (int ui = 0; ui < nen_; ++ui)
      {
        for (int idim = 0; idim < nsd_; ++idim)
        {
          const double w = temp(idim) * funct_(ui);
          for (int jdim = 0; jdim < nsd_; ++jdim)
          {
            estif_u(nsd_ * vi + idim, nsd_ * ui + jdim) +=
                lin_resM_Du(nsd_ * idim + jdim, ui) * supg_test(vi) + derxy_(jdim, vi) * w;
          }  // jdim
        }    // vi
      }      // ui
    }        // idim
  }          // end if (fldpara_->IsNewton())
  else if (is_higher_order_ele_)
  {
    for (int vi = 0; vi < nen_; ++vi)
    {
      for (int ui = 0; ui < nen_; ++ui)
      {
        for (int idim = 0; idim < nsd_; ++idim)
        {
          for (int jdim = 0; jdim < nsd_; ++jdim)
          {
            estif_u(nsd_ * vi + idim, nsd_ * ui + jdim) +=
                lin_resM_Du(nsd_ * idim + jdim, ui) * supg_test(vi);
          }
        }
      }
    }
  }  // end if (is_higher_order_ele_)
  else
  {
    for (int vi = 0; vi < nen_; ++vi)
    {
      for (int ui = 0; ui < nen_; ++ui)
      {
        for (int idim = 0; idim < nsd_; ++idim)
        {
          estif_u(nsd_ * vi + idim, nsd_ * ui + idim) +=
              lin_resM_Du(nsd_ * idim + idim, ui) * supg_test(vi);
        }  // ui
      }    // idim
    }      // vi
  }        // end if not (is_higher_order_ele_) nor (newton_)

  /* supg stabilisation: pressure part  ( L_pres_p) */
  /*
           /                                    \
          |              /       n+1       \     |
          |  nabla Dp , |   rho*u   o nabla | v  |
          |              \       (i)       /     |
           \                                    /
  */
  for (int vi = 0; vi < nen_; ++vi)
  {
    //       const double v = timefacfacpre*supg_test(vi);
    double v;
    if (not fldparatimint_->IsFullImplPressureAndCont())
    {
      v = timefacfacpre * supg_test(vi);
    }
    else
    {
      v = fldparatimint_->Dt() * fac_ * supg_test(vi);
    }

    for (int ui = 0; ui < nen_; ++ui)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        estif_p_v(nsd_ * vi + idim, ui) += v * derxy_(idim, ui);
      }
    }
  }  // end for(idim)

  if (fldpara_->Tds() == INPAR::FLUID::subscales_quasistatic)
  {
    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      temp(jdim) = rhsfac * momres_old_(jdim);
    }
  }
  else
  {
    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      temp(jdim) = -rhsfac * densaf_ * sgvelint_(jdim) / (fac3 * fldparatimint_->AlphaF());
    }
  }

  for (int idim = 0; idim < nsd_; ++idim)
  {
    for (int vi = 0; vi < nen_; ++vi)
    {
      // supg stabilisation
      velforce(idim, vi) -= temp(idim) * supg_test(vi);
    }
  }  // end for(idim)

  // SUPG and Reynolds-stress term on right-hand side of
  // continuity equation for low-Mach-number flow
  if (fldpara_->PhysicalType() == INPAR::FLUID::loma and
      (fldpara_->ContiSUPG() or
          fldpara_->ContiReynolds() != INPAR::FLUID::reynolds_stress_stab_none))
  {
    const double temp_supg = rhsfac * scaconvfacaf_ * sgscaint_;

    if (fldpara_->ContiSUPG())
    {
      for (int vi = 0; vi < nen_; ++vi)
      {
        preforce(vi) -= temp_supg * conv_c_(vi);
      }
    }

    if (fldpara_->ContiReynolds() != INPAR::FLUID::reynolds_stress_stab_none)
    {
      for (int vi = 0; vi < nen_; ++vi)
      {
        preforce(vi) -= temp_supg * sgconv_c_(vi);
      }
    }
  }  // loma

  return;
}

/*---------------------------------------------------------------------------*
 | get ALE grid displacements and grid velocity for element     schott 11/14 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::GetGridDispVelALEOSTNew(
    DRT::Discretization& discretization, const std::vector<int>& lm,
    LINALG::Matrix<nsd_, nen_>& edispnp, LINALG::Matrix<nsd_, nen_>& egridvnp,
    LINALG::Matrix<nsd_, nen_>& egridvn)
{
  switch (fldpara_->PhysicalType())
  {
    case INPAR::FLUID::oseen:
    case INPAR::FLUID::stokes:
    {
      dserror(
          "ALE with Oseen or Stokes seems to be a tricky combination. Think deep before removing "
          "dserror!");
      break;
    }
    default:
    {
      GetGridDispALE(discretization, lm, edispnp);
      ExtractValuesFromGlobalVector(discretization, lm, *rotsymmpbc_, &egridvnp, NULL, "gridv");
      ExtractValuesFromGlobalVector(discretization, lm, *rotsymmpbc_, &egridvn, NULL, "gridvn");
      break;
    }
  }
}

/*---------------------------------------------------------------------------*
 |  set the (relative) convective velocity at integration point schott 11/14 |
 *---------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype, enrtype>::SetConvectiveVelintN(const bool isale)
{
  // get convective velocity at integration point
  switch (fldpara_->PhysicalType())
  {
    case INPAR::FLUID::incompressible:
    case INPAR::FLUID::artcomp:
    case INPAR::FLUID::varying_density:
    case INPAR::FLUID::loma:
    case INPAR::FLUID::tempdepwater:
    case INPAR::FLUID::boussinesq:
    case INPAR::FLUID::topopt:
    {
      convvelintn_.Update(velintn_);
      break;
    }
    case INPAR::FLUID::oseen:
    {
      dserror("not supported for new ost up to now");
      break;
    }
    case INPAR::FLUID::stokes:
    {
      convvelintn_.Clear();
      break;
    }
    default:
      dserror(
          "Physical type not implemented here. For Poro-problems see derived class "
          "FluidEleCalcPoro.");
      break;
  }

  // (ALE case handled implicitly here using the (potential
  //  mesh-movement-dependent) convective velocity, avoiding
  //  various ALE terms used to be calculated before)
  if (isale)
  {
    convvelintn_.Update(-1.0, gridvelintn_, 1.0);
  }
}

// template classes
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex8, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex8, DRT::ELEMENTS::Fluid::xwall>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet4, DRT::ELEMENTS::Fluid::xwall>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex20, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex27, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet4, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet10, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::wedge6, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::wedge15, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::pyramid5, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad4, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad8, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad9, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tri3, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tri6, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::nurbs9, DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::nurbs27, DRT::ELEMENTS::Fluid::none>;
