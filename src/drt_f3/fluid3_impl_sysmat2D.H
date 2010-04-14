/*----------------------------------------------------------------------*/
/*!
\file fluid2_impl.cpp

\brief Internal implementation of Fluid2 element

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#ifdef D_FLUID3
#ifdef CCADISCRET

#include "fluid3_impl.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/mixfrac.H"
#include "../drt_mat/sutherland.H"
#include "../drt_mat/arrhenius_pv.H"
#include "../drt_mat/ferech_pv.H"
#include "../drt_mat/carreauyasuda.H"
#include "../drt_mat/modpowerlaw.H"
#include "../drt_lib/drt_timecurve.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_lib/drt_condition_utils.H"

#include <Epetra_SerialDenseSolver.h>

#ifdef DEBUG
#endif

/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)   g.bau 03/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::Sysmat2D(
  Fluid3*                                 ele,
  const LINALG::Matrix<nsd_,nen_>&            evelaf,
  const LINALG::Matrix<nsd_,nen_>&            eveln,
  const LINALG::Matrix<nsd_,nen_>&            fsevelaf,
  const LINALG::Matrix<nen_,1>&            epreaf,
  const LINALG::Matrix<nsd_,nen_>&            eaccam,
  const LINALG::Matrix<nen_,1>&            escaaf,
  const LINALG::Matrix<nen_,1>&            escaam,
  const LINALG::Matrix<nen_,1>&            escadtam,
  const LINALG::Matrix<nsd_,nen_>&            emhist,
  const LINALG::Matrix<nsd_,nen_>&            edispnp,
  const LINALG::Matrix<nsd_,nen_>&            egridv,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&            estif,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&            emesh,
  LINALG::Matrix<(nsd_+1)*nen_,    1>&            eforce,
  const double                                   thermpressaf,
  const double                                   thermpressam,
  const double                                   thermpressdtam,
  Teuchos::RCP<const MAT::Material>       material,
  const double                                  timefac,
  double&                                 Cs
  )
{
  // get node coordinates and number of elements per node
  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<nen_; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
  }

  // add displacement when fluid nodes move in the ALE case
  if (ele->is_ale_) xyze_ += edispnp;

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  BodyForce(ele);

  // evaluate shape functions and derivatives at element center
  EvalShapeFuncAndDerivsAtEleCenter(ele->Id());

  // element volume
  const double area = fac_;

  // in case of viscous stabilization decide whether to use GLS or USFEM
  double vstabfac= 0.0;
  if (vstab_ == Fluid3::viscous_stab_usfem or
      vstab_ == Fluid3::viscous_stab_usfem_only_rhs)   vstabfac =  1.0;
  else if(vstab_ == Fluid3::viscous_stab_gls or
          vstab_ == Fluid3::viscous_stab_gls_only_rhs) vstabfac = -1.0;

  //----------------------------------------------------------------------
  // get material parameters at element center
  //----------------------------------------------------------------------
  if (not mat_gp_ or not tau_gp_)
    GetMaterialParams(material,evelaf,escaaf,escaam,thermpressaf,thermpressam,thermpressdtam);

  if (not tau_gp_)
  {
    // ---------------------------------------------------------------------
    // calculate all-scale or fine-scale subgrid viscosity at element center
    // ---------------------------------------------------------------------
    visceff_ = visc_;
    if (turb_mod_action_ == Fluid3::smagorinsky)
    {
      CalcSubgrVisc_2D(evelaf,area,Cs);

      // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
      visceff_ += sgvisc_;
    }
    else if (fssgv_ != Fluid3::no_fssgv)
      CalcFineScaleSubgrVisc_2D(evelaf,fsevelaf,area,Cs);

    // get velocity at element center
    velint_.Multiply(evelaf,funct_);

    // ---------------------------------------------------------------------
    // calculate stabilization parameter at element center
    // ---------------------------------------------------------------------
    CalcStabParameter_2D(timefac,area);
  }

  // Gaussian integration points
  // const DRT::UTILS::IntegrationPoints2D intpoints(ele->gaussrule_);
  const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // integration loop
  for (int iquad=0; iquad<intpoints.IP().nquad; ++iquad)
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());

    //----------------------------------------------------------------------
    // get material parameters (evaluation at integration point)
    //----------------------------------------------------------------------
    if (mat_gp_) GetMaterialParams(material,evelaf,escaaf,escaam,thermpressaf,thermpressam,thermpressdtam);

    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    velint_.Multiply(evelaf,funct_);

    if (tau_gp_)
    {
      // ---------------------------------------------------------------------
      // calculate all-scale or fine-scale subgrid viscosity at element center
      // ---------------------------------------------------------------------
      visceff_ = visc_;
      if (turb_mod_action_ == Fluid3::smagorinsky)
      {
        CalcSubgrVisc_2D(evelaf,area,Cs);

        // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
        visceff_ += sgvisc_;
      }
      else if (fssgv_ != Fluid3::no_fssgv)
        CalcFineScaleSubgrVisc_2D(evelaf,fsevelaf,area,Cs);

      // ---------------------------------------------------------------------
      // calculate stabilization parameter at element center
      // ---------------------------------------------------------------------
      CalcStabParameter_2D(timefac,area);
    }

    // get momentum history data at integration point
    histmom_.Multiply(emhist,funct_);

    // get velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    vderxy_.MultiplyNT(evelaf,derxy_);

    // get fine-scale velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    if (fssgv_ != Fluid3::no_fssgv) fsvderxy_.MultiplyNT(fsevelaf,derxy_);
    else                           fsvderxy_.Clear();

    // get convective velocity at integration point
    // We handle the ale case very implicitly here using the (possible mesh
    // movement dependent) convective velocity. This avoids a lot of ale terms
    // we used to calculate.
    convvelint_.Update(velint_);
    if (ele->is_ale_) convvelint_.Multiply(-1.0, egridv, funct_, 1.0);

    // get pressure gradient at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    gradp_.Multiply(derxy_,epreaf);

    // get pressure at integration point
    // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double press = funct_.Dot(epreaf);

    // get bodyforce at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    bodyforce_.Multiply(edeadaf_,funct_);

    //--------------------------------------------------------------------
    // get numerical representation of some single operators
    //--------------------------------------------------------------------
    if (is_higher_order_ele_)
    {
      /*--- viscous term: div(epsilon(u)) -------------------------------*/
      /*     /                              \
           1 |  2 N_x,xx + N_x,yy + N_y,xy  |    with N_x .. x-line of N
           - |                              |         N_y .. y-line of N
           2 |  N_y,xx + N_x,yx + 2 N_y,yy  |
             \                              /                            */

      /*--- subtraction for low-Mach-number flow: div((1/3)*(div u)*I) */
      /*    /                   \
          1 |  N_x,xx + N_y,yx  |
       -  - |                   |
          3 |  N_x,xy + N_y,yy  |
            \                   /

               with N_x .. x-line of N
               N_y .. y-line of N                                      */

      // varying density: div u = 0
      // therefore, no stabilization of the second term
      double prefac;
      if (physicaltype_ == INPAR::FLUID::loma)
      {
        prefac = 1.0/3.0;
        derxy2_.Scale(prefac);
      }
      else prefac = 1.0;

      double sum = (derxy2_(0,0)+derxy2_(1,0))/prefac;

      viscs2_(0,0) = 0.5 * (sum + derxy2_(0,0));
      viscs2_(1,0) = 0.5 * derxy2_(2,0);
      viscs2_(3,0) = 0.5 * (sum + derxy2_(1,0));

      /* viscous term  div epsilon(u_old) */
      visc_old_(0) = viscs2_(0,0)*evelaf(0,0)+viscs2_(1,0)*evelaf(1,0);
      visc_old_(1) = viscs2_(1,0)*evelaf(0,0)+viscs2_(3,0)*evelaf(1,0);

      for (int i=1; i<nen_; ++i)
      {
        double sum = (derxy2_(0,i)+derxy2_(1,i))/prefac;

        viscs2_(0,i) = 0.5 * (sum + derxy2_(0,i));
        viscs2_(1,i) = 0.5 * derxy2_(2,i);
        viscs2_(3,i) = 0.5 * (sum + derxy2_(1,i));

        /* viscous term  div epsilon(u_old) */
        visc_old_(0) += viscs2_(0,i)*evelaf(0,i)+viscs2_(1,i)*evelaf(1,i);
        visc_old_(1) += viscs2_(1,i)*evelaf(0,i)+viscs2_(3,i)*evelaf(1,i);
      }
    }
    else
    {
      viscs2_.Clear();
      visc_old_.Clear();
    }

    // convective term from previous iteration
    conv_old_.Multiply(vderxy_,convvelint_);

    // compute convective operator
    conv_c_.MultiplyTN(derxy_,convvelint_);

    // velocity divergence from previous iteration
    vdiv_ = vderxy_(0, 0) + vderxy_(1, 1);

    //--------------------------------------------------------------------
    // factors for stabilization, time integration
    // and fine-scale subgrid-viscosity
    //--------------------------------------------------------------------
    const double tau_M       = tau_(0)*fac_;
    const double tau_Mp      = tau_(1)*fac_;
    const double tau_C       = tau_(2)*fac_;

    const double timefacfac  = timefac * fac_;
    const double timetauM    = timefac * tau_M;
    const double timetauMp   = timefac * tau_Mp;

    double rhsfac            = fac_;

    const double fssgviscfac = fssgvisc_*timefacfac;

    //--------------------------------------------------------------------
    // The following computations are performed depending on
    // time-integration, that is, whether it is generalized-alpha or not,
    // since several terms differ with respect to the scheme:
    //
    // 1) calculation of rhs for momentum equation and momentum residual
    // -> different for generalized-alpha and other schemes
    //
    // 2) calculation of additional subgrid-scale velocity when cross-
    //    and Reynolds-stress are included:
    // - Cross- and Reynolds-stress are always included simultaneously.
    // - They are included in a complete form on left- and right-hand side.
    // - For this purpose, a subgrid-scale convective term is computed.
    // - Within a Newton linearization, the present formulation is not
    //   consistent for the reactive terms.
    // - To turn them off, both flags must be "no".
    //
    // 3) calculation of convective scalar term, rhs for continuity
    //    equation and residual of continuity equations
    // -> only required for low-Mach-number flow
    // -> for incompressible flow, residual is velocity divergence only
    //--------------------------------------------------------------------
    if (is_genalpha_)
    {
      // rhs of momentum equation: density*bodyforce at n+alpha_F
      rhsmom_.Update(densaf_,bodyforce_,0.0);

      // get acceleration at time n+alpha_M at integration point
      accint_.Multiply(eaccam,funct_);

      // evaluate momentum residual once for all stabilization right hand sides
      for (int rr=0;rr<2;++rr)
      {
        momres_old_(rr) = densam_*accint_(rr)+densaf_*conv_old_(rr)+gradp_(rr)-2*visceff_*visc_old_(rr)-densaf_*bodyforce_(rr);
      }

      if (cross_    != Fluid3::cross_stress_stab_none or
          reynolds_ != Fluid3::reynolds_stress_stab_none)
      {
        // compute subgrid-scale velocity
        sgvelint_.Update(-tau_(1),momres_old_,0.0);

        // compute subgrid-scale convective operator
        sgconv_c_.MultiplyTN(derxy_,sgvelint_);

        // re-calculate convective term from previous iteration if cross-stress
        // is included
        convvelint_.Update(1.0,sgvelint_,1.0);
        conv_old_.Multiply(vderxy_,convvelint_);
      }
      else sgconv_c_.Clear();

      // "incompressible" part of continuity residual: velocity divergence
      conres_old_ = vdiv_;

      if (physicaltype_ == INPAR::FLUID::loma)
      {
        // time derivative of scalar at n+alpha_M
        const double tder_sca = funct_.Dot(escadtam);

        // gradient of scalar value at n+alpha_F
        grad_scaaf_.Multiply(derxy_,escaaf);

        // convective scalar term at n+alpha_F
        conv_scaaf_ = velint_.Dot(grad_scaaf_);

        // add subgrid-scale velocity part also to convective scalar term
        // -> currently not taken into account
        /*if (cross    != Fluid3::cross_stress_stab_none or
            reynolds != Fluid3::reynolds_stress_stab_none)
          conv_scaaf_ += sgvelint_.Dot(grad_scaaf_);*/

        // rhs of continuity equation (only relevant for low-Mach-number flow)
        rhscon_ = scadtfac_*tder_sca + scaconvfacaf_*conv_scaaf_ + thermpressadd_;

        // residual of continuity equation
        conres_old_ -= rhscon_;
      }
    }
    else
    {
      // rhs of momentum equation:
      // density*timefac*bodyforce at n+1 + density*histmom at n

      // in the case of a Boussinesq approximation: f = (rho - rho_0)/rho_0 *g
      // else:                    f = rho * g
      if(physicaltype_ == INPAR::FLUID::boussinesq)
      {
        rhsmom_.Update(densn_,histmom_,deltadens_*timefac,bodyforce_);
      }
      else rhsmom_.Update(densn_,histmom_,densaf_*timefac,bodyforce_);

      // modify integration factor for Galerkin rhs
      rhsfac *= timefac;

      // evaluate momentum residual once for all stabilization right hand sides
      for (int rr=0;rr<2;++rr)
      {
        momres_old_(rr) = densaf_*velint_(rr)+timefac*(densaf_*conv_old_(rr)+gradp_(rr)-2*visceff_*visc_old_(rr))-rhsmom_(rr);
      }

      if (cross_    != Fluid3::cross_stress_stab_none or
          reynolds_ != Fluid3::reynolds_stress_stab_none)
      {
        // compute subgrid-scale velocity
        sgvelint_.Update(-(tau_(1)/dt_),momres_old_,0.0);

        // compute subgrid-scale convective operator
        sgconv_c_.MultiplyTN(derxy_,sgvelint_);

        // re-calculate convective term from previous iteration if cross-stress
        // is included
        convvelint_.Update(1.0,sgvelint_,1.0);
        conv_old_.Multiply(vderxy_,convvelint_);
      }
      else sgconv_c_.Clear();

      // "incompressible" part of continuity residual: velocity divergence
      conres_old_ = timefac*vdiv_;

      if (physicaltype_ == INPAR::FLUID::loma or physicaltype_ == INPAR::FLUID::varying_density)
      {
        // get velocity derivatives at n
        vderxyn_.MultiplyNT(eveln,derxy_);

        // velocity divergence at n
        const double vdivn = vderxyn_(0, 0) + vderxyn_(1, 1);

        // time derivative of scalar
        const double tder_sca = funct_.Dot(escadtam);

        // gradient of scalar value at n+1
        grad_scaaf_.Multiply(derxy_,escaaf);

        // convective scalar term at n+1
        conv_scaaf_ = velint_.Dot(grad_scaaf_);

        // gradient of scalar value at n
        grad_scan_.Multiply(derxy_,escaam);

        // convective scalar term at n
        conv_scan_ = velint_.Dot(grad_scan_);

        // add subgrid-scale velocity part also to convective scalar term
        // (subgrid-scale velocity at n+1 also approximately used at n)
        // -> currently not taken into account
        /*if (cross    != Fluid3::cross_stress_stab_none or
            reynolds != Fluid3::reynolds_stress_stab_none)
        {
          conv_scaaf_ += sgvelint_.Dot(grad_scaaf_);
          conv_scan_  += sgvelint_.Dot(grad_scan_);
        }*/

        // rhs of continuity equation (only relevant for low-Mach-number flow)
        rhscon_ = scadtfac_*tder_sca + timefac*scaconvfacaf_*conv_scaaf_ + omtheta_*dt_*(scaconvfacn_*conv_scan_-vdivn) + thermpressadd_;

        // residual of continuity equation
        conres_old_ -= rhscon_;
      }
    }

    //------------------------------------------------------------------------
    // perform integration for element matrix and right hand side
    //------------------------------------------------------------------------
    {
      //----------------------------------------------------------------------
      //                            GALERKIN PART

      //----------------------------------------------------------------------
      // computation of inertia term and convection term (convective and
      // reactive part) for convective form of convection term including
      // right-hand-side contribution and potential cross-stress term
      //----------------------------------------------------------------------
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui   = 3*ui;
        const int fuip  = fui+1;
        const double v = fac_*densam_*funct_(ui)
#if 1
                         + timefacfac*densaf_*(conv_c_(ui)+sgconv_c_(ui))
#endif
                         ;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 3*vi;
          const int fvip  = fvi+1;
          /* inertia (contribution to mass matrix) */
          /*
          /                \
          |                |
          |    rho*Du , v  |
          |                |
          \                /
          */
          /* convection, convective part (convective form) */
          /*
          /                               \
          |  /       n+1        \         |
          | |   rho*u   o nabla | Du , v  |
          |  \      (i)        /          |
          \                              /
          */
          double v2 = v*funct_(vi) ;
          estif(fvi  , fui  ) += v2;
          estif(fvip , fuip ) += v2;
        }
      }

      if (is_newton_)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 3*vi;
          const int fvip  = fvi+1;
          const double v = timefacfac*densaf_*funct_(vi);
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = 3*ui;
            const int fuip  = fui+1;
            const double v2 = v*funct_(ui);
            /*  convection, reactive part (convective form)
            /                                 \
            |  /                \   n+1       |
            | |  rho*Du o nabla | u     , v   |
            |  \                /   (i)       |
            \                                /
            */
            estif(fvi,   fui)   += v2*vderxy_(0, 0) ;
            estif(fvi,   fuip)  += v2*vderxy_(0, 1) ;
            estif(fvip,  fui)   += v2*vderxy_(1, 0) ;
            estif(fvip,  fuip)  += v2*vderxy_(1, 1) ;
          }
        }
      }

      if (is_genalpha_)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = 3*vi;
          /* inertia term on right-hand side for generalized-alpha scheme */
          const double v = -fac_*densam_*funct_(vi);
          eforce(fvi    ) += v*accint_(0) ;
          eforce(fvi + 1) += v*accint_(1) ;
        }
      }
      else
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = 3*vi;
          /* inertia term on right-hand side for one-step-theta/BDF2 schem */
          const double v = -fac_*densaf_*funct_(vi);
          eforce(fvi    ) += v*velint_(0) ;
          eforce(fvi + 1) += v*velint_(1) ;
        }
      }

#if 1
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi   = 3*vi;
        /* convection (convective form) on right-hand side */
        double v = -rhsfac*densaf_*funct_(vi);
        eforce(fvi    ) += v*conv_old_(0) ;
        eforce(fvi + 1) += v*conv_old_(1) ;
      }
#endif

      //----------------------------------------------------------------------
      // computation of additions to convection term (convective and
      // reactive part) for conservative form of convection term including
      // right-hand-side contribution
      //----------------------------------------------------------------------
      if (is_conservative_)
      {
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 3*ui;
          const int fuip  = fui+1;
          double v = timefacfac*densaf_*funct_(ui)*vdiv_;
          if (physicaltype_ == INPAR::FLUID::loma) v -= timefacfac*densaf_*scaconvfacaf_*conv_scaaf_;
          // only with linear density-concentration correlation
          else if(physicaltype_ == INPAR::FLUID::varying_density) v += timefacfac*conv_scaaf_;
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 3*vi;
            const int fvip  = fvi+1;
            /* convection, convective part (conservative addition) */
            /*
            /                                                   \
            |      /              n+1    n+1           \      |
            |  Du | rho*nabla o u    +  u   *nabla rho | , v  |
            |      \             (i)     (i)          /       |
            \                                                 /
            */
            double v2 = v*funct_(vi) ;
            estif(fvi  , fui  ) += v2;
            estif(fvip , fuip ) += v2;
          }
        }

        if (is_newton_)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 3*vi;
            const int fvip  = fvi+1;
            const double v0 = timefacfac*densaf_*velint_(0)*funct_(vi);
            const double v1 = timefacfac*densaf_*velint_(1)*funct_(vi);
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 3*ui;
              const int fuip  = fui+1;
              /*  convection, reactive part (conservative addition) */
              /*
              /                              \
              |  n+1  /               \      |
              | u    | rho*nabla o Du | , v  |
              |  (i)  \              /       |
              \                             /
              */
              estif(fvi,  fui  ) += v0*derxy_(0, ui) ;
              estif(fvi,  fuip ) += v0*derxy_(1, ui) ;
              estif(fvip, fui  ) += v1*derxy_(0, ui) ;
              estif(fvip, fuip ) += v1*derxy_(1, ui) ;
            }
          }

          if (physicaltype_ == INPAR::FLUID::loma)
          {
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = 3*vi;
              const int fvip  = fvi+1;
              const double v0 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(0)*velint_(0)*funct_(vi);
              const double v1 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(1)*velint_(1)*funct_(vi);
              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui   = 3*ui;
                const int fuip  = fui+1;
                /*  convection, reactive part (conservative addition) */
                /*
                /                           \
                |  n+1  /             \      |
                | u    | Du*nabla rho | , v  |
                |  (i)  \            /       |
                \                           /
                */
                estif(fvi,  fui  ) += v0*funct_(ui) ;
                estif(fvi,  fuip ) += v0*funct_(ui) ;
                estif(fvip, fui  ) += v1*funct_(ui) ;
                estif(fvip, fuip ) += v1*funct_(ui) ;
              }
            }
          }
          if (physicaltype_ == INPAR::FLUID::varying_density)
          {
          for (int vi=0; vi<nen_; ++vi)
            {
        const int fvi   = 3*vi;
        const int fvip  = fvi+1;
        const double v0 = +timefacfac*grad_scaaf_(0)*velint_(0)*funct_(vi);
        const double v1 = +timefacfac*grad_scaaf_(1)*velint_(1)*funct_(vi);
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 3*ui;
          const int fuip  = fui+1;
          /*  convection, reactive part (conservative addition) */
          /*
          /                           \
          |  n+1  /             \      |
          | u    | Du*nabla rho | , v  |
          |  (i)  \            /       |
          \                           /
          */
          estif(fvi,  fui  ) += v0*funct_(ui) ;
          estif(fvi,  fuip ) += v0*funct_(ui) ;
          estif(fvip, fui  ) += v1*funct_(ui) ;
          estif(fvip, fuip ) += v1*funct_(ui) ;
        }
            }
          }
        }

        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 3*vi;
          /* convection (conservative addition) on right-hand side */
          double v = -rhsfac*densaf_*funct_(vi)*vdiv_;
          eforce(fvi    ) += v*velint_(0) ;
          eforce(fvi + 1) += v*velint_(1) ;
        }

        if (physicaltype_ == INPAR::FLUID::loma)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 3*vi;
            /* convection (conservative addition) on rhs for low-Mach-number flow */
            double v = rhsfac*densaf_*scaconvfacaf_*conv_scaaf_*funct_(vi);
            eforce(fvi    ) += v*velint_(0) ;
            eforce(fvi + 1) += v*velint_(1) ;
          }
        }
        if (physicaltype_ == INPAR::FLUID::varying_density)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 3*vi;
            /* convection (conservative addition) on rhs for varying_density flow */
            double v = -rhsfac*conv_scaaf_*funct_(vi);
            eforce(fvi    ) += v*velint_(0) ;
            eforce(fvi + 1) += v*velint_(1) ;
          }
        }

      }

      //----------------------------------------------------------------------
      // computation of viscosity term including right-hand-side contribution
      //----------------------------------------------------------------------
      const double visceff_timefacfac = visceff_*timefacfac;
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui   = 3*ui;
        const int fuip  = fui+1;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 3*vi;
          const int fvip  = fvi+1;

          const double derxy_0ui_0vi = derxy_(0, ui)*derxy_(0, vi);
          const double derxy_1ui_1vi = derxy_(1, ui)*derxy_(1, vi);
          /* viscosity term */
          /*
                /                          \
                |       /  \         / \   |
          2 mu  |  eps | Du | , eps | v |  |
                |       \  /         \ /   |
                \                          /
          */
          estif(fvi, fui)     += visceff_timefacfac*(2.0*derxy_0ui_0vi
                                                     +
                                                     derxy_1ui_1vi) ;
          estif(fvi , fuip)   += visceff_timefacfac*derxy_(0, ui)*derxy_(1, vi) ;
          estif(fvip, fui)    += visceff_timefacfac*derxy_(0, vi)*derxy_(1, ui) ;
          estif(fvip, fuip)   += visceff_timefacfac*(derxy_0ui_0vi
                                                     +
                                                     2.0*derxy_1ui_1vi) ;
        }
      }

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = 3*vi;
        const double v = -visceff_*rhsfac;
        /* viscosity term on right-hand side */
        eforce(fvi)     += v*(2.0*derxy_(0, vi)*vderxy_(0, 0)
                              +
                              derxy_(1, vi)*vderxy_(0, 1)
                              +
                              derxy_(1, vi)*vderxy_(1, 0)) ;
        eforce(fvi + 1) += v*(derxy_(0, vi)*vderxy_(0, 1)
                              +
                              derxy_(0, vi)*vderxy_(1, 0)
                              +
                              2.0*derxy_(1, vi)*vderxy_(1, 1)) ;
      }

      //----------------------------------------------------------------------
      // computation of pressure term including right-hand-side contribution
      //----------------------------------------------------------------------
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fuippp = 3*ui+2;
        const double v = -timefacfac*funct_(ui);
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 3*vi;
          /* pressure term */
          /*
          /                  \
          |                  |
          |  Dp , nabla o v  |
          |                  |
          \                  /
          */
          estif(fvi,     fuippp) += v*derxy_(0, vi) ;
          estif(fvi + 1, fuippp) += v*derxy_(1, vi) ;
        }
      }

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = 3*vi;
        /* pressure term on right-hand side */
        const double v = press*rhsfac;
        eforce(fvi    ) += v*derxy_(0, vi) ;
        eforce(fvi + 1) += v*derxy_(1, vi) ;
      }

      //----------------------------------------------------------------------
      // computation of continuity term including right-hand-side contribution
      //----------------------------------------------------------------------
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvippp = 3*vi+2;
        const double v = timefacfac*funct_(vi);
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 3*ui;
          /* continuity term */
          /*
            /                  \
            |                  |
            | nabla o Du  , q  |
            |                  |
            \                  /
          */
          estif(fvippp, fui)     += v*derxy_(0, ui) ;
          estif(fvippp, fui + 1) += v*derxy_(1, ui) ;
        }
      }

      const double rhsfac_vdiv = -rhsfac * vdiv_;
      for (int vi=0; vi<nen_; ++vi)
      {
        // continuity term on right-hand side
        eforce(vi*3 + 2) += rhsfac_vdiv*funct_(vi) ;
      }

      //----------------------------------------------------------------------
      // computation of body-force term on right-hand side
      //----------------------------------------------------------------------
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = 3*vi;
        const double v = fac_*funct_(vi);
        eforce(fvi    ) += v*rhsmom_(0) ;
        eforce(fvi + 1) += v*rhsmom_(1) ;
      }

      //----------------------------------------------------------------------
      // computation of additional terms for low-Mach-number flow:
      // 1) subtracted viscosity term including right-hand-side contribution
      // 2) additional rhs term of continuity equation
      //----------------------------------------------------------------------

      if (physicaltype_ == INPAR::FLUID::loma)
      {
        const double v = -(2.0/3.0)*visceff_*timefacfac ;
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 3*ui;
          const int fuip  = fui+1;
          const double v0 = v*derxy_(0,ui);
          const double v1 = v*derxy_(1,ui);
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 3*vi;
            const int fvip  = fvi+1;
            /* viscosity term - subtraction for low-Mach-number flow */
            /*
                  /                               \
                  |  1                      / \   |
           - 2 mu |  - (nabla o u) I , eps | v |  |
                  |  3                      \ /   |
                  \                               /
            */
            estif(fvi,   fui  ) += v0*derxy_(0, vi) ;
            estif(fvi,   fuip ) += v1*derxy_(0, vi) ;
            estif(fvip,  fui  ) += v0*derxy_(1, vi) ;
            estif(fvip,  fuip ) += v1*derxy_(1, vi) ;
          }
        }

        const double v_div = (2.0/3.0)*visceff_*rhsfac*vdiv_ ;
        const double fac_rhscon = fac_*rhscon_;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = 3*vi;
          /* viscosity term on rhs - subtraction for low-Mach-number flow */
          eforce(fvi    ) += derxy_(0, vi)*v_div ;
          eforce(fvi + 1) += derxy_(1, vi)*v_div ;

          /* additional rhs term of continuity equation */
          eforce(fvi + 2) += fac_rhscon*funct_(vi) ;
        }
      }

      //----------------------------------------------------------------------
      //                 PRESSURE STABILISATION PART

      if (pspg_ == Fluid3::pstab_use_pspg)
      {
        for (int ui=0; ui<nen_; ++ui)
        {
          const int tui  = 3*ui;
          const int tuip = tui+1;
          const double v = tau_Mp*densam_*funct_(ui)
#if 1
                           + timetauMp*densaf_*conv_c_(ui)
#endif
                           ;
          for (int vi=0; vi<nen_; ++vi)
          {
            const int tvipp  = 3*vi + 2;

            /* pressure stabilisation: inertia */
            /*
              /                    \
              |                    |
              |  rho*Du , nabla q  |
              |                    |
              \                   /
            */
            /* pressure stabilisation: convection, convective part */
            /*

            /                                     \
            |  /       n+1        \               |
            | |   rho*u   o nabla | Du , nabla q  |
            |  \      (i)         /               |
            \                                    /

            */

            estif(tvipp, tui ) += v*derxy_(0, vi) ;
            estif(tvipp, tuip) += v*derxy_(1, vi) ;
          }
        }

        if (is_higher_order_ele_)
        {
          const double v = -2.0*visceff_*timetauMp;
          for (int ui=0; ui<nen_; ++ui)
          {
            const int tui  = 3*ui;
            const int tuip = tui+1;
            for (int vi=0; vi<nen_; ++vi)
            {
              const int tvipp  = 3*vi + 2;

              /* pressure stabilisation: viscosity (-L_visc_u) */
              /*
                /                              \
                |               /  \             |
                |  nabla o eps | Du | , nabla q  |
                |               \  /             |
                \                              /
              */
              estif(tvipp, tui ) += v*(derxy_(0, vi)*viscs2_(0, ui)
                                       +
                                       derxy_(1, vi)*viscs2_(1, ui)) ;
              estif(tvipp, tuip) += v*(derxy_(0, vi)*viscs2_(1, ui)
                                       +
                                       derxy_(1, vi)*viscs2_(3, ui)) ;
            }
          }
        }

        for (int ui=0; ui<nen_; ++ui)
        {
          const int tuipp = 3*ui + 2;
          for (int vi=0; vi<nen_; ++vi)
          {
            /* pressure stabilisation: pressure( L_pres_p) */
            /*
              /                    \
              |                      |
              |  nabla Dp , nabla q  |
              |                      |
              \                    /
            */
            estif(vi*3 + 2, tuipp) += timetauMp*(derxy_(0, ui)*derxy_(0, vi)
                                                 +
                                                 derxy_(1, ui)*derxy_(1, vi)) ;

          } // vi
        } // ui

        if (is_newton_)
        {
          for (int ui=0; ui<nen_; ++ui)
          {
            const int tui  = 3*ui;
            const int tuip = tui+1;
            const double v = timetauMp*densaf_*funct_(ui);
            for (int vi=0; vi<nen_; ++vi)
            {
              const int tvipp = 3*vi + 2;
              /*  pressure stabilisation: convection, reactive part

              /                                     \
              |  /                 \  n+1           |
              | |   rho*Du o nabla | u     , grad q |
              |  \                /   (i)           |
              \                                     /

              */
              estif(tvipp, tui ) += v*(derxy_(0, vi)*vderxy_(0, 0)
                                       +
                                       derxy_(1, vi)*vderxy_(1, 0)) ;
              estif(tvipp, tuip) += v*(derxy_(0, vi)*vderxy_(0, 1)
                                       +
                                       derxy_(1, vi)*vderxy_(1, 1)) ;

            } // vi
          } // ui
        } // if newton

        for (int vi=0; vi<nen_; ++vi)
        {
          // pressure stabilisation
          eforce(vi*3 + 2) -= tau_Mp*(momres_old_(0)*derxy_(0, vi)
                                      +
                                      momres_old_(1)*derxy_(1, vi)) ;
        }
      }

      //----------------------------------------------------------------------
      //                     SUPG STABILISATION PART

      if(supg_ == Fluid3::convective_stab_supg)
      {
#if 1
        for (int ui=0; ui<nen_; ++ui)
        {
          const int tui  = 3*ui;
          const int tuip = tui+1;
          const double v = densaf_*(tau_M*densam_*funct_(ui) + timetauM*densaf_*conv_c_(ui));
          for (int vi=0; vi<nen_; ++vi)
          {
            const int tvi = 3*vi;
            /* supg stabilisation: inertia  */
            /*
              /                                 \
              |           /      n+1        \    |
              |  rho*Du , | rho*u   o nabla | v  |
              |           \      (i)        /    |
              \                                 /
            */

            /* supg stabilisation: convective part ( L_conv_u) */

            /*

            /                                                      \
            |    /       n+1         \      /      n+1        \     |
            |   |   rho*u    o nabla | Du , | rho*u    o nabla | v  |
            |    \       (i)         /      \      (i)        /     |
            \                                                       /

            */

            const double v2 = v*(conv_c_(vi)+sgconv_c_(vi));
            estif(tvi,     tui ) += v2;
            estif(tvi + 1, tuip) += v2;
          }
        }

        for (int vi=0; vi<nen_; ++vi)
        {
          const int tvi  = 3*vi;
          const int tvip = tvi+1;
          const double v = timetauM*densaf_*(conv_c_(vi)+sgconv_c_(vi));
          for (int ui=0; ui<nen_; ++ui)
          {
            const int tuipp = 3*ui + 2;
            /* supg stabilisation: pressure part  ( L_pres_p) */
            /*
              /                                      \
              |              /       n+1       \     |
              |  nabla Dp , |   rho*u   o nabla | v  |
              |              \       (i)       /     |
              \                                     /
            */
            estif(tvi,  tuipp) += v*derxy_(0, ui) ;
            estif(tvip, tuipp) += v*derxy_(1, ui) ;
          }
        }

        if (is_higher_order_ele_)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int tvi  = 3*vi;
            const int tvip = tvi+1;
            const double v = -2.0*visceff_*timetauM*densaf_*(conv_c_(vi)+sgconv_c_(vi));
            for (int ui=0; ui<nen_; ++ui)
            {
              const int tui  = 3*ui;
              const int tuip = tui+1;
              /* supg stabilisation: viscous part  (-L_visc_u) */
              /*
                /                                                \
                |               /  \    /       n+1        \     |
                |  nabla o eps | Du |, |   rho*u    o nabla | v  |
                |               \  /    \       (i)        /     |
                \                                                /
              */
              estif(tvi,  tui ) += v*viscs2_(0, ui) ;
              estif(tvip, tui ) += v*viscs2_(1, ui) ;

              estif(tvi,  tuip) += v*viscs2_(1, ui) ;
              estif(tvip, tuip) += v*viscs2_(3, ui) ;
            }
          }
        }
#endif

        if (is_newton_)
        {
          for (int ui=0; ui<nen_; ++ui)
          {
            const int tui  = 3*ui;
            const int tuip = tui+1;
            const double v = tau_M*densam_*densaf_*funct_(ui);
            const double v0 = v*velint_(0);
            const double v1 = v*velint_(1);
            for (int vi=0; vi<nen_; ++vi)
            {
              const int tvi  = 3*vi;
              const int tvip = tvi+1;
              /* supg stabilisation: inertia, linearisation of testfunction  */
              /*
                /                                         \
                |         n+1      /                \     |
                |    rho*u      ,  | rho*Du o nabla | v   |
                |         (i)      \                /     |
                \                                         /

              */
              estif(tvi,  tui ) += v0*derxy_(0, vi) ;
              estif(tvip, tui ) += v1*derxy_(0, vi) ;

              estif(tvi,  tuip) += v0*derxy_(1, vi) ;
              estif(tvip, tuip) += v1*derxy_(1, vi) ;
            }
          }

#if 1
          {
            const double v0 = convvelint_(0)*vderxy_(0, 0) + convvelint_(1)*vderxy_(0, 1);
            const double v1 = convvelint_(0)*vderxy_(1, 0) + convvelint_(1)*vderxy_(1, 1);

            for (int ui=0; ui<nen_; ++ui)
            {
              const int tui  = 3*ui;
              const int tuip = tui+1;
              const double v = timetauM*densaf_*densaf_*funct_(ui);
              for (int vi=0; vi<nen_; ++vi)
              {
                const int tvi  = 3*vi;
                const int tvip = tvi+1;

                /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
                /*
                  /                                                         \
                  |    /       n+1        \   n+1    /                \     |
                  |   |   rho*u    o nabla | u    ,  | rho*Du o nabla | v   |
                  |    \       (i)        /   (i)    \                /     |
                  \                                                        /

                  /                                                         \
                  |    /                \   n+1    /       n+1        \     |
                  |    | rho*Du o nabla | u    ,   | rho*u   o nabla  | v   |
                  |    \                /   (i)    \       (i)        /     |
                  \                                                        /
                */
                estif(tvi,  tui ) += (conv_c_(vi)*vderxy_(0, 0) + v0*derxy_(0, vi))*v;
                estif(tvip, tui ) += (conv_c_(vi)*vderxy_(1, 0) + v1*derxy_(0, vi))*v;

                estif(tvi,  tuip) += (conv_c_(vi)*vderxy_(0, 1) + v0*derxy_(1, vi))*v;
                estif(tvip, tuip) += (conv_c_(vi)*vderxy_(1, 1) + v1*derxy_(1, vi))*v;
              }
            }
          }
#endif

          for (int ui=0; ui<nen_; ++ui)
          {
            const int tui  = 3*ui;
            const int tuip = tui+1;
            const double v = timetauM*densaf_*funct_(ui);
            const double v0 = v*gradp_(0);
            const double v1 = v*gradp_(1);
            for (int vi=0; vi<nen_; ++vi)
            {
              const int tvi  = 3*vi;
              const int tvip = tvi+1;
              /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
              /*
                /                                       \
                |         n+1    /                \     |
                |  nabla p    , |   rho*Du o nabla | v  |
                |         (i)    \                /     |
                \                                      /
              */
              estif(tvi,  tui ) += v0*derxy_(0, vi) ;
              estif(tvip, tui ) += v1*derxy_(0, vi) ;

              estif(tvi,  tuip) += v0*derxy_(1, vi) ;
              estif(tvip, tuip) += v1*derxy_(1, vi) ;

            }
          }

          if (is_higher_order_ele_)
          {
            for (int ui=0; ui<nen_; ++ui)
            {
              const int tui  = 3*ui;
              const int tuip = tui+1;
              const double v = -2.0*visceff_*timetauM*densaf_*funct_(ui);
              const double v0 = v*visc_old_(0);
              const double v1 = v*visc_old_(1);
              for (int vi=0; vi<nen_; ++vi)
              {
                const int tvi  = 3*vi;
                const int tvip = tvi+1;

                /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
                /*
                  /                                                 \
                  |               / n+1 \    /                \     |
                  |  nabla o eps | u     |, |  rho*Du o nabla | v   |
                  |               \ (i) /    \                /     |
                  \                                                 /
                */
                estif(tvi,  tui ) += v0*derxy_(0, vi) ;
                estif(tvip, tui ) += v1*derxy_(0, vi) ;

                estif(tvi,  tuip) += v0*derxy_(1, vi) ;
                estif(tvip, tuip) += v1*derxy_(1, vi) ;

              }
            }
          }

          for (int ui=0; ui<nen_; ++ui)
          {
            const int tui  = 3*ui;
            const int tuip = tui+1;
            const double v = -tau_M*densaf_*funct_(ui);
            const double v0 = v*rhsmom_(0);
            const double v1 = v*rhsmom_(1);
            for (int vi=0; vi<nen_; ++vi)
            {
              const int tvi  = 3*vi;
              const int tvip = tvi+1;

              /* supg stabilisation: bodyforce part, linearisation of test function */

              /*
                /                                       \
                |                 /                \     |
                |  rho*rhsint   , |  rho*Du o nabla | v  |
                |                 \                /     |
                \                                        /

              */
              estif(tvi , tui ) += v0*derxy_(0, vi) ;
              estif(tvip, tui ) += v1*derxy_(0, vi) ;

              estif(tvi , tuip) += v0*derxy_(1, vi) ;
              estif(tvip, tuip) += v1*derxy_(1, vi) ;

            } // vi
          } // ui
        } // if newton

#if 1
        // NOTE: Here we have a difference to the previous version of this
        // element!  Before we did not care for the mesh velocity in this
        // term. This seems unreasonable and wrong.
        for (int vi=0; vi<nen_; ++vi)
        {
          const int tvi = 3*vi;
          // supg stabilisation
          const double v = -tau_M*densaf_*(conv_c_(vi)+sgconv_c_(vi));
          eforce(tvi    ) += v*momres_old_(0) ;
          eforce(tvi + 1) += v*momres_old_(1) ;
        }
#endif
      }

      //----------------------------------------------------------------------
      //                       STABILISATION, VISCOUS PART

      if (is_higher_order_ele_)
      {
        if(vstab_ != Fluid3::viscous_stab_none)
        {
          const double two_visc_tauMp = vstabfac*2.0*visc_*tau_Mp;
          // viscous stabilization either on left hand side or on right hand side
          if (vstab_ == Fluid3::viscous_stab_gls || vstab_ == Fluid3::viscous_stab_usfem)
          {
            const double two_visc_timetauMp   = vstabfac*2.0*visc_*timetauMp;
            const double four_visc2_timetauMp = vstabfac*4.0*visceff_*visc_*timetauMp;

            // viscous stabilization on left hand side
            for (int ui=0; ui<nen_; ++ui)
            {
              const int tui  = 3*ui;
              const int tuip = tui+1;
              const double v = two_visc_tauMp*densam_*funct_(ui)
#if 1
                         + two_visc_timetauMp*densaf_*conv_c_(ui)
#endif
                               ;
              for (int vi=0; vi<nen_; ++vi)
              {
                const int tvi  = 3*vi;
                const int tvip = tvi+1;
                /* viscous stabilisation, inertia part */
                /*
                  /                          \
                  |                          |
              +/- |    rho*Du , div eps (v)  |
                  |                          |
                  \                          /
                */
                /* viscous stabilisation, convective part */
                /*
                  /                                        \
                  |  /       n+1       \                   |
              +/- | |   rho*u   o nabla | Du , div eps (v) |
                  |  \       (i)       /                   |
                  \                                        /
                */
                estif(tvi,  tui ) += v*viscs2_(0, vi) ;
                estif(tvip, tui ) += v*viscs2_(1, vi) ;

                estif(tvi,  tuip) += v*viscs2_(1, vi) ;
                estif(tvip, tuip) += v*viscs2_(3, vi) ;
              }
            }

            for (int ui=0; ui<nen_; ++ui)
            {
              const int tuipp = 3*ui + 2;
              for (int vi=0; vi<nen_; ++vi)
              {
                const int tvi = 3*vi;
                /* viscous stabilisation, pressure part ( L_pres_p) */
                /*
                  /                          \
                  |                          |
             +/-  |  nabla Dp , div eps (v)  |
                  |                          |
                  \                          /
                */
                estif(tvi,     tuipp) += two_visc_timetauMp*(derxy_(0, ui)*viscs2_(0, vi)
                                                            +
                                                            derxy_(1, ui)*viscs2_(1, vi)) ;
                estif(tvi + 1, tuipp) += two_visc_timetauMp*(derxy_(0, ui)*viscs2_(1, vi)
                                                            +
                                                            derxy_(1, ui)*viscs2_(3, vi)) ;

              }
            }

            for (int ui=0; ui<nen_; ++ui)
            {
              const int tui  = 3*ui;
              const int tuip = tui+1;
              for (int vi=0; vi<nen_; ++vi)
              {
                const int tvi  = 3*vi;
                const int tvip = tvi+1;
                /* viscous stabilisation, viscous part (-L_visc_u) */
                /*
                  /                                 \
                  |               /  \                |
             -/+  |  nabla o eps | Du | , div eps (v) |
                  |               \  /                |
                  \                                 /
                */
                estif(tvi,  tui ) -= four_visc2_timetauMp*(viscs2_(0,ui)*viscs2_(0,vi)+viscs2_(1,ui)*viscs2_(1,vi)) ;
                estif(tvip, tui ) -= four_visc2_timetauMp*(viscs2_(0,ui)*viscs2_(1,vi)+viscs2_(1,ui)*viscs2_(3,vi)) ;

                estif(tvi,  tuip) -= four_visc2_timetauMp*(viscs2_(0,vi)*viscs2_(1,ui)+viscs2_(1,vi)*viscs2_(3,ui)) ;
                estif(tvip, tuip) -= four_visc2_timetauMp*(viscs2_(1,ui)*viscs2_(1,vi)+viscs2_(3,ui)*viscs2_(3,vi)) ;
              } // vi
            } // ui

            if (is_newton_)
            {
              for (int ui=0; ui<nen_; ++ui)
              {
                const int tui  = 3*ui;
                const int tuip = tui+1;
                const double v = two_visc_timetauMp*densaf_*funct_(ui);
                for (int vi=0; vi<nen_; ++vi)
                {
                  const int tvi  = 3*vi;
                  const int tvip = tvi+1;
                  /* viscous stabilisation, reactive part of convection */
                  /*
                    /                                         \
                    |  /                \   n+1               |
                +/- | |   rho*Du o nabla | u    , div eps (v) |
                    |  \                /   (i)               |
                    \                                         /
                  */
                  estif(tvi,  tui ) += v*(viscs2_(0,vi)*vderxy_(0,0)+viscs2_(1,vi)*vderxy_(1,0)) ;
                  estif(tvip, tui ) += v*(viscs2_(1,vi)*vderxy_(0,0)+viscs2_(3,vi)*vderxy_(1,0)) ;

                  estif(tvi,  tuip) += v*(viscs2_(0,vi)*vderxy_(0,1)+viscs2_(1,vi)*vderxy_(1,1)) ;
                  estif(tvip, tuip) += v*(viscs2_(1,vi)*vderxy_(0,1)+viscs2_(3,vi)*vderxy_(1,1)) ;
                } // vi
              } // ui
            } // if newton
          } // end if viscous stabilization on left hand side

          for (int vi=0; vi<nen_; ++vi)
          {
            const int tvi = 3*vi;
            /* viscous stabilisation */
            eforce(tvi    ) -= two_visc_tauMp*(momres_old_(0)*viscs2_(0, vi)+momres_old_(1)*viscs2_(1, vi)) ;
            eforce(tvi + 1) -= two_visc_tauMp*(momres_old_(0)*viscs2_(1, vi)+momres_old_(1)*viscs2_(3, vi)) ;
          }
        }
      }

      //----------------------------------------------------------------------
      //                     STABILISATION, CONTINUITY PART

      if (cstab_ == Fluid3::continuity_stab_yes)
      {
        const double timetauC = timefac*tau_C;
        for (int ui=0; ui<nen_; ++ui)
        {
          const int tui  = 3*ui;
          const int tuip = tui+1;
          const double v0 = timetauC*derxy_(0, ui);
          const double v1 = timetauC*derxy_(1, ui);
          for (int vi=0; vi<nen_; ++vi)
          {
            const int tvi  = 3*vi;
            const int tvip = tvi+1;
            /* continuity stabilisation on left hand side */
            /*
              /                         \
              |                          |
              | nabla o Du  , nabla o v  |
              |                          |
              \                         /
            */
            estif(tvi,  tui ) += v0*derxy_(0, vi) ;
            estif(tvip, tui ) += v0*derxy_(1, vi) ;

            estif(tvi,  tuip) += v1*derxy_(0, vi) ;
            estif(tvip, tuip) += v1*derxy_(1, vi) ;
          }
        }

        const double tauC_conres = tau_C*conres_old_;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int tvi = 3*vi;
          /* continuity stabilisation on right hand side */
          eforce(tvi    ) -= tauC_conres*derxy_(0, vi) ;
          eforce(tvi + 1) -= tauC_conres*derxy_(1, vi) ;
        }
      }

      //----------------------------------------------------------------------
      //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)

      if(fssgv_ != Fluid3::no_fssgv)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int tvi = 3*vi;
          /* fine-scale subgrid-viscosity term on right hand side */
          /*
                              /                          \
                             |       /    \         / \   |
             - mu_art(fsu) * |  eps | Dfsu | , eps | v |  |
                             |       \    /         \ /   |
                              \                          /
          */
          eforce(tvi    ) -= fssgviscfac*(2.0*derxy_(0, vi)*fsvderxy_(0, 0)
                                         +    derxy_(1, vi)*fsvderxy_(0, 1)
                                         +    derxy_(1, vi)*fsvderxy_(1, 0)) ;
          eforce(tvi + 1) -= fssgviscfac*(    derxy_(0, vi)*fsvderxy_(0, 1)
                                         +    derxy_(0, vi)*fsvderxy_(1, 0)
                                         +2.0*derxy_(1, vi)*fsvderxy_(1, 1)) ;
        }
      }
    }

    // linearization with respect to mesh motion
    if (emesh.IsInitialized())
    {

      // xGderiv_ = sum(gridx(k,i) * deriv_(j,k), k);
      // xGderiv_ == xjm_

      // mass + rhs
      for (int vi=0; vi<nen_; ++vi)
      {
        const int tvi   = 3*vi;
        const int tvip  = tvi + 1;

        const double v = fac_*funct_(vi);
        for (int ui=0; ui<nen_; ++ui)
        {
          const int tui   = 3*ui;
          const int tuip  = tui + 1;

          emesh(tvi,   tui ) += v*(velint_(0)-rhsmom_(0))*derxy_(0, ui);
          emesh(tvi,   tuip) += v*(velint_(0)-rhsmom_(0))*derxy_(1, ui);

          emesh(tvip,  tui ) += v*(velint_(1)-rhsmom_(1))*derxy_(0, ui);
          emesh(tvip,  tuip) += v*(velint_(1)-rhsmom_(1))*derxy_(1, ui);
        }
      }

      vderiv_.MultiplyNT(evelaf, deriv_);

//#define derxjm_(r,c,d,i) derxjm_ ## r ## c ## d (i)

//#define derxjm_001(ui) (deriv_(2, ui)*xjm_(1, 2) - deriv_(1, ui)*xjm_(2, 2))
//#define derxjm_100(ui) (deriv_(1, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(1, 2))
//#define derxjm_011(ui) (deriv_(0, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(0, 2))
//#define derxjm_110(ui) (deriv_(2, ui)*xjm_(0, 2) - deriv_(0, ui)*xjm_(2, 2))

      for (int vi=0; vi<nen_; ++vi)
      {
        const int tvi  = 3*vi;
        const int tvip = tvi+1;
        const double v = timefacfac/det_*funct_(vi);
        for (int ui=0; ui<nen_; ++ui)
        {
          const int tui  = 3*ui;
          const int tuip = tui+1;

          emesh(tvi , tui ) += v*(
          + convvelint_(1)*(-vderiv_(0, 0)*deriv_(1,ui) + vderiv_(0, 1)*deriv_(0,ui))
          );

          emesh(tvi , tuip) += v*(
          + convvelint_(0)*(-vderiv_(0, 0)*deriv_(1,ui) + vderiv_(0, 1)*deriv_(0,ui))
          );

          emesh(tvip, tui ) += v*(
          + convvelint_(1)*(-vderiv_(1, 0)*deriv_(1,ui) + vderiv_(1, 1)*deriv_(0,ui))
          );

          emesh(tvip, tuip) += v*(
          + convvelint_(0)*(-vderiv_(1, 0)*deriv_(1,ui) + vderiv_(1, 1)*deriv_(0,ui))
          );
        }
      }

      // pressure
      for (int vi=0; vi<nen_; ++vi)
      {
        const int tvi  = 3*vi;
        const int tvip = tvi+1;
        const double v = press*timefacfac/det_;
        for (int ui=0; ui<nen_; ++ui)
        {
          const int tui = 3*ui;
          emesh(tvi,  tui + 1) += v*(deriv_(0, vi)*deriv_(1, ui) - deriv_(0, ui)*deriv_(1, vi)) ;
          emesh(tvip, tui    ) += v*(deriv_(0, vi)*deriv_(1, ui) - deriv_(0, ui)*deriv_(1, vi)) ;
        }
      }

      // div u
      for (int vi=0; vi<nen_; ++vi)
      {
        const int tvipp = 3*vi + 2;
        const double v = timefacfac/det_*funct_(vi);
        for (int ui=0; ui<nen_; ++ui)
        {
          const int tui = 3*ui;
          emesh(tvipp, tui) += v*(
          deriv_(0,ui)*vderiv_(1,1) - deriv_(1,ui)*vderiv_(1,0)
          ) ;

          emesh(tvipp, tui + 1) += v*(
          deriv_(0,ui)*vderiv_(0,1) - deriv_(1,ui)*vderiv_(0,0)
          ) ;
        }
      }
    }
  } // loop gausspoints
}


/*----------------------------------------------------------------------*
 |  calculation of (all-scale) subgrid viscosity               vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::CalcSubgrVisc_2D(
  const LINALG::Matrix<nsd_,nen_>&        evelaf,
  const double                        area,
  double&                             Cs
  )
{
  // get characteristic element length for Smagorinsky model
  const double hk = sqrt(area);

  //
  // SMAGORINSKY MODEL
  // -----------------
  //                                   +-                                 -+ 1
  //                               2   |          / h \           / h \    | -
  //    visc          = dens * lmix  * | 2 * eps | u   |   * eps | u   |   | 2
  //        turbulent           |      |          \   / ij        \   / ij |
  //                            |      +-                                 -+
  //                            |
  //                            |      |                                   |
  //                            |      +-----------------------------------+
  //                            |           'resolved' rate of strain
  //                    mixing length
  // -> either provided by dynamic modeling procedure and stored in Cs_delta_sq
  // -> or computed based on fixed Smagorinsky constant Cs:
  //             Cs = 0.17   (Lilly --- Determined from filter
  //                          analysis of Kolmogorov spectrum of
  //                          isotropic turbulence)
  //             0.1 < Cs < 0.24 (depending on the flow)
  //

  // compute (all-scale) rate of strain
  double rateofstrain = -1.0e30;
  rateofstrain = GetStrainRate(evelaf,derxy_,vderxy_);

  // mixing length set proportional to element length: lmix = Cs * hk
  double lmix = Cs * hk;

  // subgrid viscosity
  sgvisc_ = densaf_ * lmix * lmix * rateofstrain;

  return;
}

/*----------------------------------------------------------------------*
 |  calculation of fine-scale subgrid viscosity                vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::CalcFineScaleSubgrVisc_2D(
  const LINALG::Matrix<nsd_,nen_>&       evelaf,
  const LINALG::Matrix<nsd_,nen_>&       fsevelaf,
  const double                       area,
  double&                            Cs
  )
{
  // get characteristic element length for Smagorinsky model
  const double hk = sqrt(area);

  if (fssgv_ == Fluid3::smagorinsky_all)
  {
    //
    // ALL-SCALE SMAGORINSKY MODEL
    // ---------------------------
    //                                      +-                                 -+ 1
    //                                  2   |          / h \           / h \    | -
    //    visc          = dens * (C_S*h)  * | 2 * eps | u   |   * eps | u   |   | 2
    //        turbulent                     |          \   / ij        \   / ij |
    //                                      +-                                 -+
    //                                      |                                   |
    //                                      +-----------------------------------+
    //                                            'resolved' rate of strain
    //

    // compute (all-scale) rate of strain
    double rateofstrain = -1.0e30;
    rateofstrain = GetStrainRate(evelaf,derxy_,vderxy_);

    fssgvisc_ = densaf_ * Cs * Cs * hk * hk * rateofstrain;
  }
  else if (fssgv_ == Fluid3::smagorinsky_small)
  {
    //
    // FINE-SCALE SMAGORINSKY MODEL
    // ----------------------------
    //                                      +-                                 -+ 1
    //                                  2   |          /    \          /   \    | -
    //    visc          = dens * (C_S*h)  * | 2 * eps | fsu |   * eps | fsu |   | 2
    //        turbulent                     |          \   / ij        \   / ij |
    //                                      +-                                 -+
    //                                      |                                   |
    //                                      +-----------------------------------+
    //                                            'resolved' rate of strain
    //

    // fine-scale rate of strain
    double fsrateofstrain = -1.0e30;
    fsrateofstrain = GetStrainRate(fsevelaf,derxy_,fsvderxy_);

    fssgvisc_ = densaf_ * Cs * Cs * hk * hk * fsrateofstrain;
  }

  return;
}

/*----------------------------------------------------------------------*
 |  calculation of stabilization parameter                     vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::CalcStabParameter_2D(
  const double               timefac,
  const double               area
  )
{
  // get characteristic element length
  const double hk = sqrt(area);

  // get element-type constant for tau
  double mk=0.0;
  switch (distype)
  {
    case DRT::Element::tri3:
    case DRT::Element::quad4:
      mk = 0.333333333333333333333;
      break;
    case DRT::Element::quad8:
    case DRT::Element::quad9:
    case DRT::Element::tri6:
      mk = 0.083333333333333333333;
      break;
    default:
      dserror("type unknown!\n");
  }

  // get velocity norm
  const double vel_norm = velint_.Norm2();

  // normed velocity at element centre (currently not used)
  //if (vel_norm>=1e-6) velino_ = velint_/vel_norm;
  //else
  //{
  //  velino_.Clear();
  //  velino_(0) = 1;
  //}

  // get streamlength (currently not used)
  //const double val = blitz::sum(blitz::abs(blitz::sum(velino_(j)*derxy_(j,i),j)));
  //const double strle = 2.0/val;

  // ---------------------------------------------------------------
  // computation of stabilization parameter tau
  // ---------------------------------------------------------------
  if (whichtau_ == Fluid3::franca_barrenechea_valentin_wall)
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
    const double re1 = 4.0 * timefac * visceff_ / (mk * densaf_ * DSQR(hk));

    /* convective : viscous forces */
    const double re2 = mk * densaf_ * vel_norm * hk / (2.0 * visceff_);

    const double xi1 = DMAX(re1,1.0);
    const double xi2 = DMAX(re2,1.0);

    tau_(0) = timefac*DSQR(hk)/(DSQR(hk)*densaf_*xi1+(4.0*timefac*visceff_/mk)*xi2);
    tau_(1) = tau_(0);

    /*------------------------------------------------------ compute tau_C ---*/
    // PhD thesis Wall (1999)
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
    const double xi_tau_c = DMIN(re2,1.0);
    tau_(2) = densaf_ * vel_norm * hk * 0.5 * xi_tau_c;
  }
  else if(whichtau_ == Fluid3::bazilevs)
  {
    /*

    tau_M: Bazilevs et al. (2007)
                                                                              1.0
                 +-                                                      -+ - ---
                 |        2                                               |   2.0
                 | 4.0*rho         n+1             n+1          2         |
          tau  = | -------  + rho*u     * G * rho*u     + C * mu  * G : G |
             M   |     2                  -                I        -   - |
                 |   dt                   -                         -   - |
                 +-                                                      -+

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
    /*                                 +----
               n+1             n+1     \         n+1              n+1
          rho*u     * G * rho*u     =   +   rho*u    * G   * rho*u
                      -                /         i     -ij        j
                      -               +----        -
                                        i,j
    */
    double G;
    double normG = 0;
    double Gnormu = 0;

    const double dens_sqr = densaf_*densaf_;
    for (int nn=0;nn<2;++nn)
    {
      const double dens_sqr_velint_nn = dens_sqr*velint_(nn);
      for (int rr=0;rr<2;++rr)
      {
        G = xji_(nn,0)*xji_(rr,0);
        for (int mm=1;mm<2;++mm)
        {
          G += xji_(nn,mm)*xji_(rr,mm);
        }
        normG += G*G;
        Gnormu+=dens_sqr_velint_nn*G*velint_(rr);
      }
    }

    // definition of constant:
    // 12.0/m_k = 36.0 for linear elements and 144.0 for quadratic elements
    // (differently defined, e.g., in Akkerman et al. (2008))
    const double CI = 12.0/mk;

    tau_(0) = 1.0/(sqrt((4.0*dens_sqr)/(dt_*dt_)+Gnormu+CI*visceff_*visceff_*normG));
    tau_(1) = tau_(0);

    /*
      tau_C: Bazilevs et al. (2007), derived from fine-scale complement Shur
                                     operator of the pressure equation

                                  1.0
                    tau  = -----------------
                       C            /     \
                            tau  * | g * g |
                               M    \-   -/
    */
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
    const double g0 = xji_(0,0) + xji_(0,1);
    const double g1 = xji_(1,0) + xji_(1,1);

    const double normgsq = g0*g0+g1*g1;

    tau_(2) = 1.0/(tau_(0)*normgsq);
  }
  else if(whichtau_ == Fluid3::codina)
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
    const double re1 = 4.0 * timefac * visceff_ / (mk * densaf_ * DSQR(hk));

    /* convective : viscous forces */
    const double re2 = mk * densaf_ * vel_norm * hk / (2.0 * visceff_);

    const double xi1 = DMAX(re1,1.0);
    const double xi2 = DMAX(re2,1.0);

    tau_(0) = timefac*DSQR(hk)/(DSQR(hk)*densaf_*xi1+(4.0*timefac*visceff_/mk)*xi2);
    tau_(1) = tau_(0);

    /*------------------------------------------------------ compute tau_C ---*/
    /*-- stability parameter definition according to Codina (2002), CMAME 191
     *
     * Analysis of a stabilized finite element approximation of the transient
     * convection-diffusion-reaction equation using orthogonal subscales.
     * Ramon Codina, Jordi Blasco; Comput. Visual. Sci., 4 (3): 167-174, 2002.
     *
     * */
    tau_(2) = sqrt(DSQR(visceff_)+DSQR(0.5*densaf_*vel_norm*hk));
  }
  else dserror("unknown definition of tau\n");

  return;
}


#endif
#endif
