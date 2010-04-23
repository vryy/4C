/*----------------------------------------------------------------------*/
/*!
\file fluid3_impl.cpp

\brief Internal implementation of Fluid3 element

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
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

#include "../drt_geometry/position_array.H"

#include <Epetra_SerialDenseSolver.h>

#ifdef DEBUG
#endif

/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)   g.bau 03/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::Sysmat2D3D(
  Fluid3*                                     ele,
  const LINALG::Matrix<nsd_,nen_>&             evelaf,
  const LINALG::Matrix<nsd_,nen_>&             eveln,
  const LINALG::Matrix<nsd_,nen_>&             fsevelaf,
  const LINALG::Matrix<nen_,1>&                epreaf,
  const LINALG::Matrix<nsd_,nen_>&             eaccam,
  const LINALG::Matrix<nen_,1>&                escaaf,
  const LINALG::Matrix<nen_,1>&                escaam,
  const LINALG::Matrix<nen_,1>&                escadtam,
  const LINALG::Matrix<nsd_,nen_>&             emhist,
  const LINALG::Matrix<nsd_,nen_>&             edispnp,
  const LINALG::Matrix<nsd_,nen_>&             egridv,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  estif,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  emesh,
  LINALG::Matrix<(nsd_+1)*nen_,1>&             eforce,
  const double                                 thermpressaf,
  const double                                 thermpressam,
  const double                                 thermpressdtam,
  Teuchos::RCP<const MAT::Material>           material,
  const double                                timefac,
  double&                                     Cs,
  double&                                     Cs_delta_sq,
  double&                                     l_tau
  )
{
  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  /*
  DRT::Node** nodes = ele->Nodes();
  for (int inode=0; inode<nen_; inode++)
  {
    const double* x = nodes[inode]->X();
    xyze_(0,inode) = x[0];
    xyze_(1,inode) = x[1];
    xyze_(2,inode) = x[2];
  }
*/
  // add displacement when fluid nodes move in the ALE case
  if (ele->IsAle()) xyze_ += edispnp;

  // ---------------------------------------------------------------------
  // call routine for calculation of body force in element nodes
  // (time n+alpha_F for generalized-alpha scheme, at time n+1 otherwise)
  // ---------------------------------------------------------------------
  BodyForce(ele);

  // evaluate shape functions and derivatives at element center
  EvalShapeFuncAndDerivsAtEleCenter(ele->Id());

  // element aera or volume
  const double vol = fac_;

  // in case of viscous stabilization decide whether to use GLS or USFEM
  double vstabfac= 0.0;
  if (vstab_ == INPAR::FLUID::viscous_stab_usfem or
      vstab_ == INPAR::FLUID::viscous_stab_usfem_only_rhs)   vstabfac =  1.0;
  else if(vstab_ == INPAR::FLUID::viscous_stab_gls or
          vstab_ == INPAR::FLUID::viscous_stab_gls_only_rhs) vstabfac = -1.0;

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
    if (turb_mod_action_ != Fluid3::no_model)
    {
      CalcSubgrVisc(evelaf,vol,Cs,Cs_delta_sq,l_tau);

      // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
      visceff_ += sgvisc_;
    }
    else if (fssgv_ != Fluid3::no_fssgv)
      CalcFineScaleSubgrVisc(evelaf,fsevelaf,vol,Cs);

    // get velocity at element center
    velint_.Multiply(evelaf,funct_);

    // ---------------------------------------------------------------------
    // calculate stabilization parameter at element center
    // ---------------------------------------------------------------------
    CalcStabParameter(timefac,vol);
  }

  // Gaussian integration points
  //const DRT::UTILS::IntegrationPoints3D intpoints(ele->gaussrule_);
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
      if (turb_mod_action_ != Fluid3::no_model)
      {
        CalcSubgrVisc(evelaf,vol,Cs,Cs_delta_sq,l_tau);

        // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
        visceff_ += sgvisc_;
      }
      else if (fssgv_ != Fluid3::no_fssgv)
        CalcFineScaleSubgrVisc(evelaf,fsevelaf,vol,Cs);

      // ---------------------------------------------------------------------
      // calculate stabilization parameter at element center
      // ---------------------------------------------------------------------
      CalcStabParameter(timefac,vol);
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
    if (ele->IsAle()) convvelint_.Multiply(-1.0, egridv, funct_, 1.0);

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

      double prefac;
      if(physicaltype_ == INPAR::FLUID::loma)
      {
        prefac = 1.0/3.0;
        derxy2_.Scale(prefac);
      }
      else prefac = 1.0;

      double sum = (derxy2_(0,0)+derxy2_(1,0)+derxy2_(2,0))/prefac;

      viscs2_(0,0) = 0.5 * (sum + derxy2_(0,0));
      viscs2_(1,0) = 0.5 *  derxy2_(3,0);
      viscs2_(2,0) = 0.5 *  derxy2_(4,0);
      viscs2_(3,0) = 0.5 *  derxy2_(3,0);
      viscs2_(4,0) = 0.5 * (sum + derxy2_(1,0));
      viscs2_(5,0) = 0.5 *  derxy2_(5,0);
      viscs2_(6,0) = 0.5 *  derxy2_(4,0);
      viscs2_(7,0) = 0.5 *  derxy2_(5,0);
      viscs2_(8,0) = 0.5 * (sum + derxy2_(2,0));

      visc_old_(0) = viscs2_(0,0)*evelaf(0,0)+viscs2_(1,0)*evelaf(1,0)+viscs2_(2,0)*evelaf(2,0);
      visc_old_(1) = viscs2_(3,0)*evelaf(0,0)+viscs2_(4,0)*evelaf(1,0)+viscs2_(5,0)*evelaf(2,0);
      visc_old_(2) = viscs2_(6,0)*evelaf(0,0)+viscs2_(7,0)*evelaf(1,0)+viscs2_(8,0)*evelaf(2,0);

      for (int i=1; i<nen_; ++i)
      {
        double sum = (derxy2_(0,i)+derxy2_(1,i)+derxy2_(2,i))/prefac;

        viscs2_(0,i) = 0.5 * (sum + derxy2_(0,i));
        viscs2_(1,i) = 0.5 *  derxy2_(3,i);
        viscs2_(2,i) = 0.5 *  derxy2_(4,i);
        viscs2_(3,i) = 0.5 *  derxy2_(3,i);
        viscs2_(4,i) = 0.5 * (sum + derxy2_(1,i));
        viscs2_(5,i) = 0.5 *  derxy2_(5,i);
        viscs2_(6,i) = 0.5 *  derxy2_(4,i);
        viscs2_(7,i) = 0.5 *  derxy2_(5,i);
        viscs2_(8,i) = 0.5 * (sum + derxy2_(2,i));

        visc_old_(0) += viscs2_(0,i)*evelaf(0,i)+viscs2_(1,i)*evelaf(1,i)+viscs2_(2,i)*evelaf(2,i);
        visc_old_(1) += viscs2_(3,i)*evelaf(0,i)+viscs2_(4,i)*evelaf(1,i)+viscs2_(5,i)*evelaf(2,i);
        visc_old_(2) += viscs2_(6,i)*evelaf(0,i)+viscs2_(7,i)*evelaf(1,i)+viscs2_(8,i)*evelaf(2,i);
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
    vdiv_ = vderxy_(0, 0) + vderxy_(1, 1) + vderxy_(2, 2);

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
      for (int rr=0;rr<3;++rr)
      {
        momres_old_(rr) = densam_*accint_(rr)+densaf_*conv_old_(rr)+gradp_(rr)-2*visceff_*visc_old_(rr)-densaf_*bodyforce_(rr);
      }

      if (cross_    != INPAR::FLUID::cross_stress_stab_none or
          reynolds_ != INPAR::FLUID::reynolds_stress_stab_none)
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

      if(physicaltype_ == INPAR::FLUID::loma)
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
      // else: 										f = rho * g

      if (physicaltype_ == INPAR::FLUID::boussinesq)
        rhsmom_.Update(densn_,histmom_,deltadens_*timefac,bodyforce_);
      else rhsmom_.Update(densn_,histmom_,densaf_*timefac,bodyforce_);

      // modify integration factor for Galerkin rhs
      rhsfac *= timefac;

      // evaluate momentum residual once for all stabilization right hand sides
      if (is_stationary_ == false)
      {
        // momres_old = u_(n+1) + theta ( ... ) - histmom_ - bodyforce_
        for (int rr=0;rr<3;++rr)
        {
          momres_old_(rr) = densaf_*velint_(rr)+timefac*(densaf_*conv_old_(rr)+gradp_(rr)-2*visceff_*visc_old_(rr))-rhsmom_(rr);
        }
      }
      else
      {
        // momres_old = theta ( ... ) -  bodyforce_
        for (int rr=0;rr<3;++rr)
        {
          momres_old_(rr) = timefac*(densaf_*conv_old_(rr)+gradp_(rr)-2*visceff_*visc_old_(rr))-rhsmom_(rr);
        }
      }

      if (cross_    != INPAR::FLUID::cross_stress_stab_none or
          reynolds_ != INPAR::FLUID::reynolds_stress_stab_none)
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
        const double vdivn = vderxyn_(0, 0) + vderxyn_(1, 1) + vderxyn_(2, 2);

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

      // // The computation of the inertia and the convective term is split
      if (is_stationary_ == false)
      {
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;

          const double v = fac_*densam_*funct_(ui);

          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            /* inertia (contribution to mass matrix) */
            /*
            /                \
            |                |
            |    rho*Du , v  |
            |                |
            \                /
            */

            double v2 = v*funct_(vi) ;
            estif(fvi  , fui  ) += v2;
            estif(fvip , fuip ) += v2;
            estif(fvipp, fuipp) += v2;
          }
        }
      }

      // convective part
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui   = 4*ui;
        const int fuip  = fui+1;
        const int fuipp = fui+2;

        const double v = timefacfac*densaf_*(conv_c_(ui)+sgconv_c_(ui));

        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;

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
          estif(fvipp, fuipp) += v2;
        }
      }

      if (is_newton_)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;
          const double v = timefacfac*densaf_*funct_(vi);
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
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
            estif(fvi,   fuipp) += v2*vderxy_(0, 2) ;
            estif(fvip,  fui)   += v2*vderxy_(1, 0) ;
            estif(fvip,  fuip)  += v2*vderxy_(1, 1) ;
            estif(fvip,  fuipp) += v2*vderxy_(1, 2) ;
            estif(fvipp, fui)   += v2*vderxy_(2, 0) ;
            estif(fvipp, fuip)  += v2*vderxy_(2, 1) ;
            estif(fvipp, fuipp) += v2*vderxy_(2, 2) ;
          }
        }
      }

      // inertia terms on the right hand side for instationary fluids
      if (is_stationary_ == false)
      {
        if (is_genalpha_)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi = 4*vi;
            /* inertia term on right-hand side for generalized-alpha scheme */
            const double v = -fac_*densam_*funct_(vi);
            eforce(fvi    ) += v*accint_(0) ;
            eforce(fvi + 1) += v*accint_(1) ;
            eforce(fvi + 2) += v*accint_(2) ;
          }
        }
        else
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi = 4*vi;
            /* inertia term on right-hand side for one-step-theta/BDF2 schem */
            const double v = -fac_*densaf_*funct_(vi);
            eforce(fvi    ) += v*velint_(0) ;
            eforce(fvi + 1) += v*velint_(1) ;
            eforce(fvi + 2) += v*velint_(2) ;
          }
        }
      }

#if 1
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi   = 4*vi;
        /* convection (convective form) on right-hand side */
        double v = -rhsfac*densaf_*funct_(vi);
        eforce(fvi    ) += v*conv_old_(0) ;
        eforce(fvi + 1) += v*conv_old_(1) ;
        eforce(fvi + 2) += v*conv_old_(2) ;
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
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          double v = timefacfac*densaf_*funct_(ui)*vdiv_;
          if (physicaltype_ == INPAR::FLUID::loma) v -= timefacfac*densaf_*scaconvfacaf_*conv_scaaf_;
          // only with linear density-concentration correlation
          else if(physicaltype_ == INPAR::FLUID::varying_density) v += timefacfac*conv_scaaf_;
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
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
            estif(fvipp, fuipp) += v2;
          }
        }

        if (is_newton_)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            const double v0 = timefacfac*densaf_*velint_(0)*funct_(vi);
            const double v1 = timefacfac*densaf_*velint_(1)*funct_(vi);
            const double v2 = timefacfac*densaf_*velint_(2)*funct_(vi);
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
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
              estif(fvi,  fuipp) += v0*derxy_(2, ui) ;
              estif(fvip, fui  ) += v1*derxy_(0, ui) ;
              estif(fvip, fuip ) += v1*derxy_(1, ui) ;
              estif(fvip, fuipp) += v1*derxy_(2, ui) ;
              estif(fvipp,fui  ) += v2*derxy_(0, ui) ;
              estif(fvipp,fuip ) += v2*derxy_(1, ui) ;
              estif(fvipp,fuipp) += v2*derxy_(2, ui) ;
            }
          }

          if (physicaltype_ == INPAR::FLUID::loma)
          {
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = 4*vi;
              const int fvip  = fvi+1;
              const int fvipp = fvi+2;
              const double v0 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(0)*velint_(0)*funct_(vi);
              const double v1 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(1)*velint_(1)*funct_(vi);
              const double v2 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(2)*velint_(2)*funct_(vi);
              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui   = 4*ui;
                const int fuip  = fui+1;
                const int fuipp = fui+2;
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
                estif(fvi,  fuipp) += v0*funct_(ui) ;
                estif(fvip, fui  ) += v1*funct_(ui) ;
                estif(fvip, fuip ) += v1*funct_(ui) ;
                estif(fvip, fuipp) += v1*funct_(ui) ;
                estif(fvipp,fui  ) += v2*funct_(ui) ;
                estif(fvipp,fuip ) += v2*funct_(ui) ;
                estif(fvipp,fuipp) += v2*funct_(ui) ;
              }
            }
          }
          if (physicaltype_ == INPAR::FLUID::varying_density)
          {
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = 4*vi;
              const int fvip  = fvi+1;
              const int fvipp = fvi+2;
              const double v0 = +timefacfac*grad_scaaf_(0)*velint_(0)*funct_(vi);
              const double v1 = +timefacfac*grad_scaaf_(1)*velint_(1)*funct_(vi);
              const double v2 = +timefacfac*grad_scaaf_(2)*velint_(2)*funct_(vi);
              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui   = 4*ui;
                const int fuip  = fui+1;
                const int fuipp = fui+2;
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
                estif(fvi,  fuipp) += v0*funct_(ui) ;
                estif(fvip, fui  ) += v1*funct_(ui) ;
                estif(fvip, fuip ) += v1*funct_(ui) ;
                estif(fvip, fuipp) += v1*funct_(ui) ;
                estif(fvipp,fui  ) += v2*funct_(ui) ;
                estif(fvipp,fuip ) += v2*funct_(ui) ;
                estif(fvipp,fuipp) += v2*funct_(ui) ;
              }
            }
          }
        }

        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 4*vi;
          /* convection (conservative addition) on right-hand side */
          double v = -rhsfac*densaf_*funct_(vi)*vdiv_;
          eforce(fvi    ) += v*velint_(0) ;
          eforce(fvi + 1) += v*velint_(1) ;
          eforce(fvi + 2) += v*velint_(2) ;
        }

        if (physicaltype_ == INPAR::FLUID::loma)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            /* convection (conservative addition) on rhs for low-Mach-number flow */
            double v = rhsfac*densaf_*scaconvfacaf_*conv_scaaf_*funct_(vi);
            eforce(fvi    ) += v*velint_(0) ;
            eforce(fvi + 1) += v*velint_(1) ;
            eforce(fvi + 2) += v*velint_(2) ;
          }
        }
        if (physicaltype_ == INPAR::FLUID::varying_density)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            /* convection (conservative addition) on rhs for low-Mach-number flow */
            double v = -rhsfac*conv_scaaf_*funct_(vi);
            eforce(fvi    ) += v*velint_(0) ;
            eforce(fvi + 1) += v*velint_(1) ;
            eforce(fvi + 2) += v*velint_(2) ;
          }
        }
      }

      //----------------------------------------------------------------------
      // computation of viscosity term including right-hand-side contribution
      //----------------------------------------------------------------------
      const double visceff_timefacfac = visceff_*timefacfac;
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui   = 4*ui;
        const int fuip  = fui+1;
        const int fuipp = fui+2;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;

          const double derxy_0ui_0vi = derxy_(0, ui)*derxy_(0, vi);
          const double derxy_1ui_1vi = derxy_(1, ui)*derxy_(1, vi);
          const double derxy_2ui_2vi = derxy_(2, ui)*derxy_(2, vi);
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
                                                     derxy_1ui_1vi
                                                     +
                                                     derxy_2ui_2vi) ;
          estif(fvi , fuip)   += visceff_timefacfac*derxy_(0, ui)*derxy_(1, vi) ;
          estif(fvi , fuipp)  += visceff_timefacfac*derxy_(0, ui)*derxy_(2, vi) ;
          estif(fvip, fui)    += visceff_timefacfac*derxy_(0, vi)*derxy_(1, ui) ;
          estif(fvip, fuip)   += visceff_timefacfac*(derxy_0ui_0vi
                                                     +
                                                     2.0*derxy_1ui_1vi
                                                     +
                                                     derxy_2ui_2vi) ;
          estif(fvip , fuipp) += visceff_timefacfac*derxy_(1, ui)*derxy_(2, vi) ;
          estif(fvipp, fui)   += visceff_timefacfac*derxy_(0, vi)*derxy_(2, ui) ;
          estif(fvipp, fuip)  += visceff_timefacfac*derxy_(1, vi)*derxy_(2, ui) ;
          estif(fvipp, fuipp) += visceff_timefacfac*(derxy_0ui_0vi
                                                     +
                                                     derxy_1ui_1vi
                                                     +
                                                     2.0*derxy_2ui_2vi) ;

        }
      }

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = 4*vi;
        const double v = -visceff_*rhsfac;
        /* viscosity term on right-hand side */
        eforce(fvi)     += v*(2.0*derxy_(0, vi)*vderxy_(0, 0)
                              +
                              derxy_(1, vi)*vderxy_(0, 1)
                              +
                              derxy_(1, vi)*vderxy_(1, 0)
                              +
                              derxy_(2, vi)*vderxy_(0, 2)
                              +
                              derxy_(2, vi)*vderxy_(2, 0)) ;
        eforce(fvi + 1) += v*(derxy_(0, vi)*vderxy_(0, 1)
                              +
                              derxy_(0, vi)*vderxy_(1, 0)
                              +
                              2.0*derxy_(1, vi)*vderxy_(1, 1)
                              +
                              derxy_(2, vi)*vderxy_(1, 2)
                              +
                              derxy_(2, vi)*vderxy_(2, 1)) ;
        eforce(fvi + 2) += v*(derxy_(0, vi)*vderxy_(0, 2)
                              +
                              derxy_(0, vi)*vderxy_(2, 0)
                              +
                              derxy_(1, vi)*vderxy_(1, 2)
                              +
                              derxy_(1, vi)*vderxy_(2, 1)
                              +
                              2.0*derxy_(2, vi)*vderxy_(2, 2)) ;
      }

      //----------------------------------------------------------------------
      // computation of pressure term including right-hand-side contribution
      //----------------------------------------------------------------------
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fuippp = 4*ui+3;
        const double v = -timefacfac*funct_(ui);
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 4*vi;
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
          estif(fvi + 2, fuippp) += v*derxy_(2, vi) ;
        }
      }

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = 4*vi;
        /* pressure term on right-hand side */
        const double v = press*rhsfac;
        eforce(fvi    ) += v*derxy_(0, vi) ;
        eforce(fvi + 1) += v*derxy_(1, vi) ;
        eforce(fvi + 2) += v*derxy_(2, vi) ;
      }

      //----------------------------------------------------------------------
      // computation of continuity term including right-hand-side contribution
      //----------------------------------------------------------------------
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvippp = 4*vi+3;
        const double v = timefacfac*funct_(vi);
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 4*ui;
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
          estif(fvippp, fui + 2) += v*derxy_(2, ui) ;
        }
      }

      const double rhsfac_vdiv = -rhsfac * vdiv_;
      for (int vi=0; vi<nen_; ++vi)
      {
        // continuity term on right-hand side
        eforce(vi*4 + 3) += rhsfac_vdiv*funct_(vi) ;
      }

      //----------------------------------------------------------------------
      // computation of body-force term on right-hand side
      //----------------------------------------------------------------------
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = 4*vi;
        const double v = fac_*funct_(vi);
        eforce(fvi    ) += v*rhsmom_(0) ;
        eforce(fvi + 1) += v*rhsmom_(1) ;
        eforce(fvi + 2) += v*rhsmom_(2) ;
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
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          const double v0 = v*derxy_(0,ui);
          const double v1 = v*derxy_(1,ui);
          const double v2 = v*derxy_(2,ui);
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
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
            estif(fvi,   fuipp) += v2*derxy_(0, vi) ;
            estif(fvip,  fui  ) += v0*derxy_(1, vi) ;
            estif(fvip,  fuip ) += v1*derxy_(1, vi) ;
            estif(fvip,  fuipp) += v2*derxy_(1, vi) ;
            estif(fvipp, fui  ) += v0*derxy_(2, vi) ;
            estif(fvipp, fuip ) += v1*derxy_(2, vi) ;
            estif(fvipp, fuipp) += v2*derxy_(2, vi) ;
          }
        }

        const double v_div = (2.0/3.0)*visceff_*rhsfac*vdiv_ ;
        const double fac_rhscon = fac_*rhscon_;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = 4*vi;
          /* viscosity term on rhs - subtraction for low-Mach-number flow */
          eforce(fvi    ) += derxy_(0, vi)*v_div ;
          eforce(fvi + 1) += derxy_(1, vi)*v_div ;
          eforce(fvi + 2) += derxy_(2, vi)*v_div ;

          /* additional rhs term of continuity equation */
          eforce(fvi + 3) += fac_rhscon*funct_(vi) ;
        }
      }

      //----------------------------------------------------------------------
      //                 PRESSURE STABILISATION PART

      // conservative, ... terms are neglected (Hughes .. .)

      if (pspg_ == INPAR::FLUID::pstab_use_pspg)
      {
        // The computation of the inertia and the convective term is split
        if (is_stationary_ == false)
        {
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            double v = tau_Mp*densam_*funct_(ui);

            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvippp = 4*vi+3;

              /* pressure stabilisation: inertia */
              /*
                /                    \
                |                    |
                |  rho*Du , nabla q  |
                |                    |
                \                   /
              */

              estif(fvippp, fui)   += v*derxy_(0, vi) ;
              estif(fvippp, fuip)  += v*derxy_(1, vi) ;
              estif(fvippp, fuipp) += v*derxy_(2, vi) ;
            }
          }
        }  // end PSPG stabilization inertia

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          double v =  timetauMp*densaf_*conv_c_(ui);

          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvippp = 4*vi+3;

            /* pressure stabilisation: convection, convective part */
            /*

            /                                     \
            |  /       n+1        \               |
            | |   rho*u   o nabla | Du , nabla q  |
            |  \      (i)         /               |
            \                                    /

            */

            estif(fvippp, fui)   += v*derxy_(0, vi) ;
            estif(fvippp, fuip)  += v*derxy_(1, vi) ;
            estif(fvippp, fuipp) += v*derxy_(2, vi) ;
          }
        }

        if (is_higher_order_ele_)
        {
          const double v = -2.0*visceff_*timetauMp;
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            for (int vi=0; vi<nen_; ++vi)
            {

              const int fvippp = 4*vi+3;
              /* pressure stabilisation: viscosity (-L_visc_u) */
              /*
                /                              \
                |               /  \             |
                |  nabla o eps | Du | , nabla q  |
                |               \  /             |
                \                              /
              */
              estif(fvippp, fui)   += v*(derxy_(0, vi)*viscs2_(0, ui)
                                         +
                                         derxy_(1, vi)*viscs2_(1, ui)
                                         +
                                         derxy_(2, vi)*viscs2_(2, ui)) ;
              estif(fvippp, fuip)  += v*(derxy_(0, vi)*viscs2_(1, ui)
                                         +
                                         derxy_(1, vi)*viscs2_(4, ui)
                                         +
                                         derxy_(2, vi)*viscs2_(5, ui)) ;
              estif(fvippp, fuipp) += v*(derxy_(0, vi)*viscs2_(2, ui)
                                         +
                                         derxy_(1, vi)*viscs2_(5, ui)
                                         +
                                         derxy_(2, vi)*viscs2_(8, ui)) ;
            }
          }
        }

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fuippp = 4*ui+3;
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
            estif(vi*4 + 3, fuippp) += timetauMp*(derxy_(0, ui)*derxy_(0, vi)
                                                  +
                                                  derxy_(1, ui)*derxy_(1, vi)
                                                  +
                                                  derxy_(2, ui)*derxy_(2, vi)) ;

          } // vi
        } // ui

        if (is_newton_)
        {
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            const double v = timetauMp*densaf_*funct_(ui);
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvippp = 4*vi + 3;
              /*  pressure stabilisation: convection, reactive part

              /                                     \
              |  /                 \  n+1           |
              | |   rho*Du o nabla | u     , grad q |
              |  \                /   (i)           |
              \                                     /

              */
              estif(fvippp, fui)   += v*(derxy_(0, vi)*vderxy_(0, 0)
                                         +
                                         derxy_(1, vi)*vderxy_(1, 0)
                                         +
                                         derxy_(2, vi)*vderxy_(2, 0)) ;
              estif(fvippp, fuip)  += v*(derxy_(0, vi)*vderxy_(0, 1)
                                         +
                                         derxy_(1, vi)*vderxy_(1, 1)
                                         +
                                         derxy_(2, vi)*vderxy_(2, 1)) ;
              estif(fvippp, fuipp) += v*(derxy_(0, vi)*vderxy_(0, 2)
                                         +
                                         derxy_(1, vi)*vderxy_(1, 2)
                                         +
                                         derxy_(2, vi)*vderxy_(2, 2)) ;

            } // vi
          } // ui
        } // if newton

        for (int vi=0; vi<nen_; ++vi)
        {
          // pressure stabilisation
          eforce(vi*4 + 3) -= tau_Mp*(momres_old_(0)*derxy_(0, vi)
                                      +
                                      momres_old_(1)*derxy_(1, vi)
                                      +
                                      momres_old_(2)*derxy_(2, vi)) ;
        }
      }

      //----------------------------------------------------------------------
      //                     SUPG STABILISATION PART

      if(supg_ == INPAR::FLUID::convective_stab_supg)
      {
        // The computation of the inertia and the convective term is split
        if(is_stationary_ == false)
        {
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            const double v = densaf_*tau_M*densam_*funct_(ui);
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = 4*vi;
              const int fvip  = fvi+1;
              const int fvipp = fvi+2;
              /* supg stabilisation: inertia  */
              /*
                /                                 \
                |           /      n+1        \    |
                |  rho*Du , | rho*u   o nabla | v  |
                |           \      (i)        /    |
                \                                 /
              */

              const double v2 = v*(conv_c_(vi)+sgconv_c_(vi));
              estif(fvi,   fui)   += v2;
              estif(fvip,  fuip)  += v2;
              estif(fvipp, fuipp) += v2;
            }
          }
        }

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          const double v = densaf_*timetauM*densaf_*conv_c_(ui);
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;

            /* supg stabilisation: convective part ( L_conv_u) */

            /*

            /                                                      \
            |    /       n+1         \      /      n+1        \     |
            |   |   rho*u    o nabla | Du , | rho*u    o nabla | v  |
            |    \       (i)         /      \      (i)        /     |
            \                                                       /

            */

            const double v2 = v*(conv_c_(vi)+sgconv_c_(vi));
            estif(fvi,   fui)   += v2;
            estif(fvip,  fuip)  += v2;
            estif(fvipp, fuipp) += v2;
          }
        }

        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;
          const double v = timetauM*densaf_*(conv_c_(vi)+sgconv_c_(vi));
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fuippp = 4*ui + 3;
            /* supg stabilisation: pressure part  ( L_pres_p) */
            /*
              /                                      \
              |              /       n+1       \     |
              |  nabla Dp , |   rho*u   o nabla | v  |
              |              \       (i)       /     |
              \                                     /
            */
            estif(fvi,   fuippp) += v*derxy_(0, ui) ;
            estif(fvip,  fuippp) += v*derxy_(1, ui) ;
            estif(fvipp, fuippp) += v*derxy_(2, ui) ;
          }
        }

        if (is_higher_order_ele_)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            const double v = -2.0*visceff_*timetauM*densaf_*(conv_c_(vi)+sgconv_c_(vi));
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              /* supg stabilisation: viscous part  (-L_visc_u) */
              /*
                /                                                \
                |               /  \    /       n+1        \     |
                |  nabla o eps | Du |, |   rho*u    o nabla | v  |
                |               \  /    \       (i)        /     |
                \                                                /
              */
              estif(fvi, fui)     += v*viscs2_(0, ui) ;
              estif(fvip, fui)    += v*viscs2_(1, ui) ;
              estif(fvipp, fui)   += v*viscs2_(2, ui) ;

              estif(fvi, fuip)    += v*viscs2_(1, ui) ;
              estif(fvip, fuip)   += v*viscs2_(4, ui) ;
              estif(fvipp, fuip)  += v*viscs2_(5, ui) ;

              estif(fvi, fuipp)   += v*viscs2_(2, ui) ;
              estif(fvip, fuipp)  += v*viscs2_(5, ui) ;
              estif(fvipp, fuipp) += v*viscs2_(8, ui) ;
            }
          }
        }

        if (is_newton_)
        {
          if (is_conservative_ == false)
          {
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              const double v = tau_M*densam_*densaf_*funct_(ui);
              const double v0 = v*velint_(0);
              const double v1 = v*velint_(1);
              const double v2 = v*velint_(2);
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi   = 4*vi;
                const int fvip  = fvi+1;
                const int fvipp = fvi+2;
                /* supg stabilisation: inertia, linearisation of testfunction  */
                /*
                  /                                         \
                  |         n+1      /                \     |
                  |    rho*u      ,  | rho*Du o nabla | v   |
                  |         (i)      \                /     |
                  \                                         /

                */
                estif(fvi,  fui)     += v0*derxy_(0, vi) ;
                estif(fvip,  fui)    += v1*derxy_(0, vi) ;
                estif(fvipp, fui)    += v2*derxy_(0, vi) ;

                estif(fvi,   fuip)   += v0*derxy_(1, vi) ;
                estif(fvip,  fuip)   += v1*derxy_(1, vi) ;
                estif(fvipp, fuip)   += v2*derxy_(1, vi) ;

                estif(fvi,   fuipp)  += v0*derxy_(2, vi) ;
                estif(fvip,  fuipp)  += v1*derxy_(2, vi) ;
                estif(fvipp, fuipp)  += v2*derxy_(2, vi) ;
              }
            }
          }

          {
            const double v0 = convvelint_(0)*vderxy_(0, 0) + convvelint_(1)*vderxy_(0, 1) + convvelint_(2)*vderxy_(0, 2);
            const double v1 = convvelint_(0)*vderxy_(1, 0) + convvelint_(1)*vderxy_(1, 1) + convvelint_(2)*vderxy_(1, 2);
            const double v2 = convvelint_(0)*vderxy_(2, 0) + convvelint_(1)*vderxy_(2, 1) + convvelint_(2)*vderxy_(2, 2);

            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              const double v = timetauM*densaf_*densaf_*funct_(ui);
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi   = 4*vi;
                const int fvip  = fvi+1;
                const int fvipp = fvi+2;
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
                estif(fvi, fui)     += (conv_c_(vi,0)*vderxy_(0, 0) + v0*derxy_(0, vi))*v;
                estif(fvip, fui)    += (conv_c_(vi,0)*vderxy_(1, 0) + v1*derxy_(0, vi))*v;
                estif(fvipp, fui)   += (conv_c_(vi,0)*vderxy_(2, 0) + v2*derxy_(0, vi))*v;

                estif(fvi, fuip)    += (conv_c_(vi,0)*vderxy_(0, 1) + v0*derxy_(1, vi))*v;
                estif(fvip, fuip)   += (conv_c_(vi,0)*vderxy_(1, 1) + v1*derxy_(1, vi))*v;
                estif(fvipp, fuip)  += (conv_c_(vi,0)*vderxy_(2, 1) + v2*derxy_(1, vi))*v;

                estif(fvi, fuipp)   += (conv_c_(vi,0)*vderxy_(0, 2) + v0*derxy_(2, vi))*v;
                estif(fvip, fuipp)  += (conv_c_(vi,0)*vderxy_(1, 2) + v1*derxy_(2, vi))*v;
                estif(fvipp, fuipp) += (conv_c_(vi,0)*vderxy_(2, 2) + v2*derxy_(2, vi))*v;
              }
            }
          }

          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            const double v = timetauM*densaf_*funct_(ui);
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = 4*vi;
              const int fvip  = fvi+1;
              const int fvipp = fvi+2;
              /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
              /*
                /                                       \
                |         n+1    /                \     |
                |  nabla p    , |   rho*Du o nabla | v  |
                |         (i)    \                /     |
                \                                      /
              */
              estif(fvi, fui)     += v*gradp_(0,0)*derxy_(0, vi) ;
              estif(fvip, fui)    += v*gradp_(1,0)*derxy_(0, vi) ;
              estif(fvipp, fui)   += v*gradp_(2,0)*derxy_(0, vi) ;

              estif(fvi, fuip)    += v*gradp_(0,0)*derxy_(1, vi) ;
              estif(fvip, fuip)   += v*gradp_(1,0)*derxy_(1, vi) ;
              estif(fvipp, fuip)  += v*gradp_(2,0)*derxy_(1, vi) ;

              estif(fvi, fuipp)   += v*gradp_(0,0)*derxy_(2, vi) ;
              estif(fvip, fuipp)  += v*gradp_(1,0)*derxy_(2, vi) ;
              estif(fvipp, fuipp) += v*gradp_(2,0)*derxy_(2, vi) ;

            }
          }

          if (is_higher_order_ele_)
          {
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = 4*ui;
              const int fuip  = fui+1;
              const int fuipp = fui+2;
              const double v = -2.0*visceff_*timetauM*densaf_*funct_(ui);
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi   = 4*vi;
                const int fvip  = fvi+1;
                const int fvipp = fvi+2;
                const double v0 = v*visc_old_(0);
                const double v1 = v*visc_old_(1);
                const double v2 = v*visc_old_(2);

                /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
                /*
                  /                                                 \
                  |               / n+1 \    /                \     |
                  |  nabla o eps | u     |, |  rho*Du o nabla | v   |
                  |               \ (i) /    \                /     |
                  \                                                 /
                */
                estif(fvi, fui)     += v0*derxy_(0, vi) ;
                estif(fvip, fui)    += v1*derxy_(0, vi) ;
                estif(fvipp, fui)   += v2*derxy_(0, vi) ;

                estif(fvi, fuip)    += v0*derxy_(1, vi) ;
                estif(fvip, fuip)   += v1*derxy_(1, vi) ;
                estif(fvipp, fuip)  += v2*derxy_(1, vi) ;

                estif(fvi, fuipp)   += v0*derxy_(2, vi) ;
                estif(fvip, fuipp)  += v1*derxy_(2, vi) ;
                estif(fvipp, fuipp) += v2*derxy_(2, vi) ;
              }
            }
          }

          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = 4*ui;
            const int fuip  = fui+1;
            const int fuipp = fui+2;
            const double v = -tau_M*densaf_*funct_(ui);
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = 4*vi;
              const int fvip  = fvi+1;
              const int fvipp = fvi+2;
              const double v0 = v*rhsmom_(0);
              const double v1 = v*rhsmom_(1);
              const double v2 = v*rhsmom_(2);

              /* supg stabilisation: bodyforce part, linearisation of test function */

              /*
                /                                       \
                |                 /                \     |
                |  rho*rhsint   , |  rho*Du o nabla | v  |
                |                 \                /     |
                \                                        /

              */
              estif(fvi    , fui)   += (v0*derxy_(0, vi)) ;
              estif(fvip   , fui)   += (v1*derxy_(0, vi)) ;
              estif(fvipp  , fui)   += (v2*derxy_(0, vi)) ;

              estif(fvi    , fuip)  += (v0*derxy_(1, vi)) ;
              estif(fvip   , fuip)  += (v1*derxy_(1, vi)) ;
              estif(fvipp  , fuip)  += (v2*derxy_(1, vi)) ;

              estif(fvi    , fuipp) += (v0*derxy_(2, vi)) ;
              estif(fvip   , fuipp) += (v1*derxy_(2, vi)) ;
              estif(fvipp  , fuipp) += (v2*derxy_(2, vi)) ;

            } // vi
          } // ui
        } // end if(newton)

        // NOTE: Here we have a difference to the previous version of this
        // element!  Before we did not care for the mesh velocity in this
        // term. This seems unreasonable and wrong.
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = 4*vi;
          // supg stabilisation
          const double v = -tau_M*densaf_*(conv_c_(vi)+sgconv_c_(vi));
          eforce(fvi)     += v*momres_old_(0) ;
          eforce(fvi + 1) += v*momres_old_(1) ;
          eforce(fvi + 2) += v*momres_old_(2) ;
        }
      }

      //----------------------------------------------------------------------
      //                       STABILISATION, VISCOUS PART

      if (is_higher_order_ele_)
      {
        if(vstab_ != INPAR::FLUID::viscous_stab_none)
        {
          const double two_visc_tauMp = vstabfac*2.0*visc_*tau_Mp;
          // viscous stabilization either on left hand side or on right hand side
          if (vstab_ == INPAR::FLUID::viscous_stab_gls || vstab_ == INPAR::FLUID::viscous_stab_usfem)
          {
            const double two_visc_timetauMp   = vstabfac*2.0*visc_*timetauMp;
            const double four_visc2_timetauMp = vstabfac*4.0*visceff_*visc_*timetauMp;

            // viscous stabilization on left hand side
            if (is_stationary_ == false)
            {
              for (int ui=0; ui<nen_; ++ui)
              {
                const double v = two_visc_tauMp*densam_*funct_(ui);
                for (int vi=0; vi<nen_; ++vi)
                {
                  /* viscous stabilisation, inertia part */
                  /*
                    /                          \
                    |                          |
                +/- |    rho*Du , div eps (v)  |
                    |                          |
                    \                          /
                  */

                  estif(vi*4,     ui*4    ) += v*viscs2_(0, vi) ;
                  estif(vi*4 + 1, ui*4    ) += v*viscs2_(1, vi) ;
                  estif(vi*4 + 2, ui*4    ) += v*viscs2_(2, vi) ;

                  estif(vi*4,     ui*4 + 1) += v*viscs2_(1, vi) ;
                  estif(vi*4 + 1, ui*4 + 1) += v*viscs2_(4, vi) ;
                  estif(vi*4 + 2, ui*4 + 1) += v*viscs2_(5, vi) ;

                  estif(vi*4,     ui*4 + 2) += v*viscs2_(2, vi) ;
                  estif(vi*4 + 1, ui*4 + 2) += v*viscs2_(5, vi) ;
                  estif(vi*4 + 2, ui*4 + 2) += v*viscs2_(8, vi) ;
                }
              }
            }

            // viscous stabilization on left hand side
            for (int ui=0; ui<nen_; ++ui)
            {
              const double v = two_visc_timetauMp*densaf_*conv_c_(ui);
              for (int vi=0; vi<nen_; ++vi)
              {
                /* viscous stabilisation, inertia part */

                /* viscous stabilisation, convective part */
                /*
                  /                                        \
                  |  /       n+1       \                   |
              +/- | |   rho*u   o nabla | Du , div eps (v) |
                  |  \       (i)       /                   |
                  \                                        /
                */
                estif(vi*4,     ui*4    ) += v*viscs2_(0, vi) ;
                estif(vi*4 + 1, ui*4    ) += v*viscs2_(1, vi) ;
                estif(vi*4 + 2, ui*4    ) += v*viscs2_(2, vi) ;

                estif(vi*4,     ui*4 + 1) += v*viscs2_(1, vi) ;
                estif(vi*4 + 1, ui*4 + 1) += v*viscs2_(4, vi) ;
                estif(vi*4 + 2, ui*4 + 1) += v*viscs2_(5, vi) ;

                estif(vi*4,     ui*4 + 2) += v*viscs2_(2, vi) ;
                estif(vi*4 + 1, ui*4 + 2) += v*viscs2_(5, vi) ;
                estif(vi*4 + 2, ui*4 + 2) += v*viscs2_(8, vi) ;
              }
            }

            for (int ui=0; ui<nen_; ++ui)
            {
              for (int vi=0; vi<nen_; ++vi)
              {

                /* viscous stabilisation, pressure part ( L_pres_p) */
                /*
                  /                          \
                  |                          |
             +/-  |  nabla Dp , div eps (v)  |
                  |                          |
                  \                          /
                */
                estif(vi*4,     ui*4 + 3) += two_visc_timetauMp*(derxy_(0, ui)*viscs2_(0, vi)
                                                                +
                                                                derxy_(1, ui)*viscs2_(1, vi)
                                                                +
                                                                derxy_(2, ui)*viscs2_(2, vi)) ;
                estif(vi*4 + 1, ui*4 + 3) += two_visc_timetauMp*(derxy_(0, ui)*viscs2_(1, vi)
                                                                +
                                                                derxy_(1, ui)*viscs2_(4, vi)
                                                                +
                                                                derxy_(2, ui)*viscs2_(5, vi)) ;
                estif(vi*4 + 2, ui*4 + 3) += two_visc_timetauMp*(derxy_(0, ui)*viscs2_(2, vi)
                                                                +
                                                                derxy_(1, ui)*viscs2_(5, vi)
                                                                +
                                                                derxy_(2, ui)*viscs2_(8, vi)) ;

              }
            }

            for (int ui=0; ui<nen_; ++ui)
            {
              for (int vi=0; vi<nen_; ++vi)
              {
                /* viscous stabilisation, viscous part (-L_visc_u) */
                /*
                  /                                 \
                  |               /  \                |
             -/+  |  nabla o eps | Du | , div eps (v) |
                  |               \  /                |
                  \                                 /
                */
                estif(vi*4,     ui*4    ) -= four_visc2_timetauMp*(viscs2_(0,ui)*viscs2_(0,vi)+viscs2_(1,ui)*viscs2_(1,vi)+viscs2_(2,ui)*viscs2_(2,vi)) ;
                estif(vi*4 + 1, ui*4    ) -= four_visc2_timetauMp*(viscs2_(0,ui)*viscs2_(1,vi)+viscs2_(1,ui)*viscs2_(4,vi)+viscs2_(2,ui)*viscs2_(5,vi)) ;
                estif(vi*4 + 2, ui*4    ) -= four_visc2_timetauMp*(viscs2_(0,ui)*viscs2_(2,vi)+viscs2_(1,ui)*viscs2_(5,vi)+viscs2_(2,ui)*viscs2_(8,vi)) ;

                estif(vi*4,     ui*4 + 1) -= four_visc2_timetauMp*(viscs2_(0,vi)*viscs2_(1,ui)+viscs2_(1,vi)*viscs2_(4,ui)+viscs2_(2,vi)*viscs2_(5,ui)) ;
                estif(vi*4 + 1, ui*4 + 1) -= four_visc2_timetauMp*(viscs2_(1,ui)*viscs2_(1,vi)+viscs2_(4,ui)*viscs2_(4,vi)+viscs2_(5,ui)*viscs2_(5,vi)) ;
                estif(vi*4 + 2, ui*4 + 1) -= four_visc2_timetauMp*(viscs2_(1,ui)*viscs2_(2,vi)+viscs2_(4,ui)*viscs2_(5,vi)+viscs2_(5,ui)*viscs2_(8,vi)) ;

                estif(vi*4,     ui*4 + 2) -= four_visc2_timetauMp*(viscs2_(0,vi)*viscs2_(2,ui)+viscs2_(1,vi)*viscs2_(5,ui)+viscs2_(2,vi)*viscs2_(8,ui)) ;
                estif(vi*4 + 1, ui*4 + 2) -= four_visc2_timetauMp*(viscs2_(1,vi)*viscs2_(2,ui)+viscs2_(4,vi)*viscs2_(5,ui)+viscs2_(5,vi)*viscs2_(8,ui)) ;
                estif(vi*4 + 2, ui*4 + 2) -= four_visc2_timetauMp*(viscs2_(2,ui)*viscs2_(2,vi)+viscs2_(5,ui)*viscs2_(5,vi)+viscs2_(8,ui)*viscs2_(8,vi)) ;
              } // vi
            } // ui

            if (is_newton_)
            {
              for (int ui=0; ui<nen_; ++ui)
              {
                double v = two_visc_timetauMp*densaf_*funct_(ui);
                for (int vi=0; vi<nen_; ++vi)
                {
                  /* viscous stabilisation, reactive part of convection */
                  /*
                    /                                         \
                    |  /                \   n+1               |
                +/- | |   rho*Du o nabla | u    , div eps (v) |
                    |  \                /   (i)               |
                    \                                         /
                  */
                  estif(vi*4,     ui*4    ) += v*(viscs2_(0,vi)*vderxy_(0,0)+
                                                  viscs2_(1,vi)*vderxy_(1,0)+
                                                  viscs2_(2,vi)*vderxy_(2,0)) ;
                  estif(vi*4 + 1, ui*4    ) += v*(viscs2_(1,vi)*vderxy_(0,0)+
                                                  viscs2_(4,vi)*vderxy_(1,0)+
                                                  viscs2_(5,vi)*vderxy_(2,0)) ;
                  estif(vi*4 + 2, ui*4    ) += v*(viscs2_(2,vi)*vderxy_(0,0)+
                                                  viscs2_(5,vi)*vderxy_(1,0)+
                                                  viscs2_(8,vi)*vderxy_(2,0)) ;

                  estif(vi*4,     ui*4 + 1) += v*(viscs2_(0,vi)*vderxy_(0,1)+
                                                  viscs2_(1,vi)*vderxy_(1,1)+
                                                  viscs2_(2,vi)*vderxy_(2,1)) ;
                  estif(vi*4 + 1, ui*4 + 1) += v*(viscs2_(1,vi)*vderxy_(0,1)+
                                                  viscs2_(4,vi)*vderxy_(1,1)+
                                                  viscs2_(5,vi)*vderxy_(2,1)) ;
                  estif(vi*4 + 2, ui*4 + 1) += v*(viscs2_(2,vi)*vderxy_(0,1)+
                                                  viscs2_(5,vi)*vderxy_(1,1)+
                                                  viscs2_(8,vi)*vderxy_(2,1)) ;

                  estif(vi*4,     ui*4 + 2) += v*(viscs2_(0,vi)*vderxy_(0,2)+
                                                  viscs2_(1,vi)*vderxy_(1,2)+
                                                  viscs2_(2,vi)*vderxy_(2,2)) ;
                  estif(vi*4 + 1, ui*4 + 2) += v*(viscs2_(1,vi)*vderxy_(0,2)+
                                                  viscs2_(4,vi)*vderxy_(1,2)+
                                                  viscs2_(5,vi)*vderxy_(2,2)) ;
                  estif(vi*4 + 2, ui*4 + 2) += v*(viscs2_(2,vi)*vderxy_(0,2)+
                                                  viscs2_(5,vi)*vderxy_(1,2)+
                                                  viscs2_(8,vi)*vderxy_(2,2)) ;
                } // vi
              } // ui
            } // if newton
          } // end if viscous stabilization on left hand side

          for (int vi=0; vi<nen_; ++vi)
          {

            /* viscous stabilisation */
            eforce(vi*4    ) -= two_visc_tauMp*(momres_old_(0)*viscs2_(0, vi)+momres_old_(1)*viscs2_(1, vi)+momres_old_(2)*viscs2_(2, vi)) ;
            eforce(vi*4 + 1) -= two_visc_tauMp*(momres_old_(0)*viscs2_(1, vi)+momres_old_(1)*viscs2_(4, vi)+momres_old_(2)*viscs2_(5, vi)) ;
            eforce(vi*4 + 2) -= two_visc_tauMp*(momres_old_(0)*viscs2_(2, vi)+momres_old_(1)*viscs2_(5, vi)+momres_old_(2)*viscs2_(8, vi)) ;
          }
        }
      }

      //----------------------------------------------------------------------
      //                     STABILISATION, CONTINUITY PART

      if (cstab_ == INPAR::FLUID::continuity_stab_yes)
      {
        const double timetauC = timefac*tau_C;
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = 4*ui;
          const int fuip  = fui+1;
          const int fuipp = fui+2;
          double v0 = timetauC*derxy_(0, ui);
          double v1 = timetauC*derxy_(1, ui);
          double v2 = timetauC*derxy_(2, ui);
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = 4*vi;
            const int fvip  = fvi+1;
            const int fvipp = fvi+2;
            /* continuity stabilisation on left hand side */
            /*
              /                         \
              |                          |
              | nabla o Du  , nabla o v  |
              |                          |
              \                         /
            */
            estif(fvi,  fui  ) += v0*derxy_(0, vi) ;
            estif(fvip, fui  ) += v0*derxy_(1, vi) ;
            estif(fvipp,fui  ) += v0*derxy_(2, vi) ;

            estif(fvi,  fuip ) += v1*derxy_(0, vi) ;
            estif(fvip, fuip ) += v1*derxy_(1, vi) ;
            estif(fvipp,fuip ) += v1*derxy_(2, vi) ;

            estif(fvi,  fuipp) += v2*derxy_(0, vi) ;
            estif(fvip, fuipp) += v2*derxy_(1, vi) ;
            estif(fvipp,fuipp) += v2*derxy_(2, vi) ;
          }
        }

        const double tauC_conres = tau_C*conres_old_;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = 4*vi;
          const int fvip  = fvi+1;
          const int fvipp = fvi+2;
          /* continuity stabilisation on right hand side */
          eforce(fvi  ) -= tauC_conres*derxy_(0, vi) ;
          eforce(fvip ) -= tauC_conres*derxy_(1, vi) ;
          eforce(fvipp) -= tauC_conres*derxy_(2, vi) ;
        }
      }

      //TODO: Cross-stress and Reynolds-Stress stabilization are missing

      //----------------------------------------------------------------------
      //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)

      //TODO: Difference stationary <-> instationary

      if(fssgv_ != Fluid3::no_fssgv)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = 4*vi;
          /* fine-scale subgrid-viscosity term on right hand side */
          /*
                              /                          \
                             |       /    \         / \   |
             - mu_art(fsu) * |  eps | Dfsu | , eps | v |  |
                             |       \    /         \ /   |
                              \                          /
          */
          eforce(fvi    ) -= fssgviscfac*(2.0*derxy_(0, vi)*fsvderxy_(0, 0)
                                         +    derxy_(1, vi)*fsvderxy_(0, 1)
                                         +    derxy_(1, vi)*fsvderxy_(1, 0)
                                         +    derxy_(2, vi)*fsvderxy_(0, 2)
                                         +    derxy_(2, vi)*fsvderxy_(2, 0)) ;
          eforce(fvi + 1) -= fssgviscfac*(    derxy_(0, vi)*fsvderxy_(0, 1)
                                         +    derxy_(0, vi)*fsvderxy_(1, 0)
                                         +2.0*derxy_(1, vi)*fsvderxy_(1, 1)
                                         +    derxy_(2, vi)*fsvderxy_(1, 2)
                                         +    derxy_(2, vi)*fsvderxy_(2, 1)) ;
          eforce(fvi + 2) -= fssgviscfac*(    derxy_(0, vi)*fsvderxy_(0, 2)
                                         +    derxy_(0, vi)*fsvderxy_(2, 0)
                                         +    derxy_(1, vi)*fsvderxy_(1, 2)
                                         +    derxy_(1, vi)*fsvderxy_(2, 1)
                                         +2.0*derxy_(2, vi)*fsvderxy_(2, 2)) ;
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
        double v = fac_*funct_(vi,0);
        for (int ui=0; ui<nen_; ++ui)
        {
          emesh(vi*4    , ui*4    ) += v*(velint_(0)-rhsmom_(0))*derxy_(0, ui);
          emesh(vi*4    , ui*4 + 1) += v*(velint_(0)-rhsmom_(0))*derxy_(1, ui);
          emesh(vi*4    , ui*4 + 2) += v*(velint_(0)-rhsmom_(0))*derxy_(2, ui);

          emesh(vi*4 + 1, ui*4    ) += v*(velint_(1)-rhsmom_(1))*derxy_(0, ui);
          emesh(vi*4 + 1, ui*4 + 1) += v*(velint_(1)-rhsmom_(1))*derxy_(1, ui);
          emesh(vi*4 + 1, ui*4 + 2) += v*(velint_(1)-rhsmom_(1))*derxy_(2, ui);

          emesh(vi*4 + 2, ui*4    ) += v*(velint_(2)-rhsmom_(2))*derxy_(0, ui);
          emesh(vi*4 + 2, ui*4 + 1) += v*(velint_(2)-rhsmom_(2))*derxy_(1, ui);
          emesh(vi*4 + 2, ui*4 + 2) += v*(velint_(2)-rhsmom_(2))*derxy_(2, ui);
        }
      }

      //vderiv_  = sum(evelaf(i,k) * deriv_(j,k), k);
      vderiv_.MultiplyNT(evelaf,deriv_);

#define derxjm_(r,c,d,i) derxjm_ ## r ## c ## d (i)

#define derxjm_001(ui) (deriv_(2, ui)*xjm_(1, 2) - deriv_(1, ui)*xjm_(2, 2))
#define derxjm_002(ui) (deriv_(1, ui)*xjm_(2, 1) - deriv_(2, ui)*xjm_(1, 1))

#define derxjm_100(ui) (deriv_(1, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(1, 2))
#define derxjm_102(ui) (deriv_(2, ui)*xjm_(1, 0) - deriv_(1, ui)*xjm_(2, 0))

#define derxjm_200(ui) (deriv_(2, ui)*xjm_(1, 1) - deriv_(1, ui)*xjm_(2, 1))
#define derxjm_201(ui) (deriv_(1, ui)*xjm_(2, 0) - deriv_(2, ui)*xjm_(1, 0))

#define derxjm_011(ui) (deriv_(0, ui)*xjm_(2, 2) - deriv_(2, ui)*xjm_(0, 2))
#define derxjm_012(ui) (deriv_(2, ui)*xjm_(0, 1) - deriv_(0, ui)*xjm_(2, 1))

#define derxjm_110(ui) (deriv_(2, ui)*xjm_(0, 2) - deriv_(0, ui)*xjm_(2, 2))
#define derxjm_112(ui) (deriv_(0, ui)*xjm_(2, 0) - deriv_(2, ui)*xjm_(0, 0))

#define derxjm_210(ui) (deriv_(0, ui)*xjm_(2, 1) - deriv_(2, ui)*xjm_(0, 1))
#define derxjm_211(ui) (deriv_(2, ui)*xjm_(0, 0) - deriv_(0, ui)*xjm_(2, 0))

#define derxjm_021(ui) (deriv_(1, ui)*xjm_(0, 2) - deriv_(0, ui)*xjm_(1, 2))
#define derxjm_022(ui) (deriv_(0, ui)*xjm_(1, 1) - deriv_(1, ui)*xjm_(0, 1))

#define derxjm_120(ui) (deriv_(0, ui)*xjm_(1, 2) - deriv_(1, ui)*xjm_(0, 2))
#define derxjm_122(ui) (deriv_(1, ui)*xjm_(0, 0) - deriv_(0, ui)*xjm_(1, 0))

#define derxjm_220(ui) (deriv_(1, ui)*xjm_(0, 1) - deriv_(0, ui)*xjm_(1, 1))
#define derxjm_221(ui) (deriv_(0, ui)*xjm_(1, 0) - deriv_(1, ui)*xjm_(0, 0))

      for (int ui=0; ui<nen_; ++ui)
      {
        double v00 = + convvelint_(1)*(vderiv_(0, 0)*derxjm_(0,0,1,ui) + vderiv_(0, 1)*derxjm_(0,1,1,ui) + vderiv_(0, 2)*derxjm_(0,2,1,ui))
                     + convvelint_(2)*(vderiv_(0, 0)*derxjm_(0,0,2,ui) + vderiv_(0, 1)*derxjm_(0,1,2,ui) + vderiv_(0, 2)*derxjm_(0,2,2,ui));
        double v01 = + convvelint_(0)*(vderiv_(0, 0)*derxjm_(1,0,0,ui) + vderiv_(0, 1)*derxjm_(1,1,0,ui) + vderiv_(0, 2)*derxjm_(1,2,0,ui))
                     + convvelint_(2)*(vderiv_(0, 0)*derxjm_(1,0,2,ui) + vderiv_(0, 1)*derxjm_(1,1,2,ui) + vderiv_(0, 2)*derxjm_(1,2,2,ui));
        double v02 = + convvelint_(0)*(vderiv_(0, 0)*derxjm_(2,0,0,ui) + vderiv_(0, 1)*derxjm_(2,1,0,ui) + vderiv_(0, 2)*derxjm_(2,2,0,ui))
                     + convvelint_(1)*(vderiv_(0, 0)*derxjm_(2,0,1,ui) + vderiv_(0, 1)*derxjm_(2,1,1,ui) + vderiv_(0, 2)*derxjm_(2,2,1,ui));
        double v10 = + convvelint_(1)*(vderiv_(1, 0)*derxjm_(0,0,1,ui) + vderiv_(1, 1)*derxjm_(0,1,1,ui) + vderiv_(1, 2)*derxjm_(0,2,1,ui))
                     + convvelint_(2)*(vderiv_(1, 0)*derxjm_(0,0,2,ui) + vderiv_(1, 1)*derxjm_(0,1,2,ui) + vderiv_(1, 2)*derxjm_(0,2,2,ui));
        double v11 = + convvelint_(0)*(vderiv_(1, 0)*derxjm_(1,0,0,ui) + vderiv_(1, 1)*derxjm_(1,1,0,ui) + vderiv_(1, 2)*derxjm_(1,2,0,ui))
                     + convvelint_(2)*(vderiv_(1, 0)*derxjm_(1,0,2,ui) + vderiv_(1, 1)*derxjm_(1,1,2,ui) + vderiv_(1, 2)*derxjm_(1,2,2,ui));
        double v12 = + convvelint_(0)*(vderiv_(1, 0)*derxjm_(2,0,0,ui) + vderiv_(1, 1)*derxjm_(2,1,0,ui) + vderiv_(1, 2)*derxjm_(2,2,0,ui))
                     + convvelint_(1)*(vderiv_(1, 0)*derxjm_(2,0,1,ui) + vderiv_(1, 1)*derxjm_(2,1,1,ui) + vderiv_(1, 2)*derxjm_(2,2,1,ui));
        double v20 = + convvelint_(1)*(vderiv_(2, 0)*derxjm_(0,0,1,ui) + vderiv_(2, 1)*derxjm_(0,1,1,ui) + vderiv_(2, 2)*derxjm_(0,2,1,ui))
                     + convvelint_(2)*(vderiv_(2, 0)*derxjm_(0,0,2,ui) + vderiv_(2, 1)*derxjm_(0,1,2,ui) + vderiv_(2, 2)*derxjm_(0,2,2,ui));
        double v21 = + convvelint_(0)*(vderiv_(2, 0)*derxjm_(1,0,0,ui) + vderiv_(2, 1)*derxjm_(1,1,0,ui) + vderiv_(2, 2)*derxjm_(1,2,0,ui))
                     + convvelint_(2)*(vderiv_(2, 0)*derxjm_(1,0,2,ui) + vderiv_(2, 1)*derxjm_(1,1,2,ui) + vderiv_(2, 2)*derxjm_(1,2,2,ui));
        double v22 = + convvelint_(0)*(vderiv_(2, 0)*derxjm_(2,0,0,ui) + vderiv_(2, 1)*derxjm_(2,1,0,ui) + vderiv_(2, 2)*derxjm_(2,2,0,ui))
                     + convvelint_(1)*(vderiv_(2, 0)*derxjm_(2,0,1,ui) + vderiv_(2, 1)*derxjm_(2,1,1,ui) + vderiv_(2, 2)*derxjm_(2,2,1,ui));

        for (int vi=0; vi<nen_; ++vi)
        {
          double v = timefacfac/det_*funct_(vi);

          emesh(vi*4 + 0, ui*4 + 0) += v*v00;
          emesh(vi*4 + 0, ui*4 + 1) += v*v01;
          emesh(vi*4 + 0, ui*4 + 2) += v*v02;

          emesh(vi*4 + 1, ui*4 + 0) += v*v10;
          emesh(vi*4 + 1, ui*4 + 1) += v*v11;
          emesh(vi*4 + 1, ui*4 + 2) += v*v12;

          emesh(vi*4 + 2, ui*4 + 0) += v*v20;
          emesh(vi*4 + 2, ui*4 + 1) += v*v21;
          emesh(vi*4 + 2, ui*4 + 2) += v*v22;
        }
      }

      // viscosity

#define xji_00 xji_(0,0)
#define xji_01 xji_(0,1)
#define xji_02 xji_(0,2)
#define xji_10 xji_(1,0)
#define xji_11 xji_(1,1)
#define xji_12 xji_(1,2)
#define xji_20 xji_(2,0)
#define xji_21 xji_(2,1)
#define xji_22 xji_(2,2)

#define xjm(i,j) xjm_(i,j)

      // part 1: derivative of 1/det

      double v = visceff_*timefac*fac_;
      for (int ui=0; ui<nen_; ++ui)
      {
        double derinvJ0 = -v*(deriv_(0,ui)*xji_00 + deriv_(1,ui)*xji_01 + deriv_(2,ui)*xji_02);
        double derinvJ1 = -v*(deriv_(0,ui)*xji_10 + deriv_(1,ui)*xji_11 + deriv_(2,ui)*xji_12);
        double derinvJ2 = -v*(deriv_(0,ui)*xji_20 + deriv_(1,ui)*xji_21 + deriv_(2,ui)*xji_22);
        for (int vi=0; vi<nen_; ++vi)
        {
          double visres0 =   2.0*derxy_(0, vi)* vderxy_(0, 0)
                             +     derxy_(1, vi)*(vderxy_(0, 1) + vderxy_(1, 0))
                             +     derxy_(2, vi)*(vderxy_(0, 2) + vderxy_(2, 0)) ;
          double visres1 =         derxy_(0, vi)*(vderxy_(0, 1) + vderxy_(1, 0))
                             + 2.0*derxy_(1, vi)* vderxy_(1, 1)
                             +     derxy_(2, vi)*(vderxy_(1, 2) + vderxy_(2, 1)) ;
          double visres2 =         derxy_(0, vi)*(vderxy_(0, 2) + vderxy_(2, 0))
                             +     derxy_(1, vi)*(vderxy_(1, 2) + vderxy_(2, 1))
                             + 2.0*derxy_(2, vi)* vderxy_(2, 2) ;
          emesh(vi*4 + 0, ui*4 + 0) += derinvJ0*visres0;
          emesh(vi*4 + 1, ui*4 + 0) += derinvJ0*visres1;
          emesh(vi*4 + 2, ui*4 + 0) += derinvJ0*visres2;

          emesh(vi*4 + 0, ui*4 + 1) += derinvJ1*visres0;
          emesh(vi*4 + 1, ui*4 + 1) += derinvJ1*visres1;
          emesh(vi*4 + 2, ui*4 + 1) += derinvJ1*visres2;

          emesh(vi*4 + 0, ui*4 + 2) += derinvJ2*visres0;
          emesh(vi*4 + 1, ui*4 + 2) += derinvJ2*visres1;
          emesh(vi*4 + 2, ui*4 + 2) += derinvJ2*visres2;
        }
      }

      // part 2: derivative of viscosity residual

      v = timefacfac*visceff_/det_;
      for (int ui=0; ui<nen_; ++ui)
      {
        double v0 = - vderiv_(0,0)*(xji_10*derxjm_100(ui) + xji_10*derxjm_100(ui) + xji_20*derxjm_200(ui) + xji_20*derxjm_200(ui))
                    - vderiv_(0,1)*(xji_11*derxjm_100(ui) + xji_10*derxjm_110(ui) + xji_21*derxjm_200(ui) + xji_20*derxjm_210(ui))
                    - vderiv_(0,2)*(xji_12*derxjm_100(ui) + xji_10*derxjm_120(ui) + xji_22*derxjm_200(ui) + xji_20*derxjm_220(ui))
                    - vderiv_(1,0)*(derxjm_100(ui)*xji_00)
                    - vderiv_(1,1)*(derxjm_100(ui)*xji_01)
                    - vderiv_(1,2)*(derxjm_100(ui)*xji_02)
                    - vderiv_(2,0)*(derxjm_200(ui)*xji_00)
                    - vderiv_(2,1)*(derxjm_200(ui)*xji_01)
                    - vderiv_(2,2)*(derxjm_200(ui)*xji_02);
        double v1 = - vderiv_(0,0)*(xji_10*derxjm_110(ui) + xji_11*derxjm_100(ui) + xji_20*derxjm_210(ui) + xji_21*derxjm_200(ui))
                    - vderiv_(0,1)*(xji_11*derxjm_110(ui) + xji_11*derxjm_110(ui) + xji_21*derxjm_210(ui) + xji_21*derxjm_210(ui))
                    - vderiv_(0,2)*(xji_12*derxjm_110(ui) + xji_11*derxjm_120(ui) + xji_22*derxjm_210(ui) + xji_21*derxjm_220(ui))
                    - vderiv_(1,0)*(derxjm_110(ui)*xji_00)
                    - vderiv_(1,1)*(derxjm_110(ui)*xji_01)
                    - vderiv_(1,2)*(derxjm_110(ui)*xji_02)
                    - vderiv_(2,0)*(derxjm_210(ui)*xji_00)
                    - vderiv_(2,1)*(derxjm_210(ui)*xji_01)
                    - vderiv_(2,2)*(derxjm_210(ui)*xji_02);
        double v2 = - vderiv_(0,0)*(xji_10*derxjm_120(ui) + xji_12*derxjm_100(ui) + xji_20*derxjm_220(ui) + xji_22*derxjm_200(ui))
                    - vderiv_(0,1)*(xji_11*derxjm_120(ui) + xji_12*derxjm_110(ui) + xji_21*derxjm_220(ui) + xji_22*derxjm_210(ui))
                    - vderiv_(0,2)*(xji_12*derxjm_120(ui) + xji_12*derxjm_120(ui) + xji_22*derxjm_220(ui) + xji_22*derxjm_220(ui))
                    - vderiv_(1,0)*(derxjm_120(ui)*xji_00)
                    - vderiv_(1,1)*(derxjm_120(ui)*xji_01)
                    - vderiv_(1,2)*(derxjm_120(ui)*xji_02)
                    - vderiv_(2,0)*(derxjm_220(ui)*xji_00)
                    - vderiv_(2,1)*(derxjm_220(ui)*xji_01)
                    - vderiv_(2,2)*(derxjm_220(ui)*xji_02);

        for (int vi=0; vi<nen_; ++vi)
        {
          emesh(vi*4 + 0, ui*4 + 0) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(2*derxjm_001(ui)*xji_00 + 2*derxjm_001(ui)*xji_00 + xji_20*derxjm_201(ui) + xji_20*derxjm_201(ui))
             - vderiv_(0,1)*(2*derxjm_011(ui)*xji_00 + 2*derxjm_001(ui)*xji_01 + xji_21*derxjm_201(ui) + xji_20*derxjm_211(ui))
             - vderiv_(0,2)*(2*derxjm_021(ui)*xji_00 + 2*derxjm_001(ui)*xji_02 + xji_22*derxjm_201(ui) + xji_20*derxjm_221(ui))
             - vderiv_(1,0)*(derxjm_001(ui)*xji_10)
             - vderiv_(1,1)*(derxjm_011(ui)*xji_10)
             - vderiv_(1,2)*(derxjm_021(ui)*xji_10)
             - vderiv_(2,0)*(derxjm_201(ui)*xji_00 + derxjm_001(ui)*xji_20)
             - vderiv_(2,1)*(derxjm_201(ui)*xji_01 + derxjm_011(ui)*xji_20)
             - vderiv_(2,2)*(derxjm_201(ui)*xji_02 + derxjm_021(ui)*xji_20);
        v1 = - vderiv_(0,0)*(2*derxjm_011(ui)*xji_00 + 2*derxjm_001(ui)*xji_01 + xji_21*derxjm_201(ui) + xji_20*derxjm_211(ui))
             - vderiv_(0,1)*(2*derxjm_011(ui)*xji_01 + 2*derxjm_011(ui)*xji_01 + xji_21*derxjm_211(ui) + xji_21*derxjm_211(ui))
             - vderiv_(0,2)*(2*derxjm_011(ui)*xji_02 + 2*derxjm_021(ui)*xji_01 + xji_21*derxjm_221(ui) + xji_22*derxjm_211(ui))
             - vderiv_(1,0)*(derxjm_001(ui)*xji_11)
             - vderiv_(1,1)*(derxjm_011(ui)*xji_11)
             - vderiv_(1,2)*(derxjm_021(ui)*xji_11)
             - vderiv_(2,0)*(derxjm_211(ui)*xji_00 + derxjm_001(ui)*xji_21)
             - vderiv_(2,1)*(derxjm_211(ui)*xji_01 + derxjm_011(ui)*xji_21)
             - vderiv_(2,2)*(derxjm_211(ui)*xji_02 + derxjm_021(ui)*xji_21);
        v2 = - vderiv_(0,0)*(2*derxjm_021(ui)*xji_00 + 2*derxjm_001(ui)*xji_02 + xji_22*derxjm_201(ui) + xji_20*derxjm_221(ui))
             - vderiv_(0,1)*(2*derxjm_011(ui)*xji_02 + 2*derxjm_021(ui)*xji_01 + xji_21*derxjm_221(ui) + xji_22*derxjm_211(ui))
             - vderiv_(0,2)*(2*derxjm_021(ui)*xji_02 + 2*derxjm_021(ui)*xji_02 + xji_22*derxjm_221(ui) + xji_22*derxjm_221(ui))
             - vderiv_(1,0)*(derxjm_001(ui)*xji_12)
             - vderiv_(1,1)*(derxjm_011(ui)*xji_12)
             - vderiv_(1,2)*(derxjm_021(ui)*xji_12)
             - vderiv_(2,0)*(derxjm_221(ui)*xji_00 + derxjm_001(ui)*xji_22)
             - vderiv_(2,1)*(derxjm_221(ui)*xji_01 + derxjm_011(ui)*xji_22)
             - vderiv_(2,2)*(derxjm_221(ui)*xji_02 + derxjm_021(ui)*xji_22);

        for (int vi=0; vi<nen_; ++vi)
        {
          emesh(vi*4 + 0, ui*4 + 1) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(2*derxjm_002(ui)*xji_00 + 2*derxjm_002(ui)*xji_00 + xji_10*derxjm_102(ui) + xji_10*derxjm_102(ui))
             - vderiv_(0,1)*(2*derxjm_012(ui)*xji_00 + 2*derxjm_002(ui)*xji_01 + xji_11*derxjm_102(ui) + xji_10*derxjm_112(ui))
             - vderiv_(0,2)*(2*derxjm_022(ui)*xji_00 + 2*derxjm_002(ui)*xji_02 + xji_12*derxjm_102(ui) + xji_10*derxjm_122(ui))
             - vderiv_(1,0)*(derxjm_002(ui)*xji_10 + derxjm_102(ui)*xji_00)
             - vderiv_(1,1)*(derxjm_012(ui)*xji_10 + derxjm_102(ui)*xji_01)
             - vderiv_(1,2)*(derxjm_022(ui)*xji_10 + derxjm_102(ui)*xji_02)
             - vderiv_(2,0)*(derxjm_002(ui)*xji_20)
             - vderiv_(2,1)*(derxjm_012(ui)*xji_20)
             - vderiv_(2,2)*(derxjm_022(ui)*xji_20);
        v1 = - vderiv_(0,0)*(2*derxjm_012(ui)*xji_00 + 2*derxjm_002(ui)*xji_01 + xji_11*derxjm_102(ui) + xji_10*derxjm_112(ui))
             - vderiv_(0,1)*(2*derxjm_012(ui)*xji_01 + 2*derxjm_012(ui)*xji_01 + xji_11*derxjm_112(ui) + xji_11*derxjm_112(ui))
             - vderiv_(0,2)*(2*derxjm_012(ui)*xji_02 + 2*derxjm_022(ui)*xji_01 + xji_11*derxjm_122(ui) + xji_12*derxjm_112(ui))
             - vderiv_(1,0)*(derxjm_002(ui)*xji_11 + derxjm_112(ui)*xji_00)
             - vderiv_(1,1)*(derxjm_012(ui)*xji_11 + derxjm_112(ui)*xji_01)
             - vderiv_(1,2)*(derxjm_022(ui)*xji_11 + derxjm_112(ui)*xji_02)
             - vderiv_(2,0)*(derxjm_002(ui)*xji_21)
             - vderiv_(2,1)*(derxjm_012(ui)*xji_21)
             - vderiv_(2,2)*(derxjm_022(ui)*xji_21);
        v2 = - vderiv_(0,0)*(2*derxjm_022(ui)*xji_00 + 2*derxjm_002(ui)*xji_02 + xji_12*derxjm_102(ui) + xji_10*derxjm_122(ui))
             - vderiv_(0,1)*(2*derxjm_012(ui)*xji_02 + 2*derxjm_022(ui)*xji_01 + xji_11*derxjm_122(ui) + xji_12*derxjm_112(ui))
             - vderiv_(0,2)*(2*derxjm_022(ui)*xji_02 + 2*derxjm_022(ui)*xji_02 + xji_12*derxjm_122(ui) + xji_12*derxjm_122(ui))
             - vderiv_(1,0)*(derxjm_002(ui)*xji_12 + derxjm_122(ui)*xji_00)
             - vderiv_(1,1)*(derxjm_012(ui)*xji_12 + derxjm_122(ui)*xji_01)
             - vderiv_(1,2)*(derxjm_022(ui)*xji_12 + derxjm_122(ui)*xji_02)
             - vderiv_(2,0)*(derxjm_002(ui)*xji_22)
             - vderiv_(2,1)*(derxjm_012(ui)*xji_22)
             - vderiv_(2,2)*(derxjm_022(ui)*xji_22);

        for (int vi=0; vi<nen_; ++vi)
        {
          emesh(vi*4 + 0, ui*4 + 2) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(derxjm_100(ui)*xji_00)
             - vderiv_(0,1)*(derxjm_110(ui)*xji_00)
             - vderiv_(0,2)*(derxjm_120(ui)*xji_00)
             - vderiv_(1,0)*(2*xji_10*derxjm_100(ui) + 2*xji_10*derxjm_100(ui) + xji_20*derxjm_200(ui) + xji_20*derxjm_200(ui))
             - vderiv_(1,1)*(2*xji_11*derxjm_100(ui) + 2*xji_10*derxjm_110(ui) + xji_21*derxjm_200(ui) + xji_20*derxjm_210(ui))
             - vderiv_(1,2)*(2*xji_12*derxjm_100(ui) + 2*xji_10*derxjm_120(ui) + xji_22*derxjm_200(ui) + xji_20*derxjm_220(ui))
             - vderiv_(2,0)*(derxjm_200(ui)*xji_10 + derxjm_100(ui)*xji_20)
             - vderiv_(2,1)*(derxjm_200(ui)*xji_11 + derxjm_110(ui)*xji_20)
             - vderiv_(2,2)*(derxjm_200(ui)*xji_12 + derxjm_120(ui)*xji_20);
        v1 = - vderiv_(0,0)*(derxjm_100(ui)*xji_01)
             - vderiv_(0,1)*(derxjm_110(ui)*xji_01)
             - vderiv_(0,2)*(derxjm_120(ui)*xji_01)
             - vderiv_(1,0)*(2*xji_10*derxjm_110(ui) + 2*xji_11*derxjm_100(ui) + xji_20*derxjm_210(ui) + xji_21*derxjm_200(ui))
             - vderiv_(1,1)*(2*xji_11*derxjm_110(ui) + 2*xji_11*derxjm_110(ui) + xji_21*derxjm_210(ui) + xji_21*derxjm_210(ui))
             - vderiv_(1,2)*(2*xji_12*derxjm_110(ui) + 2*xji_11*derxjm_120(ui) + xji_22*derxjm_210(ui) + xji_21*derxjm_220(ui))
             - vderiv_(2,0)*(derxjm_210(ui)*xji_10 + derxjm_100(ui)*xji_21)
             - vderiv_(2,1)*(derxjm_210(ui)*xji_11 + derxjm_110(ui)*xji_21)
             - vderiv_(2,2)*(derxjm_210(ui)*xji_12 + derxjm_120(ui)*xji_21);
        v2 = - vderiv_(0,0)*(derxjm_100(ui)*xji_02)
             - vderiv_(0,1)*(derxjm_110(ui)*xji_02)
             - vderiv_(0,2)*(derxjm_120(ui)*xji_02)
             - vderiv_(1,0)*(2*xji_10*derxjm_120(ui) + 2*xji_12*derxjm_100(ui) + xji_20*derxjm_220(ui) + xji_22*derxjm_200(ui))
             - vderiv_(1,1)*(2*xji_11*derxjm_120(ui) + 2*xji_12*derxjm_110(ui) + xji_21*derxjm_220(ui) + xji_22*derxjm_210(ui))
             - vderiv_(1,2)*(2*xji_12*derxjm_120(ui) + 2*xji_12*derxjm_120(ui) + xji_22*derxjm_220(ui) + xji_22*derxjm_220(ui))
             - vderiv_(2,0)*(derxjm_220(ui)*xji_10 + derxjm_100(ui)*xji_22)
             - vderiv_(2,1)*(derxjm_220(ui)*xji_11 + derxjm_110(ui)*xji_22)
             - vderiv_(2,2)*(derxjm_220(ui)*xji_12 + derxjm_120(ui)*xji_22);

        for (int vi=0; vi<nen_; ++vi)
        {
          emesh(vi*4 + 1, ui*4 + 0) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(derxjm_001(ui)*xji_10)
             - vderiv_(0,1)*(derxjm_001(ui)*xji_11)
             - vderiv_(0,2)*(derxjm_001(ui)*xji_12)
             - vderiv_(1,0)*(xji_00*derxjm_001(ui) + xji_00*derxjm_001(ui) + xji_20*derxjm_201(ui) + xji_20*derxjm_201(ui))
             - vderiv_(1,1)*(xji_01*derxjm_001(ui) + xji_00*derxjm_011(ui) + xji_21*derxjm_201(ui) + xji_20*derxjm_211(ui))
             - vderiv_(1,2)*(xji_02*derxjm_001(ui) + xji_00*derxjm_021(ui) + xji_22*derxjm_201(ui) + xji_20*derxjm_221(ui))
             - vderiv_(2,0)*(derxjm_201(ui)*xji_10)
             - vderiv_(2,1)*(derxjm_201(ui)*xji_11)
             - vderiv_(2,2)*(derxjm_201(ui)*xji_12);
        v1 = - vderiv_(0,0)*(derxjm_011(ui)*xji_10)
             - vderiv_(0,1)*(derxjm_011(ui)*xji_11)
             - vderiv_(0,2)*(derxjm_011(ui)*xji_12)
             - vderiv_(1,0)*(xji_00*derxjm_011(ui) + xji_01*derxjm_001(ui) + xji_20*derxjm_211(ui) + xji_21*derxjm_201(ui))
             - vderiv_(1,1)*(xji_01*derxjm_011(ui) + xji_01*derxjm_011(ui) + xji_21*derxjm_211(ui) + xji_21*derxjm_211(ui))
             - vderiv_(1,2)*(xji_02*derxjm_011(ui) + xji_01*derxjm_021(ui) + xji_22*derxjm_211(ui) + xji_21*derxjm_221(ui))
             - vderiv_(2,0)*(derxjm_211(ui)*xji_10)
             - vderiv_(2,1)*(derxjm_211(ui)*xji_11)
             - vderiv_(2,2)*(derxjm_211(ui)*xji_12);
        v2 = - vderiv_(0,0)*(derxjm_021(ui)*xji_10)
             - vderiv_(0,1)*(derxjm_021(ui)*xji_11)
             - vderiv_(0,2)*(derxjm_021(ui)*xji_12)
             - vderiv_(1,0)*(xji_00*derxjm_021(ui) + xji_02*derxjm_001(ui) + xji_20*derxjm_221(ui) + xji_22*derxjm_201(ui))
             - vderiv_(1,1)*(xji_01*derxjm_021(ui) + xji_02*derxjm_011(ui) + xji_21*derxjm_221(ui) + xji_22*derxjm_211(ui))
             - vderiv_(1,2)*(xji_02*derxjm_021(ui) + xji_02*derxjm_021(ui) + xji_22*derxjm_221(ui) + xji_22*derxjm_221(ui))
             - vderiv_(2,0)*(derxjm_221(ui)*xji_10)
             - vderiv_(2,1)*(derxjm_221(ui)*xji_11)
             - vderiv_(2,2)*(derxjm_221(ui)*xji_12);

        for (int vi=0; vi<nen_; ++vi)
        {
          emesh(vi*4 + 1, ui*4 + 1) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(derxjm_002(ui)*xji_10 + derxjm_102(ui)*xji_00)
             - vderiv_(0,1)*(derxjm_002(ui)*xji_11 + derxjm_112(ui)*xji_00)
             - vderiv_(0,2)*(derxjm_002(ui)*xji_12 + derxjm_122(ui)*xji_00)
             - vderiv_(1,0)*(xji_00*derxjm_002(ui) + xji_00*derxjm_002(ui) + 2*xji_10*derxjm_102(ui) + 2*xji_10*derxjm_102(ui))
             - vderiv_(1,1)*(xji_01*derxjm_002(ui) + xji_00*derxjm_012(ui) + 2*xji_11*derxjm_102(ui) + 2*xji_10*derxjm_112(ui))
             - vderiv_(1,2)*(xji_02*derxjm_002(ui) + xji_00*derxjm_022(ui) + 2*xji_12*derxjm_102(ui) + 2*xji_10*derxjm_122(ui))
             - vderiv_(2,0)*(derxjm_102(ui)*xji_20)
             - vderiv_(2,1)*(derxjm_112(ui)*xji_20)
             - vderiv_(2,2)*(derxjm_122(ui)*xji_20);
        v1 = - vderiv_(0,0)*(derxjm_012(ui)*xji_10 + derxjm_102(ui)*xji_01)
             - vderiv_(0,1)*(derxjm_012(ui)*xji_11 + derxjm_112(ui)*xji_01)
             - vderiv_(0,2)*(derxjm_012(ui)*xji_12 + derxjm_122(ui)*xji_01)
             - vderiv_(1,0)*(xji_00*derxjm_012(ui) + xji_01*derxjm_002(ui) + 2*xji_10*derxjm_112(ui) + 2*xji_11*derxjm_102(ui))
             - vderiv_(1,1)*(xji_01*derxjm_012(ui) + xji_01*derxjm_012(ui) + 2*xji_11*derxjm_112(ui) + 2*xji_11*derxjm_112(ui))
             - vderiv_(1,2)*(xji_02*derxjm_012(ui) + xji_01*derxjm_022(ui) + 2*xji_12*derxjm_112(ui) + 2*xji_11*derxjm_122(ui))
             - vderiv_(2,0)*(derxjm_102(ui)*xji_21)
             - vderiv_(2,1)*(derxjm_112(ui)*xji_21)
             - vderiv_(2,2)*(derxjm_122(ui)*xji_21);
        v2 = - vderiv_(0,0)*(derxjm_022(ui)*xji_10 + derxjm_102(ui)*xji_02)
             - vderiv_(0,1)*(derxjm_022(ui)*xji_11 + derxjm_112(ui)*xji_02)
             - vderiv_(0,2)*(derxjm_022(ui)*xji_12 + derxjm_122(ui)*xji_02)
             - vderiv_(1,0)*(xji_00*derxjm_022(ui) + xji_02*derxjm_002(ui) + 2*xji_10*derxjm_122(ui) + 2*xji_12*derxjm_102(ui))
             - vderiv_(1,1)*(xji_01*derxjm_022(ui) + xji_02*derxjm_012(ui) + 2*xji_11*derxjm_122(ui) + 2*xji_12*derxjm_112(ui))
             - vderiv_(1,2)*(xji_02*derxjm_022(ui) + xji_02*derxjm_022(ui) + 2*xji_12*derxjm_122(ui) + 2*xji_12*derxjm_122(ui))
             - vderiv_(2,0)*(derxjm_102(ui)*xji_22)
             - vderiv_(2,1)*(derxjm_112(ui)*xji_22)
             - vderiv_(2,2)*(derxjm_122(ui)*xji_22);

        for (int vi=0; vi<nen_; ++vi)
        {
          emesh(vi*4 + 1, ui*4 + 2) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(derxjm_200(ui)*xji_00)
             - vderiv_(0,1)*(derxjm_210(ui)*xji_00)
             - vderiv_(0,2)*(derxjm_220(ui)*xji_00)
             - vderiv_(1,0)*(derxjm_200(ui)*xji_10 + derxjm_100(ui)*xji_20)
             - vderiv_(1,1)*(derxjm_210(ui)*xji_10 + derxjm_100(ui)*xji_21)
             - vderiv_(1,2)*(derxjm_220(ui)*xji_10 + derxjm_100(ui)*xji_22)
             - vderiv_(2,0)*(xji_10*derxjm_100(ui) + xji_10*derxjm_100(ui) + 2*xji_20*derxjm_200(ui) + 2*xji_20*derxjm_200(ui))
             - vderiv_(2,1)*(xji_11*derxjm_100(ui) + xji_10*derxjm_110(ui) + 2*xji_21*derxjm_200(ui) + 2*xji_20*derxjm_210(ui))
             - vderiv_(2,2)*(xji_12*derxjm_100(ui) + xji_10*derxjm_120(ui) + 2*xji_22*derxjm_200(ui) + 2*xji_20*derxjm_220(ui));
        v1 = - vderiv_(0,0)*(derxjm_200(ui)*xji_01)
             - vderiv_(0,1)*(derxjm_210(ui)*xji_01)
             - vderiv_(0,2)*(derxjm_220(ui)*xji_01)
             - vderiv_(1,0)*(derxjm_200(ui)*xji_11 + derxjm_110(ui)*xji_20)
             - vderiv_(1,1)*(derxjm_210(ui)*xji_11 + derxjm_110(ui)*xji_21)
             - vderiv_(1,2)*(derxjm_220(ui)*xji_11 + derxjm_110(ui)*xji_22)
             - vderiv_(2,0)*(xji_10*derxjm_110(ui) + xji_11*derxjm_100(ui) + 2*xji_20*derxjm_210(ui) + 2*xji_21*derxjm_200(ui))
             - vderiv_(2,1)*(xji_11*derxjm_110(ui) + xji_11*derxjm_110(ui) + 2*xji_21*derxjm_210(ui) + 2*xji_21*derxjm_210(ui))
             - vderiv_(2,2)*(xji_12*derxjm_110(ui) + xji_11*derxjm_120(ui) + 2*xji_22*derxjm_210(ui) + 2*xji_21*derxjm_220(ui));
        v2 = - vderiv_(0,0)*(derxjm_200(ui)*xji_02)
             - vderiv_(0,1)*(derxjm_210(ui)*xji_02)
             - vderiv_(0,2)*(derxjm_220(ui)*xji_02)
             - vderiv_(1,0)*(derxjm_200(ui)*xji_12 + derxjm_120(ui)*xji_20)
             - vderiv_(1,1)*(derxjm_210(ui)*xji_12 + derxjm_120(ui)*xji_21)
             - vderiv_(1,2)*(derxjm_220(ui)*xji_12 + derxjm_120(ui)*xji_22)
             - vderiv_(2,0)*(xji_10*derxjm_120(ui) + xji_12*derxjm_100(ui) + 2*xji_20*derxjm_220(ui) + 2*xji_22*derxjm_200(ui))
             - vderiv_(2,1)*(xji_11*derxjm_120(ui) + xji_12*derxjm_110(ui) + 2*xji_21*derxjm_220(ui) + 2*xji_22*derxjm_210(ui))
             - vderiv_(2,2)*(xji_12*derxjm_120(ui) + xji_12*derxjm_120(ui) + 2*xji_22*derxjm_220(ui) + 2*xji_22*derxjm_220(ui));

        for (int vi=0; vi<nen_; ++vi)
        {
          emesh(vi*4 + 2, ui*4 + 0) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(derxjm_201(ui)*xji_00 + derxjm_001(ui)*xji_20)
             - vderiv_(0,1)*(derxjm_211(ui)*xji_00 + derxjm_001(ui)*xji_21)
             - vderiv_(0,2)*(derxjm_221(ui)*xji_00 + derxjm_001(ui)*xji_22)
             - vderiv_(1,0)*(derxjm_201(ui)*xji_10)
             - vderiv_(1,1)*(derxjm_211(ui)*xji_10)
             - vderiv_(1,2)*(derxjm_221(ui)*xji_10)
             - vderiv_(2,0)*(xji_00*derxjm_001(ui) + xji_00*derxjm_001(ui) + 2*xji_20*derxjm_201(ui) + 2*xji_20*derxjm_201(ui))
             - vderiv_(2,1)*(xji_01*derxjm_001(ui) + xji_00*derxjm_011(ui) + 2*xji_21*derxjm_201(ui) + 2*xji_20*derxjm_211(ui))
             - vderiv_(2,2)*(xji_02*derxjm_001(ui) + xji_00*derxjm_021(ui) + 2*xji_22*derxjm_201(ui) + 2*xji_20*derxjm_221(ui));
        v1 = - vderiv_(0,0)*(derxjm_201(ui)*xji_01 + derxjm_011(ui)*xji_20)
             - vderiv_(0,1)*(derxjm_211(ui)*xji_01 + derxjm_011(ui)*xji_21)
             - vderiv_(0,2)*(derxjm_221(ui)*xji_01 + derxjm_011(ui)*xji_22)
             - vderiv_(1,0)*(derxjm_201(ui)*xji_11)
             - vderiv_(1,1)*(derxjm_211(ui)*xji_11)
             - vderiv_(1,2)*(derxjm_221(ui)*xji_11)
             - vderiv_(2,0)*(xji_00*derxjm_011(ui) + xji_01*derxjm_001(ui) + 2*xji_20*derxjm_211(ui) + 2*xji_21*derxjm_201(ui))
             - vderiv_(2,1)*(xji_01*derxjm_011(ui) + xji_01*derxjm_011(ui) + 2*xji_21*derxjm_211(ui) + 2*xji_21*derxjm_211(ui))
             - vderiv_(2,2)*(xji_02*derxjm_011(ui) + xji_01*derxjm_021(ui) + 2*xji_22*derxjm_211(ui) + 2*xji_21*derxjm_221(ui));
        v2 = - vderiv_(0,0)*(derxjm_201(ui)*xji_02 + derxjm_021(ui)*xji_20)
             - vderiv_(0,1)*(derxjm_211(ui)*xji_02 + derxjm_021(ui)*xji_21)
             - vderiv_(0,2)*(derxjm_221(ui)*xji_02 + derxjm_021(ui)*xji_22)
             - vderiv_(1,0)*(derxjm_201(ui)*xji_12)
             - vderiv_(1,1)*(derxjm_211(ui)*xji_12)
             - vderiv_(1,2)*(derxjm_221(ui)*xji_12)
             - vderiv_(2,0)*(xji_00*derxjm_021(ui) + xji_02*derxjm_001(ui) + 2*xji_20*derxjm_221(ui) + 2*xji_22*derxjm_201(ui))
             - vderiv_(2,1)*(xji_01*derxjm_021(ui) + xji_02*derxjm_011(ui) + 2*xji_21*derxjm_221(ui) + 2*xji_22*derxjm_211(ui))
             - vderiv_(2,2)*(xji_02*derxjm_021(ui) + xji_02*derxjm_021(ui) + 2*xji_22*derxjm_221(ui) + 2*xji_22*derxjm_221(ui));

        for (int vi=0; vi<nen_; ++vi)
        {
          emesh(vi*4 + 2, ui*4 + 1) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }

        ////////////////////////////////////////////////////////////////

        v0 = - vderiv_(0,0)*(derxjm_002(ui)*xji_20)
             - vderiv_(0,1)*(derxjm_002(ui)*xji_21)
             - vderiv_(0,2)*(derxjm_002(ui)*xji_22)
             - vderiv_(1,0)*(derxjm_102(ui)*xji_20)
             - vderiv_(1,1)*(derxjm_102(ui)*xji_21)
             - vderiv_(1,2)*(derxjm_102(ui)*xji_22)
             - vderiv_(2,0)*(xji_00*derxjm_002(ui) + xji_00*derxjm_002(ui) + xji_10*derxjm_102(ui) + xji_10*derxjm_102(ui))
             - vderiv_(2,1)*(xji_01*derxjm_002(ui) + xji_00*derxjm_012(ui) + xji_11*derxjm_102(ui) + xji_10*derxjm_112(ui))
             - vderiv_(2,2)*(xji_02*derxjm_002(ui) + xji_00*derxjm_022(ui) + xji_12*derxjm_102(ui) + xji_10*derxjm_122(ui));
        v1 = - vderiv_(0,0)*(derxjm_012(ui)*xji_20)
             - vderiv_(0,1)*(derxjm_012(ui)*xji_21)
             - vderiv_(0,2)*(derxjm_012(ui)*xji_22)
             - vderiv_(1,0)*(derxjm_112(ui)*xji_20)
             - vderiv_(1,1)*(derxjm_112(ui)*xji_21)
             - vderiv_(1,2)*(derxjm_112(ui)*xji_22)
             - vderiv_(2,0)*(xji_00*derxjm_012(ui) + xji_01*derxjm_002(ui) + xji_10*derxjm_112(ui) + xji_11*derxjm_102(ui))
             - vderiv_(2,1)*(xji_01*derxjm_012(ui) + xji_01*derxjm_012(ui) + xji_11*derxjm_112(ui) + xji_11*derxjm_112(ui))
             - vderiv_(2,2)*(xji_02*derxjm_012(ui) + xji_01*derxjm_022(ui) + xji_12*derxjm_112(ui) + xji_11*derxjm_122(ui));
        v2 = - vderiv_(0,0)*(derxjm_022(ui)*xji_20)
             - vderiv_(0,1)*(derxjm_022(ui)*xji_21)
             - vderiv_(0,2)*(derxjm_022(ui)*xji_22)
             - vderiv_(1,0)*(derxjm_122(ui)*xji_20)
             - vderiv_(1,1)*(derxjm_122(ui)*xji_21)
             - vderiv_(1,2)*(derxjm_122(ui)*xji_22)
             - vderiv_(2,0)*(xji_00*derxjm_022(ui) + xji_02*derxjm_002(ui) + xji_10*derxjm_122(ui) + xji_12*derxjm_102(ui))
             - vderiv_(2,1)*(xji_01*derxjm_022(ui) + xji_02*derxjm_012(ui) + xji_11*derxjm_122(ui) + xji_12*derxjm_112(ui))
             - vderiv_(2,2)*(xji_02*derxjm_022(ui) + xji_02*derxjm_022(ui) + xji_12*derxjm_122(ui) + xji_12*derxjm_122(ui));

        for (int vi=0; vi<nen_; ++vi)
        {
          emesh(vi*4 + 2, ui*4 + 2) += v*(deriv_(0,vi)*v0 + deriv_(1,vi)*v1 + deriv_(2,vi)*v2);
        }
      }


      // pressure
      for (int vi=0; vi<nen_; ++vi)
      {
        double v = press*timefacfac/det_;
        for (int ui=0; ui<nen_; ++ui)
        {
          emesh(vi*4    , ui*4 + 1) += v*(deriv_(0, vi)*derxjm_(0,0,1,ui) + deriv_(1, vi)*derxjm_(0,1,1,ui) + deriv_(2, vi)*derxjm_(0,2,1,ui)) ;
          emesh(vi*4    , ui*4 + 2) += v*(deriv_(0, vi)*derxjm_(0,0,2,ui) + deriv_(1, vi)*derxjm_(0,1,2,ui) + deriv_(2, vi)*derxjm_(0,2,2,ui)) ;

          emesh(vi*4 + 1, ui*4 + 0) += v*(deriv_(0, vi)*derxjm_(1,0,0,ui) + deriv_(1, vi)*derxjm_(1,1,0,ui) + deriv_(2, vi)*derxjm_(1,2,0,ui)) ;
          emesh(vi*4 + 1, ui*4 + 2) += v*(deriv_(0, vi)*derxjm_(1,0,2,ui) + deriv_(1, vi)*derxjm_(1,1,2,ui) + deriv_(2, vi)*derxjm_(1,2,2,ui)) ;

          emesh(vi*4 + 2, ui*4 + 0) += v*(deriv_(0, vi)*derxjm_(2,0,0,ui) + deriv_(1, vi)*derxjm_(2,1,0,ui) + deriv_(2, vi)*derxjm_(2,2,0,ui)) ;
          emesh(vi*4 + 2, ui*4 + 1) += v*(deriv_(0, vi)*derxjm_(2,0,1,ui) + deriv_(1, vi)*derxjm_(2,1,1,ui) + deriv_(2, vi)*derxjm_(2,2,1,ui)) ;
        }
      }

      // div u
      for (int vi=0; vi<nen_; ++vi)
      {
        double v = timefacfac/det_*funct_(vi,0);
        for (int ui=0; ui<nen_; ++ui)
        {
          emesh(vi*4 + 3, ui*4 + 0) += v*(
            + vderiv_(1, 0)*derxjm_(0,0,1,ui) + vderiv_(1, 1)*derxjm_(0,1,1,ui) + vderiv_(1, 2)*derxjm_(0,2,1,ui)
            + vderiv_(2, 0)*derxjm_(0,0,2,ui) + vderiv_(2, 1)*derxjm_(0,1,2,ui) + vderiv_(2, 2)*derxjm_(0,2,2,ui)
            ) ;

          emesh(vi*4 + 3, ui*4 + 1) += v*(
            + vderiv_(0, 0)*derxjm_(1,0,0,ui) + vderiv_(0, 1)*derxjm_(1,1,0,ui) + vderiv_(0, 2)*derxjm_(1,2,0,ui)
            + vderiv_(2, 0)*derxjm_(1,0,2,ui) + vderiv_(2, 1)*derxjm_(1,1,2,ui) + vderiv_(2, 2)*derxjm_(1,2,2,ui)
            ) ;

          emesh(vi*4 + 3, ui*4 + 2) += v*(
            + vderiv_(0, 0)*derxjm_(2,0,0,ui) + vderiv_(0, 1)*derxjm_(2,1,0,ui) + vderiv_(0, 2)*derxjm_(2,2,0,ui)
            + vderiv_(1, 0)*derxjm_(2,0,1,ui) + vderiv_(1, 1)*derxjm_(2,1,1,ui) + vderiv_(1, 2)*derxjm_(2,2,1,ui)
            ) ;
        }
      }

    }
  } // loop gausspoints
}

#endif
#endif

