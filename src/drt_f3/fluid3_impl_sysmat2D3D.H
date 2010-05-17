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
  Fluid3*                                       ele,
  const LINALG::Matrix<nsd_,nen_>&              evelaf,
  const LINALG::Matrix<nsd_,nen_>&              eveln,
  const LINALG::Matrix<nsd_,nen_>&              fsevelaf,
  const LINALG::Matrix<nen_,1>&                 epreaf,
  const LINALG::Matrix<nsd_,nen_>&              eaccam,
  const LINALG::Matrix<nen_,1>&                 escaaf,
  const LINALG::Matrix<nen_,1>&                 escaam,
  const LINALG::Matrix<nen_,1>&                 escadtam,
  const LINALG::Matrix<nsd_,nen_>&              emhist,
  const LINALG::Matrix<nsd_,nen_>&              edispnp,
  const LINALG::Matrix<nsd_,nen_>&              egridv,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  estif,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  emesh,
  LINALG::Matrix<(nsd_+1)*nen_,1>&              eforce,
  const double                                  thermpressaf,
  const double                                  thermpressam,
  const double                                  thermpressdtam,
  Teuchos::RCP<const MAT::Material>             material,
  const double                                  timefac,
  double&                                       Cs,
  double&                                       Cs_delta_sq,
  double&                                       l_tau
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

     // set visc_old to zero
     visc_old_.Clear();

     double prefac;
     if(physicaltype_ == INPAR::FLUID::loma)
     {
       prefac = 1.0/3.0;
       derxy2_.Scale(prefac);
     }
     else prefac = 1.0;

     for (int inode=0; inode<nen_; ++inode)
     {
       if (nsd_ == 3)
       {
         double sum = (derxy2_(0,inode)+derxy2_(1,inode)+derxy2_(2,inode))/prefac;
         viscs2_(0,inode) = 0.5 * (sum + derxy2_(0,inode));
         viscs2_(1,inode) = 0.5 *  derxy2_(3,inode);
         viscs2_(2,inode) = 0.5 *  derxy2_(4,inode);
         viscs2_(3,inode) = 0.5 *  derxy2_(3,inode);
         viscs2_(4,inode) = 0.5 * (sum + derxy2_(1,inode));
         viscs2_(5,inode) = 0.5 *  derxy2_(5,inode);
         viscs2_(6,inode) = 0.5 *  derxy2_(4,inode);
         viscs2_(7,inode) = 0.5 *  derxy2_(5,inode);
         viscs2_(8,inode) = 0.5 * (sum + derxy2_(2,inode));
       }
       else if (nsd_==2)
       {
         double sum = (derxy2_(0,inode)+derxy2_(1,inode))/prefac;
         viscs2_(0,inode) = 0.5 * (sum + derxy2_(0,inode));
         viscs2_(1,inode) = 0.5 * derxy2_(2,inode);
         viscs2_(2,inode) = 0.5 * derxy2_(2,inode);
         viscs2_(3,inode) = 0.5 * (sum + derxy2_(1,inode));
       }
       else dserror("Epsilon(N) is not implemented for the 1D case");

       for (int idim=0; idim<nsd_; ++idim)
       {
         for (int jdim=0; jdim<nsd_; ++jdim)
         {
           visc_old_(idim) += viscs2_(jdim+(idim*nsd_),inode)*evelaf(jdim,inode);
         }
         //visc_old_(0) += viscs2_(0,i)*evelaf(0,i)+viscs2_(1,i)*evelaf(1,i)+viscs2_(2,i)*evelaf(2,i);
         //visc_old_(1) += viscs2_(3,i)*evelaf(0,i)+viscs2_(4,i)*evelaf(1,i)+viscs2_(5,i)*evelaf(2,i);
         //visc_old_(2) += viscs2_(6,i)*evelaf(0,i)+viscs2_(7,i)*evelaf(1,i)+viscs2_(8,i)*evelaf(2,i);
       }
     }
    }
    else
    {
      viscs2_.Clear();
      visc_old_.Clear();
    }
  

#if 0
    if (is_higher_order_ele_)
    {
      if(nsd_ == 3)
      {

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
        viscs2_(2,0) = 0.5 * derxy2_(2,0);
        viscs2_(3,0) = 0.5 * (sum + derxy2_(1,0));
      
        /* viscous term  div epsilon(u_old) */
        visc_old_(0) = viscs2_(0,0)*evelaf(0,0)+viscs2_(1,0)*evelaf(1,0);
        visc_old_(1) = viscs2_(1,0)*evelaf(0,0)+viscs2_(3,0)*evelaf(1,0);

        for (int i=1; i<nen_; ++i)
        {
          double sum = (derxy2_(0,i)+derxy2_(1,i))/prefac;

          viscs2_(0,i) = 0.5 * (sum + derxy2_(0,i));
          viscs2_(1,i) = 0.5 * derxy2_(2,i);
          viscs2_(2,i) = 0.5 * derxy2_(2,i);
          viscs2_(3,i) = 0.5 * (sum + derxy2_(1,i));

          /* viscous term  div epsilon(u_old) */
          visc_old_(0) += viscs2_(0,i)*evelaf(0,i)+viscs2_(1,i)*evelaf(1,i);
          visc_old_(1) += viscs2_(1,i)*evelaf(0,i)+viscs2_(3,i)*evelaf(1,i);
        }
      }
    }
    else
    {
      viscs2_.Clear();
      visc_old_.Clear();
    }
#endif

    // convective term from previous iteration
    conv_old_.Multiply(vderxy_,convvelint_);

    // compute convective operator
    conv_c_.MultiplyTN(derxy_,convvelint_);

    // velocity divergence from previous iteration

    vdiv_ = 0.0;
    for (int idim = 0; idim <nsd_; ++idim)
    {
      vdiv_ += vderxy_(idim, idim);
    }
    //vdiv_ = vderxy_(0, 0) + vderxy_(1, 1) + vderxy_(2, 2);

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
      for (int rr=0;rr<nsd_;++rr)
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


        /*
                                                                             
                 /                                                dp   \     
                |         1     / dT     /         \   \     1      th  |     
                |    q , --- * | ---- + | u o nabla | T | - --- * ----  |        
                |         T     \ dt     \         /   /    p      dt   |      
                 \                                           th        /    
                        +---------------------------------------------+     
                                          rhscon_                                                                             
        */

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
        for (int rr=0;rr<nsd_;++rr)
        {
          momres_old_(rr) = densaf_*velint_(rr)+timefac*(densaf_*conv_old_(rr)+gradp_(rr)-2*visceff_*visc_old_(rr))-rhsmom_(rr);
        }
      }
      else
      {
        // momres_old = theta ( ... ) -  bodyforce_
        for (int rr=0;rr<nsd_;++rr)
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
        //const double vdivn = vderxyn_(0, 0) + vderxyn_(1, 1) + vderxyn_(2, 2);

        double vdivn = 0.0;
        for (int idim = 0; idim<nsd_; ++idim)
        {
          vdivn += vderxyn_(idim,idim);
        }

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
        const double fac_densam=fac_*densam_;

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = numdofpernode_*ui;
          //const int fui   = 4*ui;
          //const int fuip  = fui+1;
          //const int fuipp = fui+2;
          
          const double v = fac_densam*funct_(ui);
          
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = numdofpernode_*vi;
            //const int fvi   = 4*vi;
            //const int fvip  = fvi+1;
            //const int fvipp = fvi+2;
            /* inertia (contribution to mass matrix) */
            /*
              /                \
              |                |
              |    rho*Du , v  |
              |                |
              \                /
            */

            double v2 = v*funct_(vi) ;
            
            for (int idim = 0; idim <nsd_; ++idim)
            {
              estif(fvi+idim,fui+idim) += v2;
              //estif(fvi  , fui  ) += v2;
              //estif(fvip , fuip ) += v2;
              //estif(fvipp, fuipp) += v2;
            }
          }
        }  // end for (idim)
      }  // end if(stationary)


      // convective part
      {
        const double timefacfac_densaf=timefacfac*densaf_;
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = numdofpernode_*ui;
          //const int fui   = 4*ui;
          //const int fuip  = fui+1;
          //const int fuipp = fui+2;
          
          const double v = timefacfac_densaf*(conv_c_(ui)+sgconv_c_(ui));

          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = numdofpernode_*vi;
            //const int fvi   = 4*vi;
            //const int fvip  = fvi+1;
            //const int fvipp = fvi+2;

            /* convection, convective part (convective form) */
            /*
            /                               \
            |  /       n+1        \         |
            | |   rho*u   o nabla | Du , v  |
            |  \      (i)        /          |
            \                              /
            */
            double v2 = v*funct_(vi) ;

            for (int idim = 0; idim <nsd_; ++idim)
            {
              estif(fvi+idim,fui+idim) += v2;
            //estif(fvi  , fui  ) += v2;
            //estif(fvip , fuip ) += v2;
            //estif(fvipp, fuipp) += v2;
            }  // end for (idim)
          } // vi
        } // ui
      }

      if (is_newton_)
      {
        const double timefacfac_densaf=timefacfac*densaf_;
        
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = numdofpernode_*vi;
            //const int fvi   = 4*vi;
            //const int fvip  = fvi+1;
            //const int fvipp = fvi+2;
            const double v = timefacfac_densaf*funct_(vi);

            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = numdofpernode_*ui;
              //const int fui   = 4*ui;
              //const int fuip  = fui+1;
              //const int fuipp = fui+2;
              const double v2 = v*funct_(ui);
              /*  convection, reactive part (convective form)
              /                                 \
              |  /                \   n+1       |
              | |  rho*Du o nabla | u     , v   |
              |  \                /   (i)       |
              \                                /
              */
              for (int idim = 0; idim <nsd_; ++idim)
              {
                for (int jdim= 0; jdim<nsd_;++jdim)
                {
                  estif(fvi+idim,fui+jdim) += v2*vderxy_(idim, jdim);
                } // end for (jdim)
              } // end for (idim)
              //estif(fvi,   fui)   += v2*vderxy_(0, 0) ;
              //estif(fvi,   fuip)  += v2*vderxy_(0, 1) ;
              //estif(fvi,   fuipp) += v2*vderxy_(0, 2) ;
              //estif(fvip,  fui)   += v2*vderxy_(1, 0) ;
              //estif(fvip,  fuip)  += v2*vderxy_(1, 1) ;
              //estif(fvip,  fuipp) += v2*vderxy_(1, 2) ;
              //estif(fvipp, fui)   += v2*vderxy_(2, 0) ;
              //estif(fvipp, fuip)  += v2*vderxy_(2, 1) ;
              //estif(fvipp, fuipp) += v2*vderxy_(2, 2) ;
            } // ui
        } //vi
      }  // end if(Newton)

      // inertia terms on the right hand side for instationary fluids
      if (is_stationary_ == false)
      {
        if (is_genalpha_)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi = numdofpernode_*vi;
            //const int fvi = 4*vi;
            /* inertia term on right-hand side for generalized-alpha scheme */
            const double v = -fac_*densam_*funct_(vi);
            
            for (int idim = 0; idim <nsd_; ++idim)
            {
              eforce(fvi+idim) += v*accint_(idim) ;
            }
            //eforce(fvi    ) += v*accint_(0) ;
            //eforce(fvi + 1) += v*accint_(1) ;
            //eforce(fvi + 2) += v*accint_(2) ;
          }
        }
        else
        {
          for (int idim = 0; idim <nsd_; ++idim)
          {
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi = numdofpernode_*vi + idim;
              //const int fvi = 4*vi;
              /* inertia term on right-hand side for one-step-theta/BDF2 schem */
              const double v = -fac_*densaf_*funct_(vi);
              eforce(fvi    ) += v*velint_(idim) ;
              //eforce(fvi    ) += v*velint_(0) ;
              //eforce(fvi + 1) += v*velint_(1) ;
              //eforce(fvi + 2) += v*velint_(2) ;
            }
          }  // end for (idim)
        }
      }  // end if(stationary)

      for (int idim = 0; idim <nsd_; ++idim)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = numdofpernode_*vi + idim;
          //const int fvi   = 4*vi;
          /* convection (convective form) on right-hand side */
          double v = -rhsfac*densaf_*funct_(vi);
          eforce(fvi    ) += v*conv_old_(idim) ;
          //eforce(fvi    ) += v*conv_old_(0) ;
          //eforce(fvi + 1) += v*conv_old_(1) ;
          //eforce(fvi + 2) += v*conv_old_(2) ;
        }
      }  // end for(idim)

      //----------------------------------------------------------------------
      // computation of additions to convection term (convective and
      // reactive part) for conservative form of convection term including
      // right-hand-side contribution
      //----------------------------------------------------------------------

      if (is_conservative_)
      {
        for (int idim = 0; idim <nsd_; ++idim)
        {
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = numdofpernode_*ui + idim;
            //const int fui   = 4*ui;
            //const int fuip  = fui+1;
            //const int fuipp = fui+2;
            double v = timefacfac*densaf_*funct_(ui)*vdiv_;
            if (physicaltype_ == INPAR::FLUID::loma) v -= timefacfac*densaf_*scaconvfacaf_*conv_scaaf_;
            // only with linear density-concentration correlation
            else if(physicaltype_ == INPAR::FLUID::varying_density) v += timefacfac*conv_scaaf_;
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = numdofpernode_*vi + idim;
              //const int fvi   = 4*vi;
              //const int fvip  = fvi+1;
              //const int fvipp = fvi+2;
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
              //estif(fvi  , fui  ) += v2;
              //estif(fvip , fuip ) += v2;
              //estif(fvipp, fuipp) += v2;
            }
          }

          if (is_newton_)
          {
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = numdofpernode_*vi + idim;
              //const int fvi   = 4*vi;
              //const int fvip  = fvi+1;
              //const int fvipp = fvi+2;
              const double v_idim = timefacfac*densaf_*velint_(idim)*funct_(vi);
              //const double v0 = timefacfac*densaf_*velint_(0)*funct_(vi);
              //const double v1 = timefacfac*densaf_*velint_(1)*funct_(vi);
              //const double v2 = timefacfac*densaf_*velint_(2)*funct_(vi);
              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui   = numdofpernode_*ui;
                //const int fui   = 4*ui;
                //const int fuip  = fui+1;
                //const int fuipp = fui+2;
                /*  convection, reactive part (conservative addition) */
                /*
                /                              \
                |  n+1  /               \      |
                | u    | rho*nabla o Du | , v  |
                |  (i)  \              /       |
                \                             /
                */
                for(int jdim=0; jdim<nsd_;++jdim)
                estif(fvi,  fui+jdim  ) += v_idim*derxy_(jdim, ui) ;

                //estif(fvi,  fui  ) += v0*derxy_(0, ui) ;
                //estif(fvi,  fuip ) += v0*derxy_(1, ui) ;
                //estif(fvi,  fuipp) += v0*derxy_(2, ui) ;
                //estif(fvip, fui  ) += v1*derxy_(0, ui) ;
                //estif(fvip, fuip ) += v1*derxy_(1, ui) ;
                //estif(fvip, fuipp) += v1*derxy_(2, ui) ;
                //estif(fvipp,fui  ) += v2*derxy_(0, ui) ;
                //estif(fvipp,fuip ) += v2*derxy_(1, ui) ;
                //estif(fvipp,fuipp) += v2*derxy_(2, ui) ;
              }
            }

            if (physicaltype_ == INPAR::FLUID::loma)
            {
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi   = numdofpernode_*vi + idim;
                //const int fvi   = 4*vi;
                //const int fvip  = fvi+1;
                //const int fvipp = fvi+2;
                const double v_idim = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(idim)*velint_(idim)*funct_(vi);
                //const double v0 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(0)*velint_(0)*funct_(vi);
                //const double v1 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(1)*velint_(1)*funct_(vi);
                //const double v2 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(2)*velint_(2)*funct_(vi);
                for (int ui=0; ui<nen_; ++ui)
                {
                  const int fui   = numdofpernode_*ui;
                  //const int fui   = 4*ui;
                  //const int fuip  = fui+1;
                  //const int fuipp = fui+2;
                  /*  convection, reactive part (conservative addition) */
                  /*
                  /                           \
                  |  n+1  /             \      |
                  | u    | Du*nabla rho | , v  |
                  |  (i)  \            /       |
                  \                           /
                  */
                  for(int jdim=0;jdim<nsd_;++jdim)
                    estif(fvi,  fui +jdim  ) += v_idim*funct_(ui) ;

                  //estif(fvi,  fui  ) += v0*funct_(ui) ;
                  //estif(fvi,  fuip ) += v0*funct_(ui) ;
                  //estif(fvi,  fuipp) += v0*funct_(ui) ;
                  //estif(fvip, fui  ) += v1*funct_(ui) ;
                  //estif(fvip, fuip ) += v1*funct_(ui) ;
                  //estif(fvip, fuipp) += v1*funct_(ui) ;
                  //estif(fvipp,fui  ) += v2*funct_(ui) ;
                  //estif(fvipp,fuip ) += v2*funct_(ui) ;
                  //estif(fvipp,fuipp) += v2*funct_(ui) ;
                }
              }
            }
            if (physicaltype_ == INPAR::FLUID::varying_density)
            {
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi   = numdofpernode_*vi + idim;
                //const int fvi   = 4*vi;
                //const int fvip  = fvi+1;
                //const int fvipp = fvi+2;
                const double v_idim = +timefacfac*grad_scaaf_(idim)*velint_(idim)*funct_(vi);
                //const double v0 = +timefacfac*grad_scaaf_(0)*velint_(0)*funct_(vi);
                //const double v1 = +timefacfac*grad_scaaf_(1)*velint_(1)*funct_(vi);
                //const double v2 = +timefacfac*grad_scaaf_(2)*velint_(2)*funct_(vi);
                for (int ui=0; ui<nen_; ++ui)
                {
                  const int fui   = numdofpernode_*ui;
                  //const int fui   = 4*ui;
                  //const int fuip  = fui+1;
                  //const int fuipp = fui+2;
                  /*  convection, reactive part (conservative addition) */
                  /*
                  /                           \
                  |  n+1  /             \      |
                  | u    | Du*nabla rho | , v  |
                  |  (i)  \            /       |
                  \                           /
                  */
                  for(int jdim=0;jdim<nsd_;++jdim)
                    estif(fvi,  fui+jdim  ) += v_idim*funct_(ui) ;

                  //estif(fvi,  fui  ) += v0*funct_(ui) ;
                  //estif(fvi,  fuip ) += v0*funct_(ui) ;
                  //estif(fvi,  fuipp) += v0*funct_(ui) ;
                  //estif(fvip, fui  ) += v1*funct_(ui) ;
                  //estif(fvip, fuip ) += v1*funct_(ui) ;
                  //estif(fvip, fuipp) += v1*funct_(ui) ;
                  //estif(fvipp,fui  ) += v2*funct_(ui) ;
                  //estif(fvipp,fuip ) += v2*funct_(ui) ;
                  //estif(fvipp,fuipp) += v2*funct_(ui) ;
                }
              }
            }
          }

          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = numdofpernode_*vi +idim;
            //const int fvi   = 4*vi;
            /* convection (conservative addition) on right-hand side */
            double v = -rhsfac*densaf_*funct_(vi)*vdiv_;
            eforce(fvi    ) += v*velint_(idim) ;
            //eforce(fvi    ) += v*velint_(0) ;
            //eforce(fvi + 1) += v*velint_(1) ;
            //eforce(fvi + 2) += v*velint_(2) ;
          }

          if (physicaltype_ == INPAR::FLUID::loma)
          {
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = numdofpernode_*vi+ idim;
              //const int fvi   = 4*vi;
              /* convection (conservative addition) on rhs for low-Mach-number flow */
              double v = rhsfac*densaf_*scaconvfacaf_*conv_scaaf_*funct_(vi);
              eforce(fvi    ) += v*velint_(idim) ;
              //eforce(fvi    ) += v*velint_(0) ;
              //eforce(fvi + 1) += v*velint_(1) ;
              //eforce(fvi + 2) += v*velint_(2) ;
            }
          }
          if (physicaltype_ == INPAR::FLUID::varying_density)
          {
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = numdofpernode_*vi + idim;
              //const int fvi   = 4*vi;
              /* convection (conservative addition) on rhs for low-Mach-number flow */
              double v = -rhsfac*conv_scaaf_*funct_(vi);
              eforce(fvi    ) += v*velint_(idim) ;
              //eforce(fvi    ) += v*velint_(0) ;
              //eforce(fvi + 1) += v*velint_(1) ;
              //eforce(fvi + 2) += v*velint_(2) ;
            }
          }
        }  // end for(idim)
      }  // end if(is_conservative)

      //----------------------------------------------------------------------
      // computation of viscosity term including right-hand-side contribution
      //----------------------------------------------------------------------
#if 0
      const double visceff_timefacfac = visceff_*timefacfac;
      for (int idim = 0; idim < nsd_; ++idim)
      {
        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = numdofpernode_*ui;
          //const int fui   = 4*ui;
          //const int fuip  = fui+1;
          //const int fuipp = fui+2;
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = numdofpernode_*vi;
            //const int fvi   = 4*vi;
            //const int fvip  = fvi+1;
            //const int fvipp = fvi+2;

            //const double derxy_0ui_0vi = derxy_(0, ui)*derxy_(0, vi);
            //const double derxy_1ui_1vi = derxy_(1, ui)*derxy_(1, vi);
            //const double derxy_2ui_2vi = derxy_(2, ui)*derxy_(2, vi);

            for(int jdim=0; jdim<nsd_; ++jdim)
            {
              const double derxy_0ui_0vi = derxy_(jdim, ui)*derxy_(jdim, vi);
              //const double derxy_0ui_0vi = derxy_(0, ui)*derxy_(0, vi);
              //const double derxy_1ui_1vi = derxy_(1, ui)*derxy_(1, vi);
              //const double derxy_2ui_2vi = derxy_(2, ui)*derxy_(2, vi);
              /* viscosity term */
              /*
                    /                          \
                    |       /  \         / \   |
              2 mu  |  eps | Du | , eps | v |  |
                    |       \  /         \ /   |
                    \                          /
              */

              estif(fvi +idim, fui+jdim)     += visceff_timefacfac*derxy_(idim, ui)*derxy_(jdim, vi);
              estif(fvi + idim, fui+idim)     += visceff_timefacfac*derxy_0ui_0vi;
            }

            /*estif(fvi, fui)     += visceff_timefacfac*(2.0*derxy_0ui_0vi
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
                                                       2.0*derxy_2ui_2vi) ; */
          }
        }


        //computation of right-hand-side viscosity term
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = numdofpernode_*vi;
          //const int fvi = 4*vi;
          const double v = -visceff_*rhsfac;
          /* viscosity term on right-hand side */
          for (int jdim = 0; jdim < nsd_; ++jdim)
          {
             eforce(fvi + idim) += v*derxy_(jdim,vi)*(vderxy_(jdim,idim)+vderxy_(idim,jdim));
          }

          /*eforce(fvi)     += v*(2.0*derxy_(0, vi)*vderxy_(0, 0)
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
                                2.0*derxy_(2, vi)*vderxy_(2, 2)) ;*/
        }
      } // end for(idim) viscosity term

#endif
      const double visceff_timefacfac = visceff_*timefacfac;

      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = numdofpernode_*vi;
          //const int fvi   = 4*vi;
          //const int fvip  = fvi+1;
          //const int fvipp = fvi+2;
          
          //const double derxy_0ui_0vi = derxy_(0, ui)*derxy_(0, vi);
          //const double derxy_1ui_1vi = derxy_(1, ui)*derxy_(1, vi);
          //const double derxy_2ui_2vi = derxy_(2, ui)*derxy_(2, vi);

          for(int jdim=0; jdim<nsd_; ++jdim)
          {
            const double visceff_timefacfac_derxy_jdim_vi=visceff_timefacfac*derxy_(jdim, vi);
          
            for (int ui=0; ui<nen_; ++ui)
            {

              const int fui   = numdofpernode_*ui;
              //const int fui   = 4*ui;
              //const int fuip  = fui+1;
              //const int fuipp = fui+2;

              for (int idim = 0; idim < nsd_; ++idim)
              {
                const int fvi_p_idim=fvi+idim;

                //const double derxy_0ui_0vi = derxy_(0, ui)*derxy_(0, vi);
                //const double derxy_1ui_1vi = derxy_(1, ui)*derxy_(1, vi);
                //const double derxy_2ui_2vi = derxy_(2, ui)*derxy_(2, vi);
                /* viscosity term */
                /*
                    /                          \
                    |       /  \         / \   |
              2 mu  |  eps | Du | , eps | v |  |
                    |       \  /         \ /   |
                    \                          /
                */
                
                estif(fvi_p_idim, fui+jdim) += visceff_timefacfac_derxy_jdim_vi*derxy_(idim,ui);
                estif(fvi_p_idim, fui+idim) += visceff_timefacfac_derxy_jdim_vi*derxy_(jdim,ui);
              }
              
              /*estif(fvi, fui)     += visceff_timefacfac*(2.0*derxy_0ui_0vi
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
                                                       2.0*derxy_2ui_2vi) ; */
            }
          }
        }

        for (int idim = 0; idim < nsd_; ++idim)
        {
          //computation of right-hand-side viscosity term
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi = numdofpernode_*vi;
            //const int fvi = 4*vi;
            const double v = -visceff_*rhsfac;
            /* viscosity term on right-hand side */
            for (int jdim = 0; jdim < nsd_; ++jdim)
            {
              eforce(fvi + idim) += v*derxy_(jdim,vi)*(vderxy_(jdim,idim)+vderxy_(idim,jdim));
            }

          /*eforce(fvi)     += v*(2.0*derxy_(0, vi)*vderxy_(0, 0)
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
                                2.0*derxy_(2, vi)*vderxy_(2, 2)) ;*/
          }
        }  // end for(idim) viscosity term
      }

      //----------------------------------------------------------------------
      // computation of pressure term including right-hand-side contribution
      //----------------------------------------------------------------------
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fuippp = numdofpernode_*ui+nsd_;
        //const int fuippp = 4*ui+3;
        const double v = -timefacfac*funct_(ui);
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = numdofpernode_*vi;
            //const int fvi   = 4*vi;
            /* pressure term */
            /*
            /                  \
            |                  |
            |  Dp , nabla o v  |
            |                  |
            \                  /
            */
          for (int idim = 0; idim <nsd_; ++idim)
          {
            estif(fvi + idim, fuippp) += v*derxy_(idim, vi) ;
            //estif(fvi,     fuippp) += v*derxy_(0, vi) ;
            //estif(fvi + 1, fuippp) += v*derxy_(1, vi) ;
            //estif(fvi + 2, fuippp) += v*derxy_(2, vi) ;
          }
        }
      }

      {
        const double v = press*rhsfac;

        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = numdofpernode_*vi;
          //const int fvi = 4*vi;
          /* pressure term on right-hand side */
          for (int idim = 0; idim <nsd_; ++idim)
          {
            eforce(fvi+idim) += v*derxy_(idim, vi) ;
            //eforce(fvi    ) += v*derxy_(0, vi) ;
            //eforce(fvi + 1) += v*derxy_(1, vi) ;
            //eforce(fvi + 2) += v*derxy_(2, vi) ;
          }
        }  //end for(idim)
      }

      //----------------------------------------------------------------------
      // computation of continuity term including right-hand-side contribution
      //----------------------------------------------------------------------
      for (int idim = 0; idim <nsd_; ++idim)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvippp = numdofpernode_*vi+nsd_;
          //const int fvippp = 4*vi+3;
          const double v = timefacfac*funct_(vi);
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = numdofpernode_*ui +idim;
            //const int fui   = 4*ui;
            /* continuity term */
            /*
              /                  \
              |                  |
              | nabla o Du  , q  |
              |                  |
              \                  /
            */
            estif(fvippp, fui)     += v*derxy_(idim, ui) ;
            //estif(fvippp, fui)     += v*derxy_(0, ui) ;
            //estif(fvippp, fui + 1) += v*derxy_(1, ui) ;
            //estif(fvippp, fui + 2) += v*derxy_(2, ui) ;
          }
        }
      }  // end for(idim)

      const double rhsfac_vdiv = -rhsfac * vdiv_;
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvippp = numdofpernode_*vi+nsd_;
        // continuity term on right-hand side
        eforce(fvippp) += rhsfac_vdiv*funct_(vi) ;
        //eforce(vi*4 + 3) += rhsfac_vdiv*funct_(vi) ;
      }


      //----------------------------------------------------------------------
      // computation of body-force term on right-hand side
      //----------------------------------------------------------------------
      for (int idim = 0; idim <nsd_; ++idim)
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = numdofpernode_*vi + idim;
          //const int fvi = 4*vi;
          const double v = fac_*funct_(vi);
          eforce(fvi    ) += v*rhsmom_(idim) ;
          //eforce(fvi    ) += v*rhsmom_(0) ;
          //eforce(fvi + 1) += v*rhsmom_(1) ;
          //eforce(fvi + 2) += v*rhsmom_(2) ;
        }
      }  // end for(idim)
      
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
          const int fui   = numdofpernode_*ui;
          //const int fui   = 4*ui;
          //const int fuip  = fui+1;
          //const int fuipp = fui+2;

          for (int idim = 0; idim <nsd_; ++idim)
          {
            const double v0 = v*derxy_(idim,ui);
            //const double v0 = v*derxy_(0,ui);
            //const double v1 = v*derxy_(1,ui);
            //const double v2 = v*derxy_(2,ui);
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = numdofpernode_*vi + idim;
              //const int fvi   = 4*vi;
              //const int fvip  = fvi+1;
              //const int fvipp = fvi+2;
              /* viscosity term - subtraction for low-Mach-number flow */
              /*
                    /                               \
                    |  1                      / \   |
             - 2 mu |  - (nabla o u) I , eps | v |  |
                    |  3                      \ /   |
                    \                               /
              */
              for(int jdim=0;jdim<nsd_;++jdim)
                estif(fvi,   fui +jdim  ) += v0*derxy_(jdim, vi) ;

              //estif(fvi,   fui  ) += v0*derxy_(0, vi) ;
              //estif(fvi,   fuip ) += v1*derxy_(0, vi) ;
              //estif(fvi,   fuipp) += v2*derxy_(0, vi) ;
              //estif(fvip,  fui  ) += v0*derxy_(1, vi) ;
              //estif(fvip,  fuip ) += v1*derxy_(1, vi) ;
              //estif(fvip,  fuipp) += v2*derxy_(1, vi) ;
              //estif(fvipp, fui  ) += v0*derxy_(2, vi) ;
              //estif(fvipp, fuip ) += v1*derxy_(2, vi) ;
              //estif(fvipp, fuipp) += v2*derxy_(2, vi) ;
            }
          }
        }  // end for(idim)
        
        if (is_newton_)
        {
          const double timefacfac_scaconvfacaf=timefacfac*scaconvfacaf_;

          LINALG::Matrix<nsd_,1> temp;

          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvippp= numdofpernode_*vi+nsd_;
            const double timefacfac_scaconvfacaf_funct_vi=timefacfac_scaconvfacaf*funct_(vi);
            
            for(int jdim=0;jdim<nsd_;++jdim)
            {
              temp(jdim)=timefacfac_scaconvfacaf_funct_vi*grad_scaaf_(jdim);
            }
            
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui=numdofpernode_*ui;

              for(int jdim=0;jdim<nsd_;++jdim)
              {
                /*
                  factor afgtd/am

                          /                    \                                      
                    1    |       /         \    |                                      
                   --- * |  q , | Du o grad | T |                                           
                    T    |       \         /    |
                          \                    /                                  
                */
                estif(fvippp,fui+jdim) -= temp(jdim)*funct_(ui);  
              }
            }
          }
        } // end if (is_newton_)
        
        const double v_div = (2.0/3.0)*visceff_*rhsfac*vdiv_ ;
        const double fac_rhscon = fac_*rhscon_;
        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi = numdofpernode_*vi;
          
          for(int idim=0;idim<nsd_;++idim)
          {
            //const int fvi = 4*vi;
            /* viscosity term on rhs - subtraction for low-Mach-number flow */
            eforce(fvi +idim) += derxy_(idim, vi)*v_div ;
            //eforce(fvi    ) += derxy_(0, vi)*v_div ;
            //eforce(fvi + 1) += derxy_(1, vi)*v_div ;
            //eforce(fvi + 2) += derxy_(2, vi)*v_div ;
          }

          /* additional rhs term of continuity equation */
          eforce(fvi + nsd_) += fac_rhscon*funct_(vi) ;
          //eforce(fvi + 3) += fac_rhscon*funct_(vi) ;
        }

      }  // end if(loma)

      //----------------------------------------------------------------------
      //  PROVIDE LINEARISATION OF MOMENTUM RESIDUAL WRT VELOCITIES FOR ALL              
      //                     STABILISATION PARTS

      //                                                             
      //           /       n+1        \       /                \  n+1 
      // rho*Du + |   rho*u   o nabla | Du + |   rho*Du o nabla | u   +
      //           \      (i)         /       \                /   (i)
      //
      //                       /  \ 
      //        + nabla o eps | Du |
      //                       \  / 

      LINALG::Matrix<nen_,nsd_*nsd_> lin_resM_Du;

      lin_resM_Du.Clear();

      const double timefac_densaf=timefac*densaf_;

      for (int ui=0; ui<nen_; ++ui)
      {  
        for (int idim = 0; idim <nsd_; ++idim)
        {
          lin_resM_Du(ui,idim*nsd_+idim)+=timefac_densaf*conv_c_(ui);
        }
      }

      if(!is_stationary_)
      {
        for (int ui=0; ui<nen_; ++ui)
        {  
          for (int idim = 0; idim <nsd_; ++idim)
          {
            lin_resM_Du(ui,idim*nsd_+idim)+=densam_*funct_(ui);
          }
        }
      }
      
      if (is_higher_order_ele_)
      {
        const double v = -2.0*visceff_*timefac;
        
        for (int ui=0; ui<nen_; ++ui)
        {  
          for (int idim = 0; idim <nsd_; ++idim)
          {
            for(int jdim=0;jdim<nsd_;++jdim)
            {
              lin_resM_Du(ui,idim*nsd_+jdim)+=v*viscs2_(jdim+nsd_*idim, ui);
            }
          }
        }
      }

      if (is_newton_)
      {
        for (int idim = 0; idim <nsd_; ++idim)
        {
          for(int jdim=0;jdim<nsd_;++jdim)
          {   
            const double temp=timefac*densaf_*vderxy_(jdim,idim);
            
            for (int ui=0; ui<nen_; ++ui)
            {  
              lin_resM_Du(ui,idim*nsd_+jdim)+=temp*funct_(ui);
            }
          }
        }
      }

      //----------------------------------------------------------------------
      //                 PRESSURE STABILISATION PART

      // conservative, stabilization terms are neglected (Hughes)

      if (pspg_ == INPAR::FLUID::pstab_use_pspg)
      {
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
              |  /       n+1        \               |
              | |   rho*u   o nabla | Du , nabla q  |
              |  \      (i)         /               |
               \                                   /
        */
        /* pressure stabilisation: convection, reactive part if Newton */
        /*
                 /                                   \
                |  /                 \  n+1           |
                | |   rho*Du o nabla | u     , grad q |
                |  \                /   (i)           |
                 \                                   /
        */
        /* pressure stabilisation: viscosity (-L_visc_u) */
        /*
                 /                              \
                |               /  \             |
                |  nabla o eps | Du | , nabla q  |
                |               \  /             |
                 \                              /
        */

        if (is_higher_order_ele_ || is_newton_)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvippp = numdofpernode_*vi+nsd_;

            for(int jdim=0;jdim<nsd_;++jdim)
            {
              const double temp_vi_jdim=derxy_(jdim, vi)*tau_Mp;
              
              for(int idim=0;idim<nsd_;++idim)
              {
                const int nsd_idim=nsd_*idim;
                
                for (int ui=0; ui<nen_; ++ui)
                {
                  const int fui_p_idim   = numdofpernode_*ui + idim;
                  
                  estif(fvippp,fui_p_idim) += lin_resM_Du(ui,nsd_idim+jdim)*temp_vi_jdim;
                } // jdim
              } // vi
            } // ui
          } //idim
        } // end if (is_higher_order_ele_) or (newton_)
        else
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvippp = numdofpernode_*vi+nsd_;
            
            for(int idim=0;idim<nsd_;++idim)
            {
              const int nsd_idim=nsd_*idim;
              
              const double temp_vi_idim=derxy_(idim, vi)*tau_Mp;
              
              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui_p_idim   = numdofpernode_*ui + idim;
                
                estif(fvippp,fui_p_idim) += lin_resM_Du(ui,nsd_idim+idim)*temp_vi_idim;
              } // vi
            } // ui
          } //idim
        } // end if not (is_higher_order_ele_) nor (newton_)

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fuippp = numdofpernode_*ui+nsd_;
          
          for (int idim = 0; idim <nsd_; ++idim)
          {
            const double timetauMp_derxy_idim_ui=timetauMp*derxy_(idim, ui);
            
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
              estif(numdofpernode_*vi+nsd_,fuippp)
                += timetauMp_derxy_idim_ui*derxy_(idim, vi);
            } // vi
          } // end for(idim)
        }  // ui
              
        for (int idim = 0; idim <nsd_; ++idim)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            // pressure stabilisation
            eforce(numdofpernode_*vi+nsd_) -= tau_Mp*momres_old_(idim)*derxy_(idim, vi);
          }
        } // end for(idim)
      }  // end PSPG stabilization

      //----------------------------------------------------------------------
      //                     SUPG STABILISATION PART

      if(supg_ == INPAR::FLUID::convective_stab_supg)
      {
        LINALG::Matrix<nen_,1> supg_test;
              
        for (int vi=0; vi<nen_; ++vi)
        {
          supg_test(vi)=(conv_c_(vi)+sgconv_c_(vi))*tau_M;
        }

        /* supg stabilisation: inertia if not stationary */
        /*
                   /                                \
                  |           /      n+1        \    |
                  |  rho*Du , | rho*u   o nabla | v  |
                  |           \      (i)        /    |
                   \                                /
        */
        /* supg stabilisation: convective part ( L_conv_u) , convective term */
        /*
               /                                                     \
              |    /       n+1         \      /      n+1        \     |
              |   |   rho*u    o nabla | Du , | rho*u    o nabla | v  |
              |    \       (i)         /      \      (i)        /     |
               \                                                     /
        */
        /* supg stabilisation: viscous part  (-L_visc_u) if is_higher_order_ele_ */
        /*
                 /                                              \
                |               /  \    /       n+1        \     |
                |  nabla o eps | Du |, |   rho*u    o nabla | v  |
                |               \  /    \       (i)        /     |
                 \                                              /
        */
        /* supg stabilisation: convective part ( L_conv_u) , reactive term if Newton */
        /*
               /                                                     \
              |    /       n+1         \      /      n+1        \     |
              |   |   rho*u    o nabla | Du , | rho*u    o nabla | v  |
              |    \       (i)         /      \      (i)        /     |
               \                                                     /
        */
        if (is_higher_order_ele_ || is_newton_)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            for(int jdim=0;jdim<nsd_;++jdim)
            {
              const int fvi_p_jdim = numdofpernode_*vi+jdim;

              for(int idim=0;idim<nsd_;++idim)
              {
                const int nsd_idim=nsd_*idim;
                
                for (int ui=0; ui<nen_; ++ui)
                {
                  const int fui_p_idim   = numdofpernode_*ui + idim;
                  
                  estif(fvi_p_jdim,fui_p_idim) += lin_resM_Du(ui,nsd_idim+jdim)*supg_test(vi);
                } // jdim
              } // vi
            } // ui
          } //idim
        } // end if (is_higher_order_ele_) or (newton_)
        else
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            for(int idim=0;idim<nsd_;++idim)
            {
              const int fvi_p_idim = numdofpernode_*vi+idim;

              const int nsd_idim=nsd_*idim;
                
              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui_p_idim   = numdofpernode_*ui + idim;
                  
                estif(fvi_p_idim,fui_p_idim) += lin_resM_Du(ui,nsd_idim+idim)*supg_test(vi);
              } // ui
            } //idim
          } // vi
        } // end if not (is_higher_order_ele_) nor (newton_)

        /* supg stabilisation: pressure part  ( L_pres_p) */
        /*
                 /                                    \
                |              /       n+1       \     |
                |  nabla Dp , |   rho*u   o nabla | v  |
                |              \       (i)       /     |
                 \                                    /
        */
        for (int vi=0; vi<nen_; ++vi)
        {
          const double v = timetauM*densaf_*supg_test(vi);

          for (int idim = 0; idim <nsd_; ++idim)
          {
            const int fvi   = numdofpernode_*vi + idim;
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fuippp = numdofpernode_*ui + nsd_;

              estif(fvi,fuippp) += v*derxy_(idim, ui);
            }
          }
        }  // end for(idim)

        /* supg stabilisation: inertia, linearisation of testfunction  */
        /*
                     /                                       \
                    |         n+1      /                \     |
                    |    rho*u      ,  | rho*Du o nabla | v   |
                    |         (i)      \                /     |
                     \                                       /
        */
        /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
        /*
                     /                                                       \
                    |    /       n+1        \   n+1    /                \     |
                    |   |   rho*u    o nabla | u    ,  | rho*Du o nabla | v   |
                    |    \       (i)        /   (i)    \                /     |
                     \                                                       /
        */
        /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
        /*
                   /                                     \
                  |         n+1    /                \     |
                  |  nabla p    , |   rho*Du o nabla | v  |
                  |         (i)    \                /     |
                   \                                     /
        */
        /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
        /*
                   /                                               \
                  |               / n+1 \    /                \     |
                  |  nabla o eps | u     |, |  rho*Du o nabla | v   |
                  |               \ (i) /    \                /     |
                   \                                               /
        */
        /* supg stabilisation: bodyforce part, linearisation of test function */
        /*
                   /                                      \
                  |                 /                \     |
                  |  rho*rhsint   , |  rho*Du o nabla | v  |
                  |                 \                /     |
                   \                                      /
        */
        if (is_newton_)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            for(int jdim=0;jdim<nsd_;++jdim)
            {
              const int fvi_p_jdim = numdofpernode_*vi+jdim;

              for(int idim=0;idim<nsd_;++idim)
              {
                for (int ui=0; ui<nen_; ++ui)
                {
                  const int fui_p_idim   = numdofpernode_*ui + idim;
                  
                  estif(fvi_p_jdim,fui_p_idim) += timetauM*densaf_*momres_old_(jdim)*derxy_(idim,vi)*funct_(ui);
                } // jdim
              } // vi
            } // ui
          } //idim
        }

        // NOTE: Here we have a difference to the previous version of this
        // element!  Before we did not care for the mesh velocity in this
        // term. This seems unreasonable and wrong.

        for (int idim = 0; idim <nsd_; ++idim)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi = numdofpernode_*vi +idim;
            //const int fvi = 4*vi;
            // supg stabilisation
            const double v = -tau_M*densaf_*(conv_c_(vi)+sgconv_c_(vi));
            eforce(fvi)     += v*momres_old_(idim) ;
          }
        }  // end for(idim)
      }  // end SUPG stabilization

      //----------------------------------------------------------------------
      //                       STABILISATION, VISCOUS PART
      
      if (is_higher_order_ele_ && (vstab_ != INPAR::FLUID::viscous_stab_none))
      {
        const double two_visc_tauMp = vstabfac*2.0*visc_*tau_Mp;
        // viscous stabilization either on left hand side or on right hand side
        if (vstab_ == INPAR::FLUID::viscous_stab_gls || vstab_ == INPAR::FLUID::viscous_stab_usfem)
        {
          const double two_visc_timetauMp   = vstabfac*2.0*visc_*timetauMp;
          
          /* viscous stabilisation, inertia part if not stationary */
          /*
                    /                          \
                    |                          |
                +/- |    rho*Du , div eps (v)  |
                    |                          |
                    \                          /
          */
          /* viscous stabilisation, convective part, convective type */
          /*
                  /                                        \
                  |  /       n+1       \                   |
              +/- | |   rho*u   o nabla | Du , div eps (v) |
                  |  \       (i)       /                   |
                  \                                        /
          */
          /* viscous stabilisation, reactive part of convection */
          /*
                    /                                         \
                    |  /                \   n+1               |
                +/- | |   rho*Du o nabla | u    , div eps (v) |
                    |  \                /   (i)               |
                    \                                         /
          */
          /* viscous stabilisation, viscous part (-L_visc_u) */
          /*
                  /                                 \
                  |               /  \                |
             -/+  |  nabla o eps | Du | , div eps (v) |
                  |               \  /                |
                  \                                 /
          */
          for (int vi=0; vi<nen_; ++vi)
          {
            for(int jdim=0;jdim<nsd_;++jdim)
            {
              const int fvi_p_jdim = numdofpernode_*vi+jdim;
                
              for(int idim=0;idim<nsd_;++idim)
              {
                const int nsd_idim=nsd_*idim;

                for (int ui=0; ui<nen_; ++ui)
                {
                  const int fui_p_idim   = numdofpernode_*ui + idim;
                  estif(fvi_p_jdim,fui_p_idim) += two_visc_timetauMp*lin_resM_Du(ui,nsd_idim+jdim)*viscs2_(jdim+nsd_idim,vi);
                } // jdim
              } // vi
            } // ui
          } //idim
         

          /* viscous stabilisation, pressure part ( L_pres_p) */
          /*
                   /                        \
                  |                          |
             +/-  |  nabla Dp , div eps (v)  |
                  |                          |
                   \                        /
          */
          for (int idim=0;idim<nsd_; ++idim)
          {
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui = ui*numdofpernode_ + nsd_;
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi = vi*numdofpernode_ + idim;

                for(int jdim=0;jdim<nsd_;++jdim)
                {
                  estif(fvi, fui) += two_visc_timetauMp*(derxy_(jdim, ui))*viscs2_(jdim+(idim*nsd_), vi);
                }
              }
            }
          } // end for(idim)


          for (int idim =0;idim<nsd_;++idim)
          {
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi=vi*numdofpernode_ +idim;
              /* viscous stabilisation */
              for (int jdim=0;jdim<nsd_;++jdim)
              eforce(fvi) -= two_visc_tauMp*momres_old_(jdim)*viscs2_(jdim+(idim*nsd_), vi);
            }
          } // end for(idim)

        } // end if viscous stabilization on left hand side
      } //end if(is_higherOrder)




#if 0
        {
          if (is_conservative_ == false)
          {
            if(is_genalpha_)
            {
              const double timetauM_densam_densaf=timetauM*densam_*densaf_;

              LINALG::Matrix<nsd_,1> timetauM_densam_densaf_accint;
              for (int idim = 0; idim <nsd_; ++idim)
              {
                timetauM_densam_densaf_accint(idim)=timetauM_densam_densaf*accint_(idim);
              }

              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui   = numdofpernode_*ui;

                for (int idim = 0; idim <nsd_; ++idim)
                {
                  //const int fui   = 4*ui;
                  //const int fuip  = fui+1;
                  //const int fuipp = fui+2;
                  const double v = timetauM_densam_densaf_accint(idim)*funct_(ui);
                  //const double v0 = v*velint_(0);
                  //const double v1 = v*velint_(1);
                  //const double v2 = v*velint_(2);
                  for (int vi=0; vi<nen_; ++vi)
                  {
                    const int fvi_p_idim = numdofpernode_*vi + idim;
                    //const int fvi   = 4*vi;
                    //const int fvip  = fvi+1;
                    //const int fvipp = fvi+2;
                    /* supg stabilisation: inertia, linearisation of testfunction  */
                    /*
                    /                                         \
                    |         n+1      /                \     |
                    |    rho*u      ,  | rho*Du o nabla | v   |
                    |         (i)      \                /     |
                    \                                         /

                    */
                    for (int jdim = 0; jdim < nsd_; ++jdim)
                    {
                      estif(fvi_p_idim,fui+jdim) += v*derxy_(jdim,vi);
                    }
                    //estif(fvi,  fui)     += v0*derxy_(0, vi) ;
                    //estif(fvip,  fui)    += v1*derxy_(0, vi) ;
                    //estif(fvipp, fui)    += v2*derxy_(0, vi) ;
                    
                    //estif(fvi,   fuip)   += v0*derxy_(1, vi) ;
                    //estif(fvip,  fuip)   += v1*derxy_(1, vi) ;
                    //estif(fvipp, fuip)   += v2*derxy_(1, vi) ;
                    
                    //estif(fvi,   fuipp)  += v0*derxy_(2, vi) ;
                    //estif(fvip,  fuipp)  += v1*derxy_(2, vi) ;
                    //estif(fvipp, fuipp)  += v2*derxy_(2, vi) ;
                  }
                }
              }  // end for(idim)
            }
            else
            {
              for (int idim = 0; idim <nsd_; ++idim)
              {
                for (int ui=0; ui<nen_; ++ui)
                {
                  const int fui   = numdofpernode_*ui;
                  //const int fui   = 4*ui;
                  //const int fuip  = fui+1;
                  //const int fuipp = fui+2;
                  const double v = tau_M*densam_*densaf_*funct_(ui);
                  const double v_idim = v*velint_(idim);
                  //const double v0 = v*velint_(0);
                  //const double v1 = v*velint_(1);
                  //const double v2 = v*velint_(2);
                  for (int vi=0; vi<nen_; ++vi)
                  {
                    const int fvi   = numdofpernode_*vi + idim;
                    //const int fvi   = 4*vi;
                    //const int fvip  = fvi+1;
                    //const int fvipp = fvi+2;
                    /* supg stabilisation: inertia, linearisation of testfunction  */
                    /*
                    /                                         \
                    |         n+1      /                \     |
                    |    rho*u      ,  | rho*Du o nabla | v   |
                    |         (i)      \                /     |
                    \                                         /

                    */
                    for (int jdim = 0; jdim < nsd_; ++jdim)
                    {
                      estif(fvi,  fui + jdim)     += v_idim*derxy_(jdim, vi) ;
                    }
                    //estif(fvi,  fui)     += v0*derxy_(0, vi) ;
                    //estif(fvip,  fui)    += v1*derxy_(0, vi) ;
                    //estif(fvipp, fui)    += v2*derxy_(0, vi) ;
                    
                    //estif(fvi,   fuip)   += v0*derxy_(1, vi) ;
                    //estif(fvip,  fuip)   += v1*derxy_(1, vi) ;
                    //estif(fvipp, fuip)   += v2*derxy_(1, vi) ;
                    
                    //estif(fvi,   fuipp)  += v0*derxy_(2, vi) ;
                    //estif(fvip,  fuipp)  += v1*derxy_(2, vi) ;
                    //estif(fvipp, fuipp)  += v2*derxy_(2, vi) ;
                  }
                }
              }  // end for(idim)
            } // end if not genalpha
          }  // end if(is_conservative)
          
          {

            const double timetauM_densaf_densaf=timetauM*densaf_*densaf_;

            LINALG::Matrix<nsd_,nsd_> scaled_vderxy;

            LINALG::Matrix<nsd_,1> convvelint_vderxy(true);
            LINALG::Matrix<nsd_,1> temp;

            convvelint_vderxy.Multiply(vderxy_,convvelint_);

            //const double v0 = convvelint_(0)*vderxy_(0, 0) + convvelint_(1)*vderxy_(0, 1) + convvelint_(2)*vderxy_(0, 2);
            //const double v1 = convvelint_(0)*vderxy_(1, 0) + convvelint_(1)*vderxy_(1, 1) + convvelint_(2)*vderxy_(1, 2);
            //const double v2 = convvelint_(0)*vderxy_(2, 0) + convvelint_(1)*vderxy_(2, 1) + convvelint_(2)*vderxy_(2, 2);
            
            for (int ui=0; ui<nen_; ++ui)
            {

              const int fui   = numdofpernode_*ui;
                //const int fui   = 4*ui;
                //const int fuip  = fui+1;
                //const int fuipp = fui+2;

              const double timetauM_densaf_densaf_funct_ui=funct_(ui)*timetauM_densaf_densaf;

              for (int idim = 0; idim <nsd_; ++idim)
              {    
                temp(idim)=timetauM_densaf_densaf_funct_ui*convvelint_vderxy(idim);

                for(int jdim=0;jdim<nsd_;++jdim)
                {
                  scaled_vderxy(idim,jdim)=timetauM_densaf_densaf_funct_ui*vderxy_(idim,jdim);
                }
              }

              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi=numdofpernode_*vi;
                
                for (int idim = 0; idim <nsd_; ++idim)
                {
                  const int fvi_p_idim   = fvi +idim;

                  //const int fvi   = 4*vi;
                  //const int fvip  = fvi+1;
                  //const int fvipp = fvi+2;
                  /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
                  /*
                    /                                                         \
                    |    /       n+1        \   n+1    /                \     |
                    |   |   rho*u    o nabla | u    ,  | rho*Du o nabla | v   |
                    |    \       (i)        /   (i)    \                /     |
                    \                                                        /

                  */
                  for(int jdim=0;jdim<nsd_;++jdim)
                  {
                    estif(fvi_p_idim, fui+jdim) += temp(idim)*derxy_(jdim, vi);
                  }
                  //estif(fvi, fui)     += (conv_c_(vi,0)*vderxy_(0, 0) + v0*derxy_(0, vi))*v;
                  //estif(fvip, fui)    += (conv_c_(vi,0)*vderxy_(1, 0) + v1*derxy_(0, vi))*v;
                  //estif(fvipp, fui)   += (conv_c_(vi,0)*vderxy_(2, 0) + v2*derxy_(0, vi))*v;

                  //estif(fvi, fuip)    += (conv_c_(vi,0)*vderxy_(0, 1) + v0*derxy_(1, vi))*v;
                  //estif(fvip, fuip)   += (conv_c_(vi,0)*vderxy_(1, 1) + v1*derxy_(1, vi))*v;
                  //estif(fvipp, fuip)  += (conv_c_(vi,0)*vderxy_(2, 1) + v2*derxy_(1, vi))*v;

                  //estif(fvi, fuipp)   += (conv_c_(vi,0)*vderxy_(0, 2) + v0*derxy_(2, vi))*v;
                  //estif(fvip, fuipp)  += (conv_c_(vi,0)*vderxy_(1, 2) + v1*derxy_(2, vi))*v;
                  //estif(fvipp, fuipp) += (conv_c_(vi,0)*vderxy_(2, 2) + v2*derxy_(2, vi))*v;
                }
              }
            }  // end for(idim)
          } // end

          for (int idim = 0; idim <nsd_; ++idim)
          {
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = numdofpernode_*ui + idim;
              //const int fui   = 4*ui;
              //const int fuip  = fui+1;
              //const int fuipp = fui+2;
              const double v = timetauM*densaf_*funct_(ui);
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi   = numdofpernode_*vi;
                //const int fvi   = 4*vi;
                //const int fvip  = fvi+1;
                //const int fvipp = fvi+2;
                /* supg stabilisation: pressure part, linearisation of test function  ( L_pres_p) */
                /*
                  /                                       \
                  |         n+1    /                \     |
                  |  nabla p    , |   rho*Du o nabla | v  |
                  |         (i)    \                /     |
                  \                                      /
                */
                for (int jdim=0; jdim<nsd_;++jdim)
                  estif(fvi+jdim, fui)     += v*gradp_(jdim,0)*derxy_(idim, vi) ;

                //estif(fvi, fui)     += v*gradp_(0,0)*derxy_(0, vi) ;
                //estif(fvip, fui)    += v*gradp_(1,0)*derxy_(0, vi) ;
                //estif(fvipp, fui)   += v*gradp_(2,0)*derxy_(0, vi) ;

                //estif(fvi, fuip)    += v*gradp_(0,0)*derxy_(1, vi) ;
                //estif(fvip, fuip)   += v*gradp_(1,0)*derxy_(1, vi) ;
                //estif(fvipp, fuip)  += v*gradp_(2,0)*derxy_(1, vi) ;

                //estif(fvi, fuipp)   += v*gradp_(0,0)*derxy_(2, vi) ;
                //estif(fvip, fuipp)  += v*gradp_(1,0)*derxy_(2, vi) ;
                //estif(fvipp, fuipp) += v*gradp_(2,0)*derxy_(2, vi) ;
              }
            }
          }  // end for(idim)

          if (is_higher_order_ele_)
          {
            const double scaling=-2.0*visceff_*timetauM*densaf_;


            for(int idim =0;idim<nsd_;++idim)
            {

              const double visc_scaling=scaling*visc_old_(idim);

              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui  = numdofpernode_*ui;

                const double v = visc_scaling*funct_(ui);

                for (int vi=0; vi<nen_; ++vi)
                {
                  const int fvi   = numdofpernode_*vi + idim;

                  /* supg stabilisation: viscous part, linearisation of test function  (-L_visc_u) */
                  /*
                  /                                                 \
                  |               / n+1 \    /                \     |
                  |  nabla o eps | u     |, |  rho*Du o nabla | v   |
                  |               \ (i) /    \                /     |
                  \                                                 /
                  */
                  for(int jdim=0;jdim<nsd_;++jdim)
                  {
                    estif(fvi, fui +jdim) += v*derxy_(jdim, vi) ;
                  }
                }
              }
            }
          }

          if(is_genalpha_)
          {
            for (int idim = 0; idim <nsd_; ++idim)
            {
              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui   = numdofpernode_*ui;
                //const int fui   = 4*ui;
                //const int fuip  = fui+1;
                //const int fuipp = fui+2;
                const double v = -timetauM*densaf_*funct_(ui);
                for (int vi=0; vi<nen_; ++vi)
                {
                    const int fvi   = numdofpernode_*vi + idim;
                    //const int fvi   = 4*vi;
                    //const int fvip  = fvi+1;
                    //const int fvipp = fvi+2;
                    const double v_idim = v*rhsmom_(idim);
                    //const double v0 = v*rhsmom_(0);
                    //const double v1 = v*rhsmom_(1);
                    //const double v2 = v*rhsmom_(2);

                    /* supg stabilisation: bodyforce part, linearisation of test function */

                    /*
                  /                                       \
                  |                 /                \     |
                  |  rho*rhsint   , |  rho*Du o nabla | v  |
                  |                 \                /     |
                  \                                        /

                    */
                    for(int jdim=0; jdim < nsd_; ++jdim)
                      estif(fvi, fui + jdim) += (v_idim*derxy_(jdim, vi)) ;

                    //estif(fvi    , fui)   += (v0*derxy_(0, vi)) ;
                    //estif(fvip   , fui)   += (v1*derxy_(0, vi)) ;
                    //estif(fvipp  , fui)   += (v2*derxy_(0, vi)) ;

                    //estif(fvi    , fuip)  += (v0*derxy_(1, vi)) ;
                    //estif(fvip   , fuip)  += (v1*derxy_(1, vi)) ;
                    //estif(fvipp  , fuip)  += (v2*derxy_(1, vi)) ;

                    //estif(fvi    , fuipp) += (v0*derxy_(2, vi)) ;
                    //estif(fvip   , fuipp) += (v1*derxy_(2, vi)) ;
                    //estif(fvipp  , fuipp) += (v2*derxy_(2, vi)) ;
                    
                  } // vi
              } // ui
            }  // end for(idim)
          }
          else
          {
            for (int idim = 0; idim <nsd_; ++idim)
            {
              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui   = numdofpernode_*ui;
                //const int fui   = 4*ui;
                
                  //const int fuip  = fui+1;
                  //const int fuipp = fui+2;
                  const double v = -tau_M*densaf_*funct_(ui);
                  for (int vi=0; vi<nen_; ++vi)
                  {
                    const int fvi   = numdofpernode_*vi + idim;
                    //const int fvi   = 4*vi;
                    //const int fvip  = fvi+1;
                    //const int fvipp = fvi+2;
                    const double v_idim = v*rhsmom_(idim);
                    //const double v0 = v*rhsmom_(0);
                    //const double v1 = v*rhsmom_(1);
                    //const double v2 = v*rhsmom_(2);

                    /* supg stabilisation: bodyforce part, linearisation of test function */

                    /*
                  /                                       \
                  |                 /                \     |
                  |  rho*rhsint   , |  rho*Du o nabla | v  |
                  |                 \                /     |
                  \                                        /

                    */
                    for(int jdim=0; jdim < nsd_; ++jdim)
                      estif(fvi, fui + jdim) += (v_idim*derxy_(jdim, vi)) ;

                    //estif(fvi    , fui)   += (v0*derxy_(0, vi)) ;
                    //estif(fvip   , fui)   += (v1*derxy_(0, vi)) ;
                    //estif(fvipp  , fui)   += (v2*derxy_(0, vi)) ;

                    //estif(fvi    , fuip)  += (v0*derxy_(1, vi)) ;
                    //estif(fvip   , fuip)  += (v1*derxy_(1, vi)) ;
                    //estif(fvipp  , fuip)  += (v2*derxy_(1, vi)) ;

                    //estif(fvi    , fuipp) += (v0*derxy_(2, vi)) ;
                    //estif(fvip   , fuipp) += (v1*derxy_(2, vi)) ;
                    //estif(fvipp  , fuipp) += (v2*derxy_(2, vi)) ;
                    
                  } // vi
              } // ui
            }  // end for(idim)
          }
        } // end if(newton)
#endif


#if 0
            const double four_visc2_timetauMp = vstabfac*4.0*visceff_*visc_*timetauMp;
            // viscous stabilization on left hand side
            if (is_stationary_ == false)
            {
              for (int idim=0;idim<nsd_; ++idim)
              {
              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui = ui * numdofpernode_;
                const double v = two_visc_tauMp*densam_*funct_(ui);
                for (int vi=0; vi<nen_; ++vi)
                {
                  const int fvi = vi * numdofpernode_ +idim;

                  for(int jdim=0;jdim<nsd_;++jdim)
                  {


                  estif(fvi, fui+jdim) += v*viscs2_(jdim+(idim*nsd_), vi) ;

                  //estif(vi*4,     ui*4    ) += v*viscs2_(0, vi) ;
                  //estif(vi*4 + 1, ui*4    ) += v*viscs2_(1, vi) ;
                  //estif(vi*4 + 2, ui*4    ) += v*viscs2_(2, vi) ;

                  //estif(vi*4,     ui*4 + 1) += v*viscs2_(1, vi) ;
                  //estif(vi*4 + 1, ui*4 + 1) += v*viscs2_(4, vi) ;
                  //estif(vi*4 + 2, ui*4 + 1) += v*viscs2_(5, vi) ;

                  //estif(vi*4,     ui*4 + 2) += v*viscs2_(2, vi) ;
                  //estif(vi*4 + 1, ui*4 + 2) += v*viscs2_(5, vi) ;
                  //estif(vi*4 + 2, ui*4 + 2) += v*viscs2_(8, vi) ;
                }
                }
              }
              } // end for(idim)
            }  // end if(stationary)

            // viscous stabilization on left hand side
            for (int idim=0;idim<nsd_; ++idim)
            {
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui = ui * numdofpernode_;
              const double v = two_visc_timetauMp*densaf_*conv_c_(ui);
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi = vi * numdofpernode_ +idim;
                for(int jdim=0;jdim<nsd_;++jdim)
                {
                /* viscous stabilisation, inertia part */



                  estif(fvi, fui + jdim   ) += v*viscs2_(jdim+(idim*nsd_), vi) ;

                //estif(vi*4,     ui*4    ) += v*viscs2_(0, vi) ;
                //estif(vi*4 + 1, ui*4    ) += v*viscs2_(1, vi) ;
                //estif(vi*4 + 2, ui*4    ) += v*viscs2_(2, vi) ;

                //estif(vi*4,     ui*4 + 1) += v*viscs2_(1, vi) ;
                //estif(vi*4 + 1, ui*4 + 1) += v*viscs2_(4, vi) ;
                //estif(vi*4 + 2, ui*4 + 1) += v*viscs2_(5, vi) ;

                //estif(vi*4,     ui*4 + 2) += v*viscs2_(2, vi) ;
                //estif(vi*4 + 1, ui*4 + 2) += v*viscs2_(5, vi) ;
                //estif(vi*4 + 2, ui*4 + 2) += v*viscs2_(8, vi) ;
              }
              } // end for(jdim)
            }
            } // end for(idim)


            for(int idim =0;idim<nsd_;++idim)
            {
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui=numdofpernode_*ui;
              for (int vi=0; vi<nen_; ++vi)
              {
                const int fvi = vi*numdofpernode_+idim;
                /* viscous stabilisation, viscous part (-L_visc_u) */
                /*
                  /                                 \
                  |               /  \                |
             -/+  |  nabla o eps | Du | , div eps (v) |
                  |               \  /                |
                  \                                 /
                */
                for (int jdim=0; jdim<nsd_; ++jdim)
                  for (int rdim=0;rdim<nsd_;++rdim)
                    estif(fvi, fui+jdim) -= four_visc2_timetauMp*viscs2_(rdim+(idim*nsd_),ui)*viscs2_(rdim+(jdim*nsd_),vi);

                //estif(vi*4,     ui*4    ) -= four_visc2_timetauMp*(viscs2_(0,ui)*viscs2_(0,vi)+viscs2_(1,ui)*viscs2_(1,vi)+viscs2_(2,ui)*viscs2_(2,vi)) ;
                //estif(vi*4 + 1, ui*4    ) -= four_visc2_timetauMp*(viscs2_(0,ui)*viscs2_(1,vi)+viscs2_(1,ui)*viscs2_(4,vi)+viscs2_(2,ui)*viscs2_(5,vi)) ;
                //estif(vi*4 + 2, ui*4    ) -= four_visc2_timetauMp*(viscs2_(0,ui)*viscs2_(2,vi)+viscs2_(1,ui)*viscs2_(5,vi)+viscs2_(2,ui)*viscs2_(8,vi)) ;

                //estif(vi*4,     ui*4 + 1) -= four_visc2_timetauMp*(viscs2_(0,vi)*viscs2_(1,ui)+viscs2_(1,vi)*viscs2_(4,ui)+viscs2_(2,vi)*viscs2_(5,ui)) ;
                //estif(vi*4 + 1, ui*4 + 1) -= four_visc2_timetauMp*(viscs2_(1,ui)*viscs2_(1,vi)+viscs2_(4,ui)*viscs2_(4,vi)+viscs2_(5,ui)*viscs2_(5,vi)) ;
                //estif(vi*4 + 2, ui*4 + 1) -= four_visc2_timetauMp*(viscs2_(1,ui)*viscs2_(2,vi)+viscs2_(4,ui)*viscs2_(5,vi)+viscs2_(5,ui)*viscs2_(8,vi)) ;

                //estif(vi*4,     ui*4 + 2) -= four_visc2_timetauMp*(viscs2_(0,vi)*viscs2_(2,ui)+viscs2_(1,vi)*viscs2_(5,ui)+viscs2_(2,vi)*viscs2_(8,ui)) ;
                //estif(vi*4 + 1, ui*4 + 2) -= four_visc2_timetauMp*(viscs2_(1,vi)*viscs2_(2,ui)+viscs2_(4,vi)*viscs2_(5,ui)+viscs2_(5,vi)*viscs2_(8,ui)) ;
                //estif(vi*4 + 2, ui*4 + 2) -= four_visc2_timetauMp*(viscs2_(2,ui)*viscs2_(2,vi)+viscs2_(5,ui)*viscs2_(5,vi)+viscs2_(8,ui)*viscs2_(8,vi)) ;
              } // vi
            } // ui
            }  // end for(idim)

            if (is_newton_)
            {
              for(int idim=0;idim<nsd_;++idim)
              {
              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui = numdofpernode_*ui;
                double v = two_visc_timetauMp*densaf_*funct_(ui);
                for (int vi=0; vi<nen_; ++vi)
                {
                  const int fvi=numdofpernode_*vi+idim;
                  /* viscous stabilisation, reactive part of convection */
                  /*
                    /                                         \
                    |  /                \   n+1               |
                +/- | |   rho*Du o nabla | u    , div eps (v) |
                    |  \                /   (i)               |
                    \                                         /
                  */
                  for (int jdim=0; jdim<nsd_; ++jdim)
                    estif(fvi, fui  ) += v*viscs2_(jdim+(idim*nsd_),vi)*vderxy_(jdim,idim);

                  /*
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
                                                  viscs2_(8,vi)*vderxy_(2,2)) ;*/
                } // vi
              } // ui
              } // end for(idim)
            } // if newton
#endif
      //----------------------------------------------------------------------
      //                     STABILISATION, CONTINUITY PART

      if (cstab_ == INPAR::FLUID::continuity_stab_yes)
      {
        const double timetauC = timefac*tau_C;
        for (int idim = 0; idim < nsd_;++idim)
        {
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui   = numdofpernode_*ui + idim;
            //const int fui   = 4*ui;
            //const int fuip  = fui+1;
            //const int fuipp = fui+2;
            double v_idim = timetauC*derxy_(idim, ui);
            //double v0 = timetauC*derxy_(0, ui);
            //double v1 = timetauC*derxy_(1, ui);
            //double v2 = timetauC*derxy_(2, ui);
            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi   = numdofpernode_*vi;
              //const int fvi   = 4*vi;
              //const int fvip  = fvi+1;
              //const int fvipp = fvi+2;
              /* continuity stabilisation on left hand side */
              /*
                /                         \
                |                          |
                | nabla o Du  , nabla o v  |
                |                          |
                \                         /
              */
              for (int jdim = 0; jdim < nsd_; ++jdim)
                estif(fvi+jdim,  fui  ) += v_idim*derxy_(jdim, vi);

              //estif(fvi,  fui  ) += v0*derxy_(0, vi) ;
              //estif(fvip, fui  ) += v0*derxy_(1, vi) ;
              //estif(fvipp,fui  ) += v0*derxy_(2, vi) ;

              //estif(fvi,  fuip ) += v1*derxy_(0, vi) ;
              //estif(fvip, fuip ) += v1*derxy_(1, vi) ;
              //estif(fvipp,fuip ) += v1*derxy_(2, vi) ;

              //estif(fvi,  fuipp) += v2*derxy_(0, vi) ;
              //estif(fvip, fuipp) += v2*derxy_(1, vi) ;
              //estif(fvipp,fuipp) += v2*derxy_(2, vi) ;
            }
          }

          const double tauC_conres = tau_C*conres_old_;
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = numdofpernode_*vi+idim;
            //const int fvi   = 4*vi;
            //const int fvip  = fvi+1;
            //const int fvipp = fvi+2;
            // continuity stabilisation on right hand side
            eforce(fvi  ) -= tauC_conres*derxy_(idim, vi) ;
            //eforce(fvi  ) -= tauC_conres*derxy_(0, vi) ;
            //eforce(fvip ) -= tauC_conres*derxy_(1, vi) ;
            //eforce(fvipp) -= tauC_conres*derxy_(2, vi) ;
          }
        }  // end  for(idim)
      }  // end Continuity stabilization

      //TODO: Cross-stress and Reynolds-Stress stabilization are missing in the Sysmat implementation

      //----------------------------------------------------------------------
      //     FINE-SCALE SUBGRID-VISCOSITY TERM (ON RIGHT HAND SIDE)

      //TODO: Difference steady <-> unsteady
      // TODO: Stabilization of FINE-SCALE SUBGRID-VISCOSITY TERM not yet variable in space dimensions

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
    } // end Galerkin Part

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

  return;
}


/*!
      \brief Do a finite difference check for a given element id --- 
      this function is for debugging purposes only

      \param ele              (i) the element those matrix is calculated 
                                  (pass-through)
      \param evelaf           (i) nodal velocities at n+alpha_F/n+1 (pass-through)
      \param eveln            (i) nodal velocities at n (pass-through)
      \param fsevelaf         (i) fine-scale nodal velocities at n+alpha_F/n+1 
                                  (pass-through)
      \param epreaf           (i) nodal pressure at n+alpha_F/n+1 (pass-through)
      \param eaccam           (i) nodal accelerations at n+alpha_M (pass-through)
      \param escaaf           (i) nodal scalar at n+alpha_F/n+1 (pass-through)
      \param escaam           (i) nodal scalar at n+alpha_M/n (pass-through)
      \param escadtam         (i) nodal scalar derivatives at n+alpha_M/n+1 
                                  (pass-through)
      \param emhist           (i) time rhs for momentum equation (pass-through)
      \param edispnp          (i) nodal displacements (on moving mesh) 
                                  (pass-through)
      \param egridv           (i) grid velocity (on moving mesh) (pass-through)
      \param estif            (i) element matrix to calculate (pass-through)
      \param emesh            (i) linearization wrt mesh motion (pass-through)
      \param eforce           (i) element rhs to calculate (pass-through)
      \param material         (i) fluid material (pass-through)
      \param time             (i) current simulation time (pass-through)
      \param timefac          (i) time discretization factor (pass-through)
      \param newton           (i) boolean flag for linearisation (pass-through)
      \param loma             (i) boolean flag for potential low-Mach-number solver 
                                  (pass-through)
      \param conservative     (i) boolean flag for conservative form (pass-through)
      \param is_genalpha      (i) boolean flag for generalized-alpha time 
                                  integration (pass-through)
      \param higher_order_ele (i) keep or drop second derivatives (pass-through)
      \param fssgv            (i) flag for type of fine-scale subgrid viscosity 
                                  (pass-through)
      \param pspg             (i) boolean flag for stabilisation (pass-through)
      \param supg             (i) boolean flag for stabilisation (pass-through)
      \param vstab            (i) boolean flag for stabilisation (pass-through)
      \param cstab            (i) boolean flag for stabilisation (pass-through)
      \param cross            (i) boolean flag for stabilisation (pass-through)
      \param reynolds         (i) boolean flag for stabilisation (pass-through)
      \param whichtau         (i) selecting available tau definition (pass-through)
      \param turb_mod_action  (i) selecting turbulence model (none, Smagorisky, 
                                  dynamic Smagorinsky, Smagorinsky with van Driest
                                  damping for channel flows) (pass-through)
      \param Cs               (i) Smagorinsky model parameter (pass-through)
      \param Cs_delta_sq      (i) Model parameter computed by dynamic Smagorinsky
                                  approach (Cs*h*h) (pass-through)
      \param l_tau            (i) viscous length scale, required for van driest
                                  damping function and defined on input (pass-through)
*/            
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::Fluid3Impl<distype>::FDcheck(
  Fluid3*                                               ele,
  const LINALG::Matrix<nsd_,nen_>&                      evelaf,
  const LINALG::Matrix<nsd_,nen_>&                      eveln,
  const LINALG::Matrix<nsd_,nen_>&                      fsevelaf,
  const LINALG::Matrix<nen_,1>&                         epreaf,
  const LINALG::Matrix<nsd_,nen_>&                      eaccam,
  const LINALG::Matrix<nen_,1>&                         escaaf,
  const LINALG::Matrix<nen_,1>&                         escaam,
  const LINALG::Matrix<nen_,1>&                         escadtam,
  const LINALG::Matrix<nsd_,nen_>&                      emhist,
  const LINALG::Matrix<nsd_,nen_>&                      edispnp,
  const LINALG::Matrix<nsd_,nen_>&                      egridv,
  const LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&    estif,
  const LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&    emesh,
  const LINALG::Matrix<(nsd_+1)*nen_,    1>&            eforce,
  const double                                          thermpressaf,
  const double                                          thermpressam,
  const double                                          thermpressdtam,
  const Teuchos::RCP<const MAT::Material>               material,
  const double                                          timefac,
  const double&                                         Cs,
  const double&                                         Cs_delta_sq,
  const double&                                         l_tau)
{
  // magnitude of dof perturbation
  const double epsilon=1e-9;

  // make a copy of all input parameters potentially modified by Sysmat 
  // call --- they are not intended to be modified
  double copy_Cs         =Cs;
  double copy_Cs_delta_sq=Cs_delta_sq;
  double copy_l_tau      =l_tau;

  Teuchos::RCP<const MAT::Material> copy_material=material;

  // allocate arrays to compute element matrices and vectors at perturbed
  // positions
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> checkmat1(true);
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> checkmat2(true);
  LINALG::Matrix<(nsd_+1)*nen_,            1> checkvec1(true);
        
  // alloc the vectors that will contain the perturbed velocities or 
  // pressures 
  LINALG::Matrix<nsd_,nen_>                   checkevelaf(true);
  LINALG::Matrix<nsd_,nen_>                   checkeaccam(true);
  LINALG::Matrix<nen_,1>                      checkepreaf(true);

  // echo to screen
  printf("+-------------------------------------------+\n");
  printf("| FINITE DIFFERENCE CHECK FOR ELEMENT %5d |\n",ele->Id());
  printf("+-------------------------------------------+\n");
  printf("\n");
  // loop columns of matrix by looping nodes and then dof per nodes

  // loop nodes
  for(int nn=0;nn<nen_;++nn)
  {
    printf("-------------------------------------\n");
    printf("-------------------------------------\n");
    printf("NODE of element local id %d\n",nn);
    // loop dofs
    for(int rr=0;rr<(nsd_+1);++rr)
    {
      // number of the matrix column to check
      int dof=nn*(nsd_+1)+rr;

      // clear element matrices and vectors to assemble
      checkmat1.Clear();
      checkmat2.Clear();
      checkvec1.Clear();

      // copy velocities and pressures to perturbed arrays
      for(int mm=0;mm<nen_;++mm)
      {
        for(int dim=0;dim<nsd_;++dim)
        {
          checkevelaf(dim,mm)=evelaf(dim,mm);
                
          checkeaccam(dim,mm)=eaccam(dim,mm);
        }
              
        checkepreaf(  mm)=epreaf(  mm);
      }

      // perturb the respective elemental quantities
      if(rr==nsd_)
      {
        printf("pressure dof (%d) %f\n",nn,epsilon);

        if (is_genalpha_)
        {
          checkepreaf(nn)+=alphaF_*epsilon;
        }
        else
        {
          checkepreaf(nn)+=epsilon;
        }
      }
      else
      {
        printf("velocity dof %d (%d)\n",rr,nn);

        if (is_genalpha_)
        {
          checkevelaf(rr,nn)+=alphaF_*epsilon;
          checkeaccam(rr,nn)+=alphaM_/(gamma_*dt_)*epsilon;
        }
        else
        {
          checkevelaf(rr,nn)+=epsilon;
        }
      }
        
      // calculate the right hand side for the perturbed vector
      Sysmat2D3D(ele,
                 checkevelaf,
                 eveln,
                 fsevelaf,
                 checkepreaf,
                 checkeaccam,
                 escaaf,
                 escaam,
                 escadtam,
                 emhist,
                 edispnp,
                 egridv,
                 checkmat1,
                 checkmat2,
                 checkvec1,
                 thermpressaf,
                 thermpressam,
                 thermpressdtam,
                 copy_material,
                 timefac,
                 copy_Cs,
                 copy_Cs_delta_sq,
                 copy_l_tau);

      // compare the difference between linaer approximation and
      // (nonlinear) right hand side evaluation

      // note that it makes more sense to compare these quantities
      // than to compare the matrix entry to the difference of the
      // the right hand sides --- the latter causes numerical problems
      // do to deletion

      for(int mm=0;mm<(nsd_+1)*nen_;++mm)
      {
        double val;
        double lin;
        double nonlin;

        // For af-generalized-alpha scheme, the residual vector for the 
        // solution rhs is scaled on the time-integration level...
        if (is_genalpha_)
        {
          val   =-eforce(mm)   /(epsilon);
          lin   =-eforce(mm)   /(epsilon)+estif(mm,dof)*(alphaM_)/(gamma_*dt_);
          nonlin=-checkvec1(mm)/(epsilon);
        }
        else
        {
          val   =-eforce(mm)/epsilon;
          lin   =-eforce(mm)/epsilon+estif(mm,dof);
          nonlin=-checkvec1(mm)/epsilon;
        }

        double norm=abs(lin);
        if(norm<1e-9)
        {
          norm=1e-9;
        }
              
        // output to screen
        printf("relerr         %+12.5e ",(lin-nonlin)/norm);
        printf("abserr         %+12.5e ",lin-nonlin);
        printf("orig. value    %+12.5e ",val);
        printf("lin. approx.   %+12.5e ",lin);
        printf("nonlin. funct. %+12.5e ",nonlin);
        printf("matrix entry   %+12.5e ",estif(mm,dof));
        printf("\n");
      }
    }
  }

  return;
}

#endif
#endif

