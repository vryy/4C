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

  //LINALG::Matrix<nen_,1>      preforce;
  //LINALG::Matrix<nsd_,nen_>   velforce;

  //LINALG::Matrix<nen_,nen_>   ppmat;

  //LINALG::Matrix<nen_*nsd_,nen_*nsd_>    estif_(true);
  estif_.Clear();

  estif_p_v_.Clear();
  estif_q_u_.Clear();

  preforce.Clear();
  velforce.Clear();

  ppmat.Clear();

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

  // if not available, the arrays for the subscale quantities have to
  // be resized and initialised to zero
  if(tds_)
  {
    ele->ActivateTDS(intpoints.IP().nquad,nsd_);
  }

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

    /*-------------------------------------------------------------------*
     *                                                                   *
     *                  get residual of momentum equation                *
     *                                                                   *
     *-------------------------------------------------------------------*/   

    if (is_genalpha_)
    {
      // rhs of momentum equation: density*bodyforce at n+alpha_F
      rhsmom_.Update(densaf_,bodyforce_,0.0);

      // get acceleration at time n+alpha_M at integration point
      accint_.Multiply(eaccam,funct_);

      // evaluate momentum residual once for all stabilization right hand sides
      for (int rr=0;rr<nsd_;++rr)
      {
        momres_old_(rr) = 
          densam_*accint_(rr)
          +
          densaf_*conv_old_(rr)
          +
          gradp_(rr)
          -
          2*visceff_*visc_old_(rr)
          -
          densaf_*bodyforce_(rr);
      }
    }
    else
    {
      // rhs of momentum equation:
      // density*timefac*bodyforce at n+1 + density*histmom at n

      // in the case of a Boussinesq approximation: 
      //
      //                    f = (rho - rho_0)/rho_0 *g
      //
      // else: 		 
      //  
      //                           f = rho * g

      if (physicaltype_ == INPAR::FLUID::boussinesq)
      {
        rhsmom_.Update(densn_,histmom_,deltadens_*timefac,bodyforce_);
      }
      else 
      {
        rhsmom_.Update(densn_,histmom_,densaf_*timefac,bodyforce_);
      }

      // modify integration factor for Galerkin rhs
      rhsfac *= timefac;

      // evaluate momentum residual once for all stabilization right hand sides

      // stationary  : momres_old = theta ( ... ) -  bodyforce_
      // instationary: momres_old = u_(n+1) + theta ( ... ) - histmom_ - bodyforce_
      for (int rr=0;rr<nsd_;++rr)
      {
        momres_old_(rr) = 
          timefac*(densaf_*conv_old_(rr)+gradp_(rr)-2*visceff_*visc_old_(rr))
          -
          rhsmom_(rr);
      }

      if (is_stationary_ == false)
      {
        for (int rr=0;rr<nsd_;++rr)
        {
          momres_old_(rr) += densaf_*velint_(rr);
        }
      }
    }

    /*-------------------------------------------------------------------*
     *                                                                   *
     *                  update of SUBSCALE VELOCITY                      *
     *                                                                   *
     *-------------------------------------------------------------------*/
    double fac1   =0.0;
    double fac2   =0.0;
    double fac3   =0.0;
    double facMtau=0.0;

    /*-------------------------------------------------------------------*
     *                                                                   *
     *               time dependent subgrid scale closure                *
     *                                                                   *
     *-------------------------------------------------------------------*/
    if(tds_)
    {
      /*
                                              1.0
         facMtau =  -------------------------------------------------------
                       n+aM                      n+aF
                    rho     * alphaM * tauM + rho     * alphaF * gamma * dt
      */
      facMtau = alphaM_/(densam_*tau_(0)+densaf_*timefac);

      /* 
         factor for old subgrid velocities:

                   n+aM                      n+aF
         fac1 = rho     * alphaM * tauM + rho     * gamma * dt * (alphaF-1)
      */
      fac1=(densam_*alphaM_*tau_(0)+densaf_*gamma_*dt_*(alphaF_-1.0));
      /*
        factor for old subgrid accelerations

                   n+aM
         fac2 = rho     * tauM * dt * (alphaM-gamma)
      */
      fac2=(densam_*dt_*tau_(0)*(alphaM_-gamma_));
      /*
        factor for residual in current subgrid velocities:

         fac3 = gamma * dt * tauM
      */
      fac3=(gamma_*dt_*tau_(0));

      // if no generalised alpha time integration is used, the momres_old_
      // contains a scaled residual
      if (!is_genalpha_)
      {
        fac3/=timefac;
      }

      /*                 +-                                       -+
        ~n+1             |        ~n           ~ n            n+1  |
        u    = facMtau * | fac1 * u  + fac2 * acc  -fac3 * res     |  
         (i)             |                                    (i)  |
                         +-                                       -+
      */

      /* compute the intermediate value of subscale velocity

              ~n+af            ~n+1                   ~n
              u     = alphaF * u     + (1.0-alphaF) * u
               (i)              (i)

      */

      for (int rr=0;rr<nsd_;++rr)
      {
        ele->UpdateSvelnpInOneDirection(
          fac1           ,
          fac2           ,
          fac3           ,
          facMtau        ,
          momres_old_(rr),
          alphaF_        ,
          rr             ,
          iquad          ,
          sgvelint_(rr)  );
      }
    }
    /*-------------------------------------------------------------------*
     *                                                                   *
     *                     classical subgrid closure                     *
     *                                                                   *
     *-------------------------------------------------------------------*/
    else
    {
      if (cross_    != INPAR::FLUID::cross_stress_stab_none or
          reynolds_ != INPAR::FLUID::reynolds_stress_stab_none)
      {
        if (is_genalpha_)
        {
          // compute subgrid-scale velocity
          sgvelint_.Update(-tau_(1),momres_old_,0.0);
        }
        else
        {
          // compute subgrid-scale velocity
          sgvelint_.Update(-(tau_(1)/timefac),momres_old_,0.0);
          
          // sgvelint_.Update(-(tau_(1)/dt_),momres_old_,0.0);
        }
      }
    }

    if (cross_    != INPAR::FLUID::cross_stress_stab_none or
        reynolds_ != INPAR::FLUID::reynolds_stress_stab_none)
    {
          
          // compute subgrid-scale convective operator
          sgconv_c_.MultiplyTN(derxy_,sgvelint_);
          
          // re-calculate convective term from previous iteration if cross-stress
          // is included
          convvelint_.Update(1.0,sgvelint_,1.0);
          conv_old_.Multiply(vderxy_,convvelint_);
    }
    else 
    {
      sgconv_c_.Clear();
    }

    /*-------------------------------------------------------------------*
     *                                                                   *
     *                 get residual of continuity equation               *
     *                                                                   *
     *-------------------------------------------------------------------*/   
    if (is_genalpha_)
    {
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
    //{
    lin_resM_Du_.Clear();

    //----------------------------------------------------------------------
    //  PROVIDE LINEARISATION OF GALERKIN MOMENTUM RESIDUAL WRT VELOCITIES

    /*
                   /       n+1        \       /                \  n+1
         rho*Du + |   rho*u   o nabla | Du + |   rho*Du o nabla | u
                   \      (i)         /       \                /   (i)
    */


    int idim_nsd_p_idim[nsd_];

    for (int idim = 0; idim <nsd_; ++idim)
    {
      idim_nsd_p_idim[idim]=idim*nsd_+idim;
    }

    if (is_stationary_ == false)
    {
      const double fac_densam=fac_*densam_;

      for (int ui=0; ui<nen_; ++ui)
      {
        const double v=fac_densam*funct_(ui);

        for (int idim = 0; idim <nsd_; ++idim)
        {
          lin_resM_Du_(idim_nsd_p_idim[idim],ui)+=v;
        }
      }
    }

    const double timefacfac_densaf=timefacfac*densaf_;

    for (int ui=0; ui<nen_; ++ui)
    {
      const double v=timefacfac_densaf*(conv_c_(ui)+sgconv_c_(ui));

      for (int idim = 0; idim <nsd_; ++idim)
      {
        lin_resM_Du_(idim_nsd_p_idim[idim],ui)+=v;
      }
    }

    if(is_newton_)
    {
      for (int ui=0; ui<nen_; ++ui)
      {
        const double temp=timefacfac_densaf*funct_(ui);

        for (int idim = 0; idim <nsd_; ++idim)
        {
          const int idim_nsd=idim*nsd_;

          for(int jdim=0;jdim<nsd_;++jdim)
          {
            lin_resM_Du_(idim_nsd+jdim,ui)+=temp*vderxy_(idim,jdim);
          }
        }
      }
    }


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
    if(is_newton_)
    {
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui   = nsd_*ui;

        for (int idim = 0; idim <nsd_; ++idim)
        {
          const int idim_nsd=idim*nsd_;

          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = nsd_*vi;

            const int fvi_p_idim = fvi+idim;

            for (int jdim= 0; jdim<nsd_;++jdim)
            {
              estif_(fvi_p_idim,fui+jdim) += funct_(vi)*lin_resM_Du_(idim_nsd+jdim,ui);
            } // end for (jdim)
          } // end for (idim)
        } //vi
      } // ui
    }
    else
    {
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui   = nsd_*ui;

        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = nsd_*vi;

          for (int idim = 0; idim <nsd_; ++idim)
          {
            estif_(fvi+idim,fui+idim) += funct_(vi)*lin_resM_Du_(idim*nsd_+idim,ui);
          } // end for (idim)
        } //vi
      } // ui
    }


    // inertia terms on the right hand side for instationary fluids
    if (is_stationary_ == false)
    {
      double vel=0.0;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        if (is_genalpha_)
        {
          vel=-fac_*densam_*accint_(idim);
        }
        else
        {
          vel=-fac_*densaf_*velint_(idim);
        }

        for (int vi=0; vi<nen_; ++vi)
        {
          /* inertia term on right-hand side for generalized-alpha scheme */
          velforce(idim,vi)+= vel*funct_(vi);
        }
      }
    }  // end if(stationary)

    for (int idim = 0; idim <nsd_; ++idim)
    {
      const double scaled_conv_old=-rhsfac*densaf_*conv_old_(idim);

      for (int vi=0; vi<nen_; ++vi)
      {
        /* convection (convective form) on right-hand side */
        velforce(idim,vi) += scaled_conv_old*funct_(vi);
      }
    }  // end for(idim)


    //----------------------------------------------------------------------
    //                       VISCOUS GALERKIN PART
    const double visceff_timefacfac = visceff_*timefacfac;

    /* viscosity term */
    /*
                     /                        \
                    |       /  \         / \   |
              2 mu  |  eps | Du | , eps | v |  |
                    |       \  /         \ /   |
                     \                        /
    */

    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi   = nsd_*vi;

      for (int jdim= 0; jdim<nsd_;++jdim)
      {
        const double temp=visceff_timefacfac*derxy_(jdim,vi);

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = nsd_*ui;

          for (int idim = 0; idim <nsd_; ++idim)
          {
            const int fvi_p_idim = fvi+idim;

            estif_(fvi_p_idim,fui+jdim) += temp*derxy_(idim, ui);

          } // end for (jdim)
        } // end for (idim)
      } // ui
    } //vi
    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi   = nsd_*vi;

      for (int jdim= 0; jdim<nsd_;++jdim)
      {
        const double temp=visceff_timefacfac*derxy_(jdim,vi);

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui   = nsd_*ui;

          for (int idim = 0; idim <nsd_; ++idim)
          {
            const int fvi_p_idim = fvi+idim;

            estif_(fvi_p_idim,fui+idim) += temp*derxy_(jdim, ui);

          } // end for (jdim)
        } // end for (idim)
      } // ui
    } //vi

    const double v = -visceff_*rhsfac;

    for (int jdim = 0; jdim < nsd_; ++jdim)
    {
      for (int idim = 0; idim < nsd_; ++idim)
      {
        const double temp=v*(vderxy_(jdim,idim)+vderxy_(idim,jdim));

        //computation of right-hand-side viscosity term
        for (int vi=0; vi<nen_; ++vi)
        {
          /* viscosity term on right-hand side */
          velforce(idim,vi)+= temp*derxy_(jdim,vi);
        }
      }
    }  // end for(idim) viscosity term

    //----------------------------------------------------------------------
    //   STABILISATION, CONTINUITY PART AND GALERKIN, VISCOUS PART (LOMA)
    //----------------------------------------------------------------------
    double conti_stab_and_vol_visc_fac=0.0;
    double conti_stab_and_vol_visc_rhs=0.0;

    if(cstab_ == INPAR::FLUID::continuity_stab_yes)
    {
      conti_stab_and_vol_visc_fac+=timefac*tau_C;
      conti_stab_and_vol_visc_rhs-=tau_C*conres_old_;
    }
    if (physicaltype_ == INPAR::FLUID::loma)
    {
      conti_stab_and_vol_visc_fac-=(2.0/3.0)*visceff_*timefacfac;
      conti_stab_and_vol_visc_rhs+=(2.0/3.0)*visceff_*rhsfac*vdiv_;
    }


    /* continuity stabilisation on left hand side */
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
    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui   = nsd_*ui;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        const int fui_p_idim = fui+idim;
        const double v0 = conti_stab_and_vol_visc_fac*derxy_(idim,ui);

        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = nsd_*vi;

          for(int jdim=0;jdim<nsd_;++jdim)
          {
            estif_(fvi+jdim,fui_p_idim) += v0*derxy_(jdim, vi) ;
          }
        }
      }
    }  // end for(idim)

    for (int vi=0; vi<nen_; ++vi)
    {
      //const int fvi = numdofpernode_*vi;

      for(int idim=0;idim<nsd_;++idim)
      {
        velforce(idim,vi)+= derxy_(idim,vi)*conti_stab_and_vol_visc_rhs;
      }
    }

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
          const int fui   = nsd_*ui + idim;
          //const int fui   = 4*ui;
          //const int fuip  = fui+1;
          //const int fuipp = fui+2;
          double v = timefacfac*densaf_*funct_(ui)*vdiv_;
          if (physicaltype_ == INPAR::FLUID::loma) v -= timefacfac*densaf_*scaconvfacaf_*conv_scaaf_;
          // only with linear density-concentration correlation
          else if(physicaltype_ == INPAR::FLUID::varying_density) v += timefacfac*conv_scaaf_;
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = nsd_*vi + idim;
            //const int fvi   = 4*vi;
            //const int fvip  = fvi+1;
            //const int fvipp = fvi+2;
            /* convection, convective part (conservative addition) */
            /*
              /                                                \
              |      /              n+1    n+1           \      |
              |  Du | rho*nabla o u    +  u   *nabla rho | , v  |
              |      \             (i)     (i)          /       |
              \                                                 /
            */
            double v2 = v*funct_(vi) ;
            estif_(fvi  , fui  ) += v2;
            //estif(fvi  , fui  ) += v2;
            //estif(fvip , fuip ) += v2;
            //estif(fvipp, fuipp) += v2;
          }
        }

        if (is_newton_)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            const int fvi   = nsd_*vi + idim;
            //const int fvi   = 4*vi;
            //const int fvip  = fvi+1;
            //const int fvipp = fvi+2;
            const double v_idim = timefacfac*densaf_*velint_(idim)*funct_(vi);
            //const double v0 = timefacfac*densaf_*velint_(0)*funct_(vi);
            //const double v1 = timefacfac*densaf_*velint_(1)*funct_(vi);
            //const double v2 = timefacfac*densaf_*velint_(2)*funct_(vi);
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui   = nsd_*ui;
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
              estif_(fvi,  fui+jdim  ) += v_idim*derxy_(jdim, ui) ;

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
              const int fvi   = nsd_*vi + idim;
              //const int fvi   = 4*vi;
              //const int fvip  = fvi+1;
              //const int fvipp = fvi+2;
              const double v_idim = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(idim)*velint_(idim)*funct_(vi);
              //const double v0 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(0)*velint_(0)*funct_(vi);
              //const double v1 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(1)*velint_(1)*funct_(vi);
              //const double v2 = -timefacfac*densaf_*scaconvfacaf_*grad_scaaf_(2)*velint_(2)*funct_(vi);
              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui   = nsd_*ui;
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
                  estif_(fvi,  fui +jdim  ) += v_idim*funct_(ui) ;

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
              const int fvi   = nsd_*vi + idim;
              //const int fvi   = 4*vi;
              //const int fvip  = fvi+1;
              //const int fvipp = fvi+2;
              const double v_idim = +timefacfac*grad_scaaf_(idim)*velint_(idim)*funct_(vi);
              //const double v0 = +timefacfac*grad_scaaf_(0)*velint_(0)*funct_(vi);
              //const double v1 = +timefacfac*grad_scaaf_(1)*velint_(1)*funct_(vi);
              //const double v2 = +timefacfac*grad_scaaf_(2)*velint_(2)*funct_(vi);
              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui   = nsd_*ui;
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
                  estif_(fvi,  fui+jdim  ) += v_idim*funct_(ui) ;

              }
            }
          }
        }

        for (int vi=0; vi<nen_; ++vi)
        {
          const int fvi   = nsd_*vi +idim;
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
            const int fvi   = nsd_*vi+ idim;
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
            const int fvi   = nsd_*vi + idim;
            //const int fvi   = 4*vi;
            /* convection (conservative addition) on rhs for low-Mach-number flow */
            double v = -rhsfac*conv_scaaf_*funct_(vi);
            eforce(fvi    ) += v*velint_(idim) ;
          }
        }
      }  // end for(idim)
    }  // end if(is_conservative)

    //----------------------------------------------------------------------
    // computation of pressure term including right-hand-side contribution
    //----------------------------------------------------------------------
    for (int ui=0; ui<nen_; ++ui)
    {
      //const int fuippp = numdofpernode_*ui+nsd_;
      //const int fuippp = 4*ui+3;
      const double v = -timefacfac*funct_(ui);
      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = nsd_*vi;
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
          estif_p_v_(fvi + idim,ui) += v*derxy_(idim, vi);
        }
      }
    }

    const double pressfac = press*rhsfac;

    for (int vi=0; vi<nen_; ++vi)
    {
      /* pressure term on right-hand side */
      for (int idim = 0; idim <nsd_; ++idim)
      {
        velforce(idim,vi)+= pressfac*derxy_(idim, vi) ;
      }
    }  //end for(idim)

    //----------------------------------------------------------------------
    // computation of continuity term including right-hand-side contribution
    //----------------------------------------------------------------------
    for (int vi=0; vi<nen_; ++vi)
    {
      //const int fvippp = numdofpernode_*vi+nsd_;

      const double v = timefacfac*funct_(vi);
      for (int ui=0; ui<nen_; ++ui)
      {
        const int fui   = nsd_*ui;

        for (int idim = 0; idim <nsd_; ++idim)
        {
          /* continuity term */
          /*
               /                \
              |                  |
              | nabla o Du  , q  |
              |                  |
               \                /
          */
          estif_q_u_(vi,fui+idim) += v*derxy_(idim,ui);
        }
      }
    }  // end for(idim)

    const double rhsfac_vdiv = -rhsfac * vdiv_;
    for (int vi=0; vi<nen_; ++vi)
    {
      // continuity term on right-hand side
      preforce(vi) += rhsfac_vdiv*funct_(vi);
    }

    //----------------------------------------------------------------------
    // computation of additional terms for low-Mach-number flow:
    // 2) additional rhs term of continuity equation
    //----------------------------------------------------------------------
    if (physicaltype_ == INPAR::FLUID::loma)
    {
      if (is_newton_)
      {
        const double timefacfac_scaconvfacaf=timefacfac*scaconvfacaf_;

        for (int ui=0; ui<nen_; ++ui)
        {
          const int fui=nsd_*ui;

          const double timefacfac_scaconvfacaf_funct_ui=timefacfac_scaconvfacaf*funct_(ui);

          for(int jdim=0;jdim<nsd_;++jdim)
          {
            const double temp=timefacfac_scaconvfacaf_funct_ui*grad_scaaf_(jdim);

            for (int vi=0; vi<nen_; ++vi)
            {
              //const int fvippp= numdofpernode_*vi+nsd_;


              /*
                  factor afgtd/am

                          /                    \
                    1    |       /         \    |
                   --- * |  q , | Du o grad | T |
                    T    |       \         /    |
                          \                    /
              */
              estif_q_u_(vi,fui+jdim) -= temp*funct_(vi);
            }
          }
        }
      } // end if (is_newton_)

      const double fac_rhscon = fac_*rhscon_;
      for (int vi=0; vi<nen_; ++vi)
      {
        /* additional rhs term of continuity equation */
        preforce(vi) += fac_rhscon*funct_(vi) ;
      }
    }  // end if(loma)

    //----------------------------------------------------------------------
    // computation of body-force term on right-hand side
    //----------------------------------------------------------------------
    for (int idim = 0; idim <nsd_; ++idim)
    {
      const double scaled_rhsmom=fac_*rhsmom_(idim);

      for (int vi=0; vi<nen_; ++vi)
      {
        velforce(idim,vi)+=scaled_rhsmom*funct_(vi);
      }
    }  // end for(idim)

    //----------------------------------------------------------------------
    //  PROVIDE LINEARISATION OF MOMENTUM RESIDUAL WRT VELOCITIES FOR ALL
    //       STABILISATION PARTS (EXTEND GALERKIN RES IF NECESSARY)

    /*
                   /       n+1       \        /                \  n+1
         rho*Du + |   rho*u   o nabla | Du + |   rho*Du o nabla | u   +
                   \      (i)        /        \                /   (i)

                               /  \
                + nabla o eps | Du |
                               \  /
    */
    if (is_higher_order_ele_)
    {
      const double v = -2.0*visceff_*timefacfac;
      for (int idim = 0; idim <nsd_; ++idim)
      {
        const int nsd_idim=nsd_*idim;

        for(int jdim=0;jdim<nsd_;++jdim)
        {
          const int nsd_idim_p_jdim=nsd_idim+jdim;

          for (int ui=0; ui<nen_; ++ui)
          {
            lin_resM_Du_(nsd_idim_p_jdim,ui)+=v*viscs2_(nsd_idim_p_jdim, ui);
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
        for(int jdim=0;jdim<nsd_;++jdim)
        {
          for (int ui=0; ui<nen_; ++ui)
          {
            const int fui_p_jdim   = nsd_*ui + jdim;

            for(int idim=0;idim<nsd_;++idim)
            {
              const int nsd_idim=nsd_*idim;

              for (int vi=0; vi<nen_; ++vi)
              {
                //const int fvippp = numdofpernode_*vi+nsd_;

                const double temp_vi_idim=derxy_(idim,vi)*tau_(1);

                estif_q_u_(vi,fui_p_jdim) += lin_resM_Du_(nsd_idim+jdim,ui)*temp_vi_idim;
              } // jdim
            } // vi
          } // ui
        } //idim
      } // end if (is_higher_order_ele_) or (newton_)
      else
      {
        for (int vi=0; vi<nen_; ++vi)
        {
          //const int fvippp = numdofpernode_*vi+nsd_;

          for(int idim=0;idim<nsd_;++idim)
          {
            const int nsd_idim=nsd_*idim;

            const double temp_vi_idim=derxy_(idim, vi)*tau_(1);

            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui_p_idim   = nsd_*ui + idim;

              estif_q_u_(vi,fui_p_idim) += lin_resM_Du_(nsd_idim+idim,ui)*temp_vi_idim;
            } // vi
          } // ui
        } //idim
      } // end if not (is_higher_order_ele_) nor (newton_)


      for (int ui=0; ui<nen_; ++ui)
      {

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
            ppmat(vi,ui)+=timetauMp_derxy_idim_ui*derxy_(idim,vi);
          } // vi
        } // end for(idim)
      }  // ui

      for (int idim = 0; idim <nsd_; ++idim)
      {
        const double temp= tau_Mp*momres_old_(idim);

        for (int vi=0; vi<nen_; ++vi)
        {
          // pressure stabilisation
          preforce(vi) -= temp*derxy_(idim, vi);
        }
      } // end for(idim)
    }  // end PSPG stabilization

      //----------------------------------------------------------------------
      //                     SUPG STABILISATION PART

      if(supg_ == INPAR::FLUID::convective_stab_supg)
      {
        LINALG::Matrix<nen_,1> supg_test;

        const double supgfac=densaf_*tau_(0);
        for (int vi=0; vi<nen_; ++vi)
        {
          supg_test(vi)=supgfac*(conv_c_(vi)+sgconv_c_(vi));
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
              |    /       n+1        \        /     n+1        \     |
              |   |   rho*u    o nabla | Du , | rho*u    o nabla | v  |
              |    \       (i)        /        \     (i)        /     |
               \                                                     /
        */
        if (is_higher_order_ele_ || is_newton_)
        {
          for (int vi=0; vi<nen_; ++vi)
          {
            for(int idim=0;idim<nsd_;++idim)
            {
              const int nsd_idim=nsd_*idim;

              const int fvi_p_idim = nsd_*vi+idim;

              for(int jdim=0;jdim<nsd_;++jdim)
              {
                const int nsd_idim_p_jdim=nsd_idim+jdim;
                for (int ui=0; ui<nen_; ++ui)
                {
                  const int fui_p_jdim   = nsd_*ui + jdim;

                  estif_(fvi_p_idim,fui_p_jdim) += lin_resM_Du_(nsd_idim_p_jdim,ui)*supg_test(vi);
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
              const int fvi_p_idim = nsd_*vi+idim;

              const int nsd_idim=nsd_*idim;

              for (int ui=0; ui<nen_; ++ui)
              {
                const int fui_p_idim   = nsd_*ui + idim;

                estif_(fvi_p_idim,fui_p_idim) += lin_resM_Du_(nsd_idim+idim,ui)*supg_test(vi);
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
          // TODO
          const double v = timefacfac*supg_test(vi);

          for (int idim = 0; idim <nsd_; ++idim)
          {
            const int fvi   = nsd_*vi + idim;
            for (int ui=0; ui<nen_; ++ui)
            {
              //const int fuippp = numdofpernode_*ui + nsd_;

              estif_p_v_(fvi,ui) += v*derxy_(idim, ui);
            }
          }
        }  // end for(idim)

        /* supg stabilisation: inertia, linearisation of testfunction  */
        /*
                    /                                       \
                   |         n+1       /              \      |
                   |    rho*u      ,  | rho*Du o nabla | v   |
                   |         (i)       \              /      |
                    \                                       /
        */
        /* supg stabilisation: reactive part of convection and linearisation of testfunction ( L_conv_u) */
        /*
                    /                                                       \
                   |    /       n+1        \   n+1     /              \      |
                   |   |   rho*u    o nabla | u    ,  | rho*Du o nabla | v   |
                   |    \       (i)        /   (i)     \              /      |
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
                  |               / n+1 \    /               \      |
                  |  nabla o eps | u     |, |  rho*Du o nabla | v   |
                  |               \ (i) /    \               /      |
                   \                                               /
        */
        /* supg stabilisation: bodyforce part, linearisation of test function */
        /*
                   /                                      \
                  |                  /               \     |
                  |  rho*rhsint   , |  rho*Du o nabla | v  |
                  |                  \               /     |
                   \                                      /
        */
        if (is_newton_)
        {
          for(int jdim=0;jdim<nsd_;++jdim)
          {
            const double temp=timetauM*densaf_*momres_old_(jdim);

            for (int vi=0; vi<nen_; ++vi)
            {
              const int fvi_p_jdim = nsd_*vi+jdim;

              for(int idim=0;idim<nsd_;++idim)
              {
                const double v=temp*derxy_(idim,vi);

                for (int ui=0; ui<nen_; ++ui)
                {
                  const int fui_p_idim   = nsd_*ui + idim;

                  estif_(fvi_p_jdim,fui_p_idim) += v*funct_(ui);
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
          const double temp=momres_old_(idim)*fac_;

          for (int vi=0; vi<nen_; ++vi)
          {
            // supg stabilisation
            velforce(idim,vi) -= temp*supg_test(vi);
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
          const double two_visc_timetauMp      = vstabfac*2.0*visc_*timefac*tau_(1);
          const double two_visc_timetauMp_fac  = vstabfac*2.0*visc_*timefac*tau_(1)*fac_;

          /* viscous stabilisation, inertia part if not stationary */
          /*
                     /                        \
                    |                          |
                +/- |    rho*Du , div eps (v)  |
                    |                          |
                     \                        /
          */
          /* viscous stabilisation, convective part, convective type */
          /*
                   /                                      \
                  |  /       n+1       \                   |
              +/- | |   rho*u   o nabla | Du , div eps (v) |
                  |  \       (i)       /                   |
                   \                                      /
          */
          /* viscous stabilisation, reactive part of convection */
          /*
                   /                                       \
                  |  /                \   n+1               |
              +/- | |   rho*Du o nabla | u    , div eps (v) |
                  |  \                /   (i)               |
                   \                                       /
          */
          /* viscous stabilisation, viscous part (-L_visc_u) */
          /*
                   /                                 \
                  |               /  \                |
             -/+  |  nabla o eps | Du | , div eps (v) |
                  |               \  /                |
                   \                                 /
          */
          for(int jdim=0;jdim<nsd_;++jdim)
          {
            for (int ui=0; ui<nen_; ++ui)
            {
              const int fui_p_jdim   = nsd_*ui + jdim;

              for(int idim=0;idim<nsd_;++idim)
              {
                const int nsd_idim_p_jdim=nsd_*idim+jdim;

                for (int vi=0; vi<nen_; ++vi)
                {
                  const int fvi_p_idim = nsd_*vi+idim;

                  estif_(fvi_p_idim,fui_p_jdim) += two_visc_timetauMp*lin_resM_Du_(nsd_idim_p_jdim,ui)*viscs2_(nsd_idim_p_jdim,vi);
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
              //const int fui = ui*numdofpernode_ + nsd_;
              for (int vi=0; vi<nen_; ++vi)
              {
                //const int fvi = vi*numdofpernode_ + idim;

                for(int jdim=0;jdim<nsd_;++jdim)
                {
                  ppmat(vi, ui) += two_visc_timetauMp_fac*(derxy_(jdim, ui))*viscs2_(jdim+(idim*nsd_), vi);
                  //estif_p_q(fvi, fui) += two_visc_timetauMp_fac*(derxy_(jdim, ui))*viscs2_(jdim+(idim*nsd_), vi);
                }
              }
            }
          } // end for(idim)


          for (int idim =0;idim<nsd_;++idim)
          {
            for (int vi=0; vi<nen_; ++vi)
            {
              /* viscous stabilisation */
              for (int jdim=0;jdim<nsd_;++jdim)
              {
                velforce(idim,vi)-= two_visc_tauMp*momres_old_(jdim)*viscs2_(jdim+(idim*nsd_), vi);
              }
            }
          } // end for(idim)

        } // end if viscous stabilization on left hand side
      } //end if(is_higherOrder)


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
      //} end Galerkin Part

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


  // add pressure part of residual to force vector
  for (int vi=0; vi<nen_; ++vi)
  {
    eforce(numdofpernode_*vi+nsd_)+=preforce(vi);
  }

  // add velocity part of residual to force vector
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int idim=0; idim<nsd_; ++idim)
    {
      eforce(numdofpernode_*vi+idim)+=velforce(idim,vi);
    }
  }

  // add pressure/pressure part of matrix to the element stiffness
  for (int ui=0; ui<nen_; ++ui)
  {
    const int fuippp = numdofpernode_*ui+nsd_;


    for (int vi=0; vi<nen_; ++vi)
    {
      estif(numdofpernode_*vi+nsd_,fuippp)+=ppmat(vi,ui);
    }
  }

  // sort estif_ into the global matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    for (int vi=0; vi<nen_; ++vi)
    {
      for (int idim=0; idim <nsd_; ++idim)
      {
        for (int jdim=0; jdim < nsd_;++jdim)
          estif(numdofpernode_*vi+idim, numdofpernode_*ui+jdim) += estif_(nsd_*vi+idim, nsd_*ui+jdim);
      }
    }
  }

  // sort estif_ into the global matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    for (int vi=0; vi<nen_; ++vi)
    {
      for (int idim=0; idim <nsd_; ++idim)
        estif(numdofpernode_*vi+idim, numdofpernode_*ui+nsd_) += estif_p_v_(nsd_*vi+idim, ui);
    }
  }

  // sort estif_ into the global matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    for (int vi=0; vi<nen_; ++vi)
    {
      for (int jdim=0; jdim < nsd_;++jdim)
        estif(numdofpernode_*vi+nsd_, numdofpernode_*ui+jdim) += estif_q_u_(vi, nsd_*ui+jdim);
    }
  }

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

