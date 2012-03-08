/*!------------------------------------------------------------------------------------------------*
\file topopt_fluidAdjoint3_impl.cpp

\brief 

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET


#include "../drt_f3/fluid3.H"
#include "../drt_f3/fluid3_ele_impl_utils.H"
#include "../drt_fem_general/drt_utils_gder2.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_mat/newtonianfluid.H"

#include "topopt_fluidAdjoint3_impl.H"
#include "topopt_fluidAdjoint3_impl_parameter.H"


//----------------------------------------------------------------------*
//
//----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidAdjoint3ImplInterface* DRT::ELEMENTS::FluidAdjoint3ImplInterface::Impl(DRT::Element::DiscretizationType distype)
{
  switch(distype)
  {
  case DRT::Element::hex8:
  {
    return FluidAdjoint3Impl<DRT::Element::hex8>::Instance();
  }
  case DRT::Element::hex20:
  {
    return FluidAdjoint3Impl<DRT::Element::hex20>::Instance();
  }
  case DRT::Element::hex27:
  {
    return FluidAdjoint3Impl<DRT::Element::hex27>::Instance();
  }
//  case DRT::Element::tet4:
//  {
//    return FluidAdjoint3Impl<DRT::Element::tet4>::Instance();
//  }
//  case DRT::Element::tet10:
//  {
//    return FluidAdjoint3Impl<DRT::Element::tet10>::Instance();
//  }
//  case DRT::Element::wedge6:
//  {
//    return FluidAdjoint3Impl<DRT::Element::wedge6>::Instance();
//  }
//  case DRT::Element::pyramid5:
//  {
//    return FluidAdjoint3Impl<DRT::Element::pyramid5>::Instance();
//  }
  case DRT::Element::quad4:
  {
    return FluidAdjoint3Impl<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return FluidAdjoint3Impl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return FluidAdjoint3Impl<DRT::Element::quad9>::Instance();
  }
//  case DRT::Element::tri3:
//  {
//    return FluidAdjoint3Impl<DRT::Element::tri3>::Instance();
//  }
//  case DRT::Element::tri6:
//  {
//    return FluidAdjoint3Impl<DRT::Element::tri6>::Instance();
//  }
  // no 1D elements
  default:
    dserror("Element shape %s not activated. Just do it.",DRT::DistypeToString(distype).c_str());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidAdjoint3Impl<distype> * DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Instance( bool create )
{
  static FluidAdjoint3Impl<distype> * instance;
  if ( create )
  {
    if ( instance==NULL )
    {
      instance = new FluidAdjoint3Impl<distype>();
    }
  }
  else
  {
    if ( instance!=NULL )
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidAdjoint3Impl<distype>::FluidAdjoint3Impl()
  : xyze_(true),
    funct_(true),
    deriv_(true),
    deriv2_(true),
    xjm_(true),
    xji_(true),
    derxy_(true),
    derxy2_(true),
    velint_(true),
    vderxy_(true),
    gradp_(true),
    fluidvelint_(true),
    fluidvelxy_(true),
    vdiv_(0.0),
    conv1_(true),
    conv2_(true),
    velint_old_(true),
    vderxy_old_(true),
    gradp_old_(true),
    fluidvelint_old_(true),
    fluidvelxy_old_(true),
    vdiv_old_(true),
    conv1_old_(true),
    conv2_old_(true),
    tau_(true),
    intpoints_( distype ),
    xsi_(true),
    det_(0.0),
    fac_(0.0),
    visc_(0.0),
    reacoeff_(0.0),
    dens_(0.0),
    is_higher_order_ele_(false)
{
  // pointer to class Fluid3ImplParameter (access to the general parameter)
  fluidAdjoint3Parameter_ = DRT::ELEMENTS::FluidAdjoint3ImplParameter::Instance();
}



/*----------------------------------------------------------------------*
 * Action type: Compute Error                                 ehrl 02/11
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ComputeError(
    DRT::ELEMENTS::Fluid3*          ele,
    ParameterList&                  params,
    Teuchos::RCP<MAT::Material>&    mat,
    DRT::Discretization&            discretization,
    vector<int>&                    lm,
    Epetra_SerialDenseVector&       elevec1
    )
{
  // analytical solution
  LINALG::Matrix<nsd_,1>  u(true);
  double p = 0.0;

  // error
  LINALG::Matrix<nsd_,1> deltavel(true);
  double         deltap=0.0;

  const int calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(params,"calculate error");

  //----------------------------------------------------------------------------
  //   Extract velocity/pressure from global vectors
  //----------------------------------------------------------------------------

  // fill the local element vector/matrix with the global values
  LINALG::Matrix<nsd_,nen_> evelaf(true);
  LINALG::Matrix<nen_,1> epreaf(true);
  ExtractValuesFromGlobalVector(discretization,lm, &evelaf, &epreaf,"u and p at time n+1 (converged)");

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  if (ele->IsAle())
  {
    LINALG::Matrix<nsd_,nen_>       edispnp(true);
    ExtractValuesFromGlobalVector(discretization,lm, &edispnp, NULL,"dispnp");

    // get new node positions for isale
     xyze_ += edispnp;
  }

  // integrations points and weights
  // more GP than usual due to (possible) cos/exp fcts in analytical solutions
  // degree 5
  DRT::UTILS::GaussIntegration intpoints(distype, 5);

//------------------------------------------------------------------
//                       INTEGRATION LOOP
//------------------------------------------------------------------

  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad,ele->Id());

    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    velint_.Multiply(evelaf,funct_);

    // get pressure at integration point
    // (value at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    double preint = funct_.Dot(epreaf);

    // get coordinates at integration point
    LINALG::Matrix<nsd_,1> xyzint(true);
    xyzint.Multiply(xyze_,funct_);

    // Compute analytical solution
    switch(calcerr)
    {
    case INPAR::FLUID::beltrami_flow:
    {
      if (nsd_ == 3)
      {
         // get viscosity
        if (mat->MaterialType() == INPAR::MAT::m_fluid)
        {
          const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());

          // get constant kinematic viscosity
          visc_ = actmat->Viscosity()/actmat->Density();
        }
        else dserror("Material is not Newtonian Fluid");

         const double a      = M_PI/4.0;
         const double d      = M_PI/2.0;

         const double t = fluidAdjoint3Parameter_->time_;

         // compute analytical pressure
         p = -a*a/2.0 *
             ( exp(2.0*a*xyzint(0))
             + exp(2.0*a*xyzint(1))
             + exp(2.0*a*xyzint(2))
             + 2.0 * sin(a*xyzint(0) + d*xyzint(1)) * cos(a*xyzint(2) + d*xyzint(0)) * exp(a*(xyzint(1)+xyzint(2)))
             + 2.0 * sin(a*xyzint(1) + d*xyzint(2)) * cos(a*xyzint(0) + d*xyzint(1)) * exp(a*(xyzint(2)+xyzint(0)))
             + 2.0 * sin(a*xyzint(2) + d*xyzint(0)) * cos(a*xyzint(1) + d*xyzint(2)) * exp(a*(xyzint(0)+xyzint(1)))
             )* exp(-2.0*visc_*d*d*t);

          // compute analytical velocities
          u(0) = -a * ( exp(a*xyzint(0)) * sin(a*xyzint(1) + d*xyzint(2)) +
                       exp(a*xyzint(2)) * cos(a*xyzint(0) + d*xyzint(1)) ) * exp(-visc_*d*d*t);
          u(1) = -a * ( exp(a*xyzint(1)) * sin(a*xyzint(2) + d*xyzint(0)) +
                       exp(a*xyzint(0)) * cos(a*xyzint(1) + d*xyzint(2)) ) * exp(-visc_*d*d*t);
          u(2) = -a * ( exp(a*xyzint(2)) * sin(a*xyzint(0) + d*xyzint(1)) +
                       exp(a*xyzint(1)) * cos(a*xyzint(2) + d*xyzint(0)) ) * exp(-visc_*d*d*t);
        }
        else dserror("action 'calc_fluid_beltrami_error' is a 3D specific action");
      break;
    }
    default:
      dserror("analytical solution is not defined");
    }

    // compute difference between analytical solution and numerical solution
    deltap    = preint - p;
    deltavel.Update(1.0, velint_, -1.0, u);

    // L2 error
    // 0: vel_mag
    // 1: p
    // 2: vel_mag,analytical
    // 3: p_analytic
    // (4: vel_x)
    // (5: vel_y)
    // (6: vel_z)
    for (int isd=0;isd<nsd_;isd++)
    {
      elevec1[0] += deltavel(isd)*deltavel(isd)*fac_;
      //integrate analytical velocity (computation of relative error)
      elevec1[2] += u(isd)*u(isd)*fac_;
      // velocity components
      //elevec1[isd+4] += deltavel(isd)*deltavel(isd)*fac_;
    }
    elevec1[1] += deltap*deltap*fac_;
    //integrate analytical pressure (computation of relative error)
    elevec1[3] += p*p*fac_;
  }

  return 0;
}



/*----------------------------------------------------------------------*
 * Action type: Evaluate
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Evaluate(DRT::ELEMENTS::Fluid3*    ele,
                                                 DRT::Discretization & discretization,
                                                 const std::vector<int> & lm,
                                                 Teuchos::ParameterList&    params,
                                                 Teuchos::RCP<MAT::Material> & mat,
                                                 Epetra_SerialDenseMatrix&  elesysmat,
                                                 Epetra_SerialDenseMatrix&  elematdummy,
                                                 Epetra_SerialDenseVector&  elerhs,
                                                 Epetra_SerialDenseVector&  elevecdummy1,
                                                 Epetra_SerialDenseVector&  elevecdummy2
)
{
  return Evaluate( ele, discretization, lm, params, mat,
                   elesysmat, elerhs, intpoints_ );
}



template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Evaluate(DRT::ELEMENTS::Fluid3*    ele,
                                                 DRT::Discretization & discretization,
                                                 const std::vector<int> & lm,
                                                 Teuchos::ParameterList&    params,
                                                 Teuchos::RCP<MAT::Material> & mat,
                                                 Epetra_SerialDenseMatrix&  elesysmat,
                                                 Epetra_SerialDenseVector&  elerhs,
                                                 const DRT::UTILS::GaussIntegration & intpoints
)
{
  // construct views
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> elemat(elesysmat,true);
  LINALG::Matrix<(nsd_+1)*nen_,            1> elevec(elerhs,true);
  // elevec2 and elevec3 are currently not in use

  // ---------------------------------------------------------------------
  // get all general state vectors: fluid/adjoint velocity/pressure
  // velocity/pressure values are at time n/n+1
  // ---------------------------------------------------------------------
  // fill the local element vector/matrix with the global values
  // ost:         velocity/pressure at time n+1
  LINALG::Matrix<nsd_,nen_> eveln(true);
  LINALG::Matrix<nen_,1>    epren(true);
  ExtractValuesFromGlobalVector(discretization,lm, &eveln, &epren,"veln");

  LINALG::Matrix<nsd_,nen_> evelnp(true);
  LINALG::Matrix<nen_,1>    eprenp(true);
  ExtractValuesFromGlobalVector(discretization,lm, &evelnp, &eprenp,"velnp");

  LINALG::Matrix<nsd_,nen_> efluidveln(true);
  ExtractValuesFromGlobalVector(discretization,lm, &efluidveln, NULL,"fluidveln");

  LINALG::Matrix<nsd_,nen_> efluidvelnp(true);
  ExtractValuesFromGlobalVector(discretization,lm, &efluidvelnp, NULL,"fluidvelnp");

  // evaluate nodal porosities
  LINALG::Matrix<nen_,1> eporo(true);
  {
    // read nodal values from global vector
    RCP<Epetra_Vector> topopt_porosity = params.get<RCP<Epetra_Vector> >("topopt_porosity");
    for (int nn=0;nn<nen_;++nn)
    {
      int lid = (ele->Nodes()[nn])->LID();
      eporo(nn,0) = (*topopt_porosity)[lid];
    }
  }

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // call inner evaluate (does not know about DRT element or discretization object)
  int result = Evaluate(
    ele->Id(),
    params,
    elemat,
    elevec,
    eveln,
    evelnp,
    epren,
    eprenp,
    efluidveln,
    efluidvelnp,
    eporo,
    mat,
    intpoints);

  return result;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Evaluate(
  int                                           eid,
  Teuchos::ParameterList&                       params,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_> & elesysmat,
  LINALG::Matrix<(nsd_+1)*nen_,            1> & elerhs,
  const LINALG::Matrix<nsd_,nen_> &             eveln,
  const LINALG::Matrix<nsd_,nen_> &             evelnp,
  const LINALG::Matrix<nen_,1>    &             epren,
  const LINALG::Matrix<nen_,1>    &             eprenp,
  const LINALG::Matrix<nsd_,nen_> &             efluidveln,
  const LINALG::Matrix<nsd_,nen_> &             efluidvelnp,
  const LINALG::Matrix<nen_,1> &                eporo,
  Teuchos::RCP<MAT::Material>                   mat,
  const DRT::UTILS::GaussIntegration &          intpoints )
{
  // flag for higher order elements
  is_higher_order_ele_ = IsHigherOrder<distype>::ishigherorder;
  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if (fluidAdjoint3Parameter_->is_inconsistent_ == true) is_higher_order_ele_ = false;

  // ---------------------------------------------------------------------
  // call routine for calculating element matrix and right hand side
  // ---------------------------------------------------------------------
  Sysmat(eid,
         eveln,
         evelnp,
         epren,
         eprenp,
         efluidveln,
         efluidvelnp,
         elesysmat,
         elerhs,
         eporo,
         mat,
         intpoints);

  return 0;
}


/*----------------------------------------------------------------------*
 |  calculate element matrix and right hand side (private)   g.bau 03/07|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::Sysmat(
  int                                           eid,
  const LINALG::Matrix<nsd_,nen_>&              eveln,
  const LINALG::Matrix<nsd_,nen_>&              evelnp,
  const LINALG::Matrix<nen_,1>&                 epren,
  const LINALG::Matrix<nen_,1>&                 eprenp,
  const LINALG::Matrix<nsd_,nen_> &             efluidveln,
  const LINALG::Matrix<nsd_,nen_> &             efluidvelnp,
  LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  estif,
  LINALG::Matrix<(nsd_+1)*nen_,1>&              eforce,
  const LINALG::Matrix<nen_,1> &                eporo,
  Teuchos::RCP<const MAT::Material>             material,
  const DRT::UTILS::GaussIntegration &          intpoints
  )
{
  //------------------------------------------------------------------------
  //  preliminary definitions and evaluations
  //------------------------------------------------------------------------
  // definition of matrices
  LINALG::Matrix<nen_*nsd_,nen_*nsd_>  estif_u(true);
  LINALG::Matrix<nen_*nsd_,nen_>       estif_p_v(true);
  LINALG::Matrix<nen_, nen_*nsd_>      estif_q_u(true);
  LINALG::Matrix<nen_,nen_>            ppmat(true);

  // definition of vectors
  LINALG::Matrix<nen_,1>     preforce(true);
  LINALG::Matrix<nsd_,nen_>  velforce(true);

  // definition of velocity-based momentum residual vectors
  LINALG::Matrix<nsd_*nsd_,nen_>  lin_resM_Du(true);
  LINALG::Matrix<nsd_,1>          resM_Du(true);

  // evaluate shape functions and derivatives at element center
  EvalShapeFuncAndDerivsAtEleCenter(eid);

  //------------------------------------------------------------------------
  // potential evaluation of material parameters, subgrid viscosity
  // and/or stabilization parameters at element center
  //------------------------------------------------------------------------
  // get material parameters at element center
  if (not fluidAdjoint3Parameter_->mat_gp_ or not fluidAdjoint3Parameter_->tau_gp_)
    GetMaterialParams(material,eveln);

  // calculate subgrid viscosity and/or stabilization parameter at element center
  if (not fluidAdjoint3Parameter_->tau_gp_)
  {
    // get velocity at element center
    velint_.Multiply(eveln,funct_);

    // calculate stabilization parameters at element center
    CalcStabParameter(fac_);
  }


  //------------------------------------------------------------------------
  //  start loop over integration points
  //------------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad,eid);

    //----------------------------------------------------------------------
    //  evaluation of various values at integration point:
    //  1) velocity (including derivatives)
    //  2) fluid velocity (including derivatives)
    //  3) pressure (including derivatives)
    //  4) body-force vector
    //  5) "history" vector for momentum equation
    //----------------------------------------------------------------------
    // get velocity at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    velint_.Multiply(eveln,funct_);
    velint_old_.Multiply(evelnp,funct_);

    // get velocity derivatives at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    vderxy_.MultiplyNT(eveln,derxy_);
    vderxy_old_.MultiplyNT(evelnp,derxy_);

    // get fluid velocity at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    fluidvelint_.Multiply(efluidveln,funct_);
    fluidvelint_old_.Multiply(efluidvelnp,funct_);

    // get fluid velocity derivatives at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    fluidvelxy_.MultiplyNT(efluidveln,derxy_);
    fluidvelxy_old_.MultiplyNT(efluidvelnp,derxy_);

    // get pressure at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    double press = funct_.Dot(epren);
    double press_old = funct_.Dot(eprenp);

    // get pressure gradient at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    gradp_.Multiply(derxy_,epren);
    gradp_old_.Multiply(derxy_,eprenp);


    //----------------------------------------------------------------------
    // potential evaluation of material parameters, subgrid viscosity
    // and/or stabilization parameters at integration point
    //----------------------------------------------------------------------
    // get material parameters at integration point
    if (fluidAdjoint3Parameter_->mat_gp_)
      GetMaterialParams(material,eveln);

    // get reaction coefficient due to porosity for topology optimization
    // !do this only at gauss point!
    // TODO does it make problems to evaluate at element center? (i think it should, winklmaier)
    reacoeff_ = funct_.Dot(eporo);

    // calculate stabilization parameter at integration point
    if (fluidAdjoint3Parameter_->tau_gp_)
      CalcStabParameter(fac_);

    // get first convective value at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    conv1_.Multiply(vderxy_,fluidvelint_);
    conv1_old_.Multiply(vderxy_old_,fluidvelint_old_);

    // get second convective value at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    conv2_.Multiply(fluidvelxy_,velint_);
    conv2_old_.Multiply(fluidvelxy_old_,velint_old_);

    // get divergence at integration point
    // 1) t^n=last iteration 2) t^n+1 = last time step
    vdiv_ = vdiv_old_ = 0.0;
    for (int idim = 0; idim <nsd_; ++idim)
    {
      vdiv_ += vderxy_(idim,idim);
      vdiv_old_ += vderxy_old_(idim,idim);
    }

    //----------------------------------------------------------------------
    // set time-integration factors for left- and right-hand side
    //----------------------------------------------------------------------
    const double timefacfac       = fluidAdjoint3Parameter_->timefac_       * fac_;
    const double timefacfacrhs    = fluidAdjoint3Parameter_->timefacrhs_    * fac_;

    const double timefacfacpre    = fluidAdjoint3Parameter_->timefacpre_    * fac_;
    const double timefacfacprerhs = fluidAdjoint3Parameter_->timefacprerhs_ * fac_;

    const double timefacfacdiv    = fluidAdjoint3Parameter_->timefacdiv_    * fac_;
    const double timefacfacdivrhs = fluidAdjoint3Parameter_->timefacdivrhs_ * fac_;

    // TODO
    // o boundary conditions
    // o body force


    /* ------------------------------------------------------------------------ *
    * standard terms                                                            *
    * ------------------------------------------------------------------------- */

    // 1) mass matrix + reactive term
    MassReactionGalPart(estif_u,velforce,timefacfac,timefacfacrhs);

    // 2) convection terms
    ConvectionGalPart(estif_u,velforce,timefacfac,timefacfacrhs);

    // 3) viscous terms
    ViscousGalPart(estif_u,velforce,timefacfac,timefacfacrhs);

    // 4) pressure term
    PressureGalPart(estif_p_v,velforce,timefacfacpre,timefacfacprerhs,press,press_old);

    // 5) continuity term
    ContinuityGalPart(estif_q_u,preforce,timefacfacdiv,timefacfacdivrhs);

    // 6) standard Galerkin bodyforce term on right-hand side
    BodyForceGalPart(velforce,timefacfac,timefacfacrhs);

    /* ------------------------------------------------------------------------ *
    * standard terms done                                                       *
    * ------------------------------------------------------------------------- */



    /* ------------------------------------------------------------------------ *
    * stabilization part                                                        *
    * ------------------------------------------------------------------------- */

    if ((fluidAdjoint3Parameter_->pspg_ == INPAR::FLUID::pstab_use_pspg) or
        (fluidAdjoint3Parameter_->supg_ == INPAR::FLUID::convective_stab_supg))
    {
      // prework for supg/psgp - stabilization: evaluate strong residual

      /* order of the derivatives in GalMomResnU is:
       * from 1 to nsd:       col-dim = x, row-dim = 1-nsd
       * from nsd+1 to 2*nsd: col-dim = y, row-dim = 1-nsd
       * and so on. so the outer loop is the column dimension
       * and the inner loop the row dimension */
      LINALG::Matrix<nsd_*nsd_,nen_> GalMomResnU(true);

      // strong residual of momentum equation of last iteration, scaled with fac*dt/rho
      LINALG::Matrix<nsd_,1> StrongResMomScaled(true);

      MomRes(GalMomResnU,
          StrongResMomScaled,
          timefacfac,
          timefacfacrhs,
          timefacfacpre,
          timefacfacprerhs,
          eveln,
          evelnp,
          efluidveln,
          efluidvelnp);

      // 7) PSPG term
      if (fluidAdjoint3Parameter_->pspg_ == INPAR::FLUID::pstab_use_pspg)
      {
        PSPG(estif_q_u,
            ppmat,
            preforce,
            GalMomResnU,
            StrongResMomScaled,
            timefacfac,
            timefacfacrhs,
            timefacfacpre,
            timefacfacprerhs);
      }

      // 8) SUPG term
      if (fluidAdjoint3Parameter_->supg_ == INPAR::FLUID::convective_stab_supg)
      {
        SUPG(estif_u,
            estif_p_v,
            velforce,
            GalMomResnU,
            StrongResMomScaled,
            timefacfac,
            timefacfacrhs,
            timefacfacpre,
            timefacfacprerhs);
      }
    }

    // 9) continuity stabilization
    if (fluidAdjoint3Parameter_->cstab_ == INPAR::FLUID::continuity_stab_yes)
    {
      ContStab(estif_u,
          velforce,
          timefacfacdiv,
          timefacfacdivrhs);
    }
  }
  //------------------------------------------------------------------------
  //  end loop over integration points
  //------------------------------------------------------------------------

  //------------------------------------------------------------------------
  //  add contributions to element matrix and right-hand-side vector
  //------------------------------------------------------------------------
  // add pressure part to right-hand-side vector
  for (int vi=0; vi<nen_; ++vi)
  {
    eforce(numdofpernode_*vi+nsd_)+=preforce(vi);
  }

  // add velocity part to right-hand-side vector
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int idim=0; idim<nsd_; ++idim)
    {
      eforce(numdofpernode_*vi+idim)+=velforce(idim,vi);
    }
  }

  // add pressure-pressure part to matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int fuippp = numdofpernode_*ui+nsd_;

    for (int vi=0; vi<nen_; ++vi)
    {
      const int numdof_vi_p_nsd = numdofpernode_*vi+nsd_;

      estif(numdof_vi_p_nsd,fuippp)+=ppmat(vi,ui);
    }
  }

  // add velocity-velocity part to matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int numdof_ui = numdofpernode_*ui;
    const int nsd_ui = nsd_*ui;

    for (int jdim=0; jdim < nsd_;++jdim)
    {
      const int numdof_ui_jdim = numdof_ui+jdim;
      const int nsd_ui_jdim = nsd_ui+jdim;

      for (int vi=0; vi<nen_; ++vi)
      {
        const int numdof_vi = numdofpernode_*vi;
        const int nsd_vi = nsd_*vi;

        for (int idim=0; idim <nsd_; ++idim)
        {
          estif(numdof_vi+idim, numdof_ui_jdim) += estif_u(nsd_vi+idim, nsd_ui_jdim);
        }
      }
    }
  }

  // add velocity-pressure part to matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int numdof_ui_nsd = numdofpernode_*ui + nsd_;

    for (int vi=0; vi<nen_; ++vi)
    {
      const int nsd_vi = nsd_*vi;
      const int numdof_vi = numdofpernode_*vi;

      for (int idim=0; idim <nsd_; ++idim)
      {
        estif(numdof_vi+idim, numdof_ui_nsd) += estif_p_v(nsd_vi+idim, ui);
      }
    }
  }

  // add pressure-velocity part to matrix
  for (int ui=0; ui<nen_; ++ui)
  {
    const int numdof_ui = numdofpernode_*ui;
    const int nsd_ui = nsd_*ui;

    for (int jdim=0; jdim < nsd_;++jdim)
    {
      const int numdof_ui_jdim = numdof_ui+jdim;
      const int nsd_ui_jdim = nsd_ui+jdim;

      for (int vi=0; vi<nen_; ++vi)
        estif(numdofpernode_*vi+nsd_, numdof_ui_jdim) += estif_q_u(vi, nsd_ui_jdim);
    }
  }

  return;
}



/*----------------------------------------------------------------------*
 |  compute body force at element nodes (private)              vg 10/11 |
 *----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at element center  vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::EvalShapeFuncAndDerivsAtEleCenter(
  const int  eleid
)
{
  // use one-point Gauss rule
  DRT::UTILS::IntPointsAndWeights<nsd_> intpoints_stab(DRT::ELEMENTS::DisTypeToStabGaussRule<distype>::rule);

  // coordinates of the current integration point
  const double* gpcoord = (intpoints_stab.IP().qxg)[0];
  for (int idim=0;idim<nsd_;idim++)
  {
    xsi_(idim) = gpcoord[idim];
  }
  const double wquad = intpoints_stab.IP().qwgt[0];

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
  if (is_higher_order_ele_)
  {
    // get the second derivatives of standard element at current GP
    DRT::UTILS::shape_function_deriv2<distype>(xsi_,deriv2_);
  }


  // compute Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
   */

  // get Jacobian matrix and determinant
  xjm_.MultiplyNT(deriv_,xyze_);
  det_ = xji_.Invert(xjm_);

  // check for degenerated elements
  if (det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det_);

  // compute integration factor
  fac_ = wquad*det_;

  // compute global first derivates
  derxy_.Multiply(xji_,deriv_);

  //--------------------------------------------------------------
  //             compute global second derivatives
  //--------------------------------------------------------------
  if (is_higher_order_ele_)
  {
    DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else derxy2_.Clear();

  return;
}



/*----------------------------------------------------------------------*
 | evaluate shape functions and derivatives at integr. point   vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::EvalShapeFuncAndDerivsAtIntPoint(
    DRT::UTILS::GaussIntegration::iterator & iquad,       // actual integration point
    const int                              eleid        // element ID
)
{
  // coordinates of the current integration point
  const double* gpcoord = iquad.Point();
  for (int idim=0;idim<nsd_;idim++)
  {
     xsi_(idim) = gpcoord[idim];
  }

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);
  derxy2_.Clear();
  if (is_higher_order_ele_)
  {
    // get the second derivatives of standard element at current GP
    DRT::UTILS::shape_function_deriv2<distype>(xsi_,deriv2_);
  }

  // get Jacobian matrix and determinant
  // actually compute its transpose....
  /*
    +-            -+ T      +-            -+
    | dx   dx   dx |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dr   dr   dr |
    |              |        |              |
    | dy   dy   dy |        | dx   dy   dz |
    | --   --   -- |   =    | --   --   -- |
    | dr   ds   dt |        | ds   ds   ds |
    |              |        |              |
    | dz   dz   dz |        | dx   dy   dz |
    | --   --   -- |        | --   --   -- |
    | dr   ds   dt |        | dt   dt   dt |
    +-            -+        +-            -+
  */
  xjm_.MultiplyNT(deriv_,xyze_);
  det_ = xji_.Invert(xjm_);

  if (det_ < 1E-16)
    dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eleid, det_);

  // compute integration factor
  fac_ = iquad.Weight()*det_;

  // compute global first derivates
  derxy_.Multiply(xji_,deriv_);

  //--------------------------------------------------------------
  //             compute global second derivatives
  //--------------------------------------------------------------
  if (is_higher_order_ele_)
  {
    DRT::UTILS::gder2<distype>(xjm_,derxy_,deriv2_,xyze_,derxy2_);
  }
  else derxy2_.Clear();

  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::GetMaterialParams(
  Teuchos::RCP<const MAT::Material>  material,
  const LINALG::Matrix<nsd_,nen_>&   evelaf
)
{

if (material->MaterialType() == INPAR::MAT::m_fluid)
{
  const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());

  // get constant density
  dens_ = actmat->Density();

  // get constant dynamic viscosity
  visc_ = actmat->Viscosity();
}
else
  dserror("Material type is not supported");

// check whether there is zero or negative (physical) viscosity
// (expect for permeable fluid)
if (visc_ < EPS15)
  dserror("zero or negative (physical) diffusivity");

return;
}



/*----------------------------------------------------------------------*
 |  calculation of stabilization parameter                     vg 09/09 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::CalcStabParameter(const double vol)
{
  //---------------------------------------------------------------------
  // preliminary definition of values which will already be computed for
  // tau_M and later be used for tau_C again by some of the subsequent
  // stabilization parameter definitions
  //---------------------------------------------------------------------
  double traceG = 0.0;
  double Gnormu = 0.0;
  double Gvisc  = 0.0;

  double strle    = 0.0;
  double hk       = 0.0;
  double fluidvel_norm = 0.0;
  double re12     = 0.0;
  double c3       = 0.0;

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
  // E) definition according to Franca et al. (2005) as well as Badia
  //    and Codina (2010)
  //    -> only for Darcy or Darcy-Stokes/Brinkman flow, hence only
  //       tau_Mp for this definition
  //---------------------------------------------------------------------
  // get element-type constant for tau
  const double mk = DRT::ELEMENTS::MK<distype>();

  // computation depending on which parameter definition is used
  switch (fluidAdjoint3Parameter_->whichtau_)
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

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient (reaction coefficient
    // ensured to remain zero in GetMaterialParams for non-reactive material)
    double sigma_tot = reacoeff_;
    if (fluidAdjoint3Parameter_->whichtau_ == INPAR::FLUID::tau_taylor_hughes_zarins or
        fluidAdjoint3Parameter_->whichtau_ == INPAR::FLUID::tau_taylor_hughes_zarins_whiting_jansen or
        fluidAdjoint3Parameter_->whichtau_ == INPAR::FLUID::tau_taylor_hughes_zarins_scaled)
      sigma_tot += 1.0/fluidAdjoint3Parameter_->dt_;

    // definition of constants as described above
    const double c1 = 4.0;
    c3 = 12.0/mk;

    // computation of various values derived from covariant metric tensor
    // (trace of covariant metric tensor required for computation of tau_C below)
    double G;
    double normG = 0.0;
    const double dens_sqr = dens_*dens_;
    for (int nn=0;nn<nsd_;++nn)
    {
      const double dens_sqr_velint_nn = dens_sqr*fluidvelint_(nn);
      for (int mm=0; mm<nsd_; ++mm)
      {
        traceG += xji_(nn,mm)*xji_(nn,mm);
      }
      for (int rr=0;rr<nsd_;++rr)
      {
        G = xji_(nn,0)*xji_(rr,0);
        for (int mm=1; mm<nsd_; ++mm)
        {
          G += xji_(nn,mm)*xji_(rr,mm);
        }
        normG  += G*G;
        Gnormu += dens_sqr_velint_nn*G*fluidvelint_(rr);
      }
    }

    // compute viscous part
    Gvisc = c3*visc_*visc_*normG;

    // computation of stabilization parameters tau_Mu and tau_Mp
    // -> identical for the present definitions
    tau_(0) = 1.0/(sqrt(c1*dens_sqr*DSQR(sigma_tot) + Gnormu + Gvisc));
    tau_(1) = tau_(0);
    break;
  }

  case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall:
  {
    /*

    literature:
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

    */
    // get velocity norm
    fluidvel_norm = fluidvelint_.Norm2();

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient (reaction coefficient
    // ensured to remain zero in GetMaterialParams for non-reactive material)
    const double sigma_tot = 1.0/fluidAdjoint3Parameter_->timefac_ + reacoeff_;

    // calculate characteristic element length
    CalcCharEleLength(vol,fluidvel_norm,strle,hk);

    // various parameter computations for case with dt:
    // relating viscous to reactive part (re01: tau_Mu, re11: tau_Mp)
    const double re01 = 4.0 * visc_ / (mk * dens_ * sigma_tot * DSQR(strle));
    const double re11 = 4.0 * visc_ / (mk * dens_ * sigma_tot * DSQR(hk));

    // relating convective to viscous part (re02: tau_Mu, re12: tau_Mp)
    const double re02 = mk * dens_ * fluidvel_norm * strle / (2.0 * visc_);
                 re12 = mk * dens_ * fluidvel_norm * hk / (2.0 * visc_);

    // respective "switching" parameters
    const double xi01 = DMAX(re01,1.0);
    const double xi11 = DMAX(re11,1.0);
    const double xi02 = DMAX(re02,1.0);
    const double xi12 = DMAX(re12,1.0);

    tau_(0) = DSQR(strle)/(DSQR(strle)*dens_*sigma_tot*xi01+(4.0*visc_/mk)*xi02);
    tau_(1) = DSQR(hk)/(DSQR(hk)*dens_*sigma_tot*xi11+(4.0*visc_/mk)*xi12);
    break;
  }

  case INPAR::FLUID::tau_franca_barrenechea_valentin_frey_wall_wo_dt:
  {
    /*

     stabilization parameter as above without inclusion of dt-part

    */
    // get velocity norm
    fluidvel_norm = fluidvelint_.Norm2();

    // calculate characteristic element length
    CalcCharEleLength(vol,fluidvel_norm,strle,hk);

    // various parameter computations for case without dt:
    // relating viscous to reactive part (re01: tau_Mu, re11: tau_Mp)
    double re01 = 4.0 * visc_ / (mk * dens_ * reacoeff_ * DSQR(strle));
    double re11 = 4.0 * visc_ / (mk * dens_ * reacoeff_ * DSQR(hk));

    // relating convective to viscous part (re02: tau_Mu, re12: tau_Mp)
    const double re02 = mk * dens_ * fluidvel_norm * strle / (2.0 * visc_);
                 re12 = mk * dens_ * fluidvel_norm * hk / (2.0 * visc_);

    // respective "switching" parameters
    const double xi01 = DMAX(re01,1.0);
    const double xi11 = DMAX(re11,1.0);
    const double xi02 = DMAX(re02,1.0);
    const double xi12 = DMAX(re12,1.0);

    tau_(0) = DSQR(strle)/(DSQR(strle)*dens_*reacoeff_*xi01+(4.0*visc_/mk)*xi02);
    tau_(1) = DSQR(hk)/(DSQR(hk)*dens_*reacoeff_*xi11+(4.0*visc_/mk)*xi12);
    break;
  }

  case INPAR::FLUID::tau_shakib_hughes_codina:
  case INPAR::FLUID::tau_shakib_hughes_codina_wo_dt:
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

       constants defined as in Shakib (1989) / Shakib and Hughes (1991),
       merely slightly different with respect to c_3:

       c_1 = 4.0 (for version with dt), 0.0 (for version without dt),
       c_2 = 4.0,
       c_3 = 4.0/(m_k*m_k) (36.0 for linear, 576.0 for quadratic ele.)

       Codina (2002) proposed present version without dt and explicit
       definition of constants
       (condition for constants as defined here: c_2 <= sqrt(c_3)).

    */
    // get velocity norm
    fluidvel_norm = fluidvelint_.Norm2();

    // calculate characteristic element length
    CalcCharEleLength(vol,fluidvel_norm,strle,hk);

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient (reaction coefficient
    // ensured to remain zero in GetMaterialParams for non-reactive material)
    double sigma_tot = reacoeff_;
    if (fluidAdjoint3Parameter_->whichtau_ == INPAR::FLUID::tau_shakib_hughes_codina)
      sigma_tot += 1.0/fluidAdjoint3Parameter_->dt_;

    // definition of constants as described above
    const double c1 = 4.0;
    const double c2 = 4.0;
    c3 = 4.0/(mk*mk);
    // alternative value as proposed in Shakib (1989): c3 = 16.0/(mk*mk);

    tau_(0) = 1.0/(sqrt(c1*DSQR(dens_)*DSQR(sigma_tot)
                      + c2*DSQR(dens_)*DSQR(fluidvel_norm)/DSQR(strle)
                      + c3*DSQR(visc_)/(DSQR(strle)*DSQR(strle))));
    tau_(1) = 1.0/(sqrt(c1*DSQR(dens_)*DSQR(sigma_tot)
                      + c2*DSQR(dens_)*DSQR(fluidvel_norm)/DSQR(hk)
                      + c3*DSQR(visc_)/(DSQR(hk)*DSQR(hk))));
    break;
  }
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
    // get velocity norm
    fluidvel_norm = fluidvelint_.Norm2();

    // calculate characteristic element length
    CalcCharEleLength(vol,fluidvel_norm,strle,hk);

    // total reaction coefficient sigma_tot: sum of "artificial" reaction
    // due to time factor and reaction coefficient (reaction coefficient
    // ensured to remain zero in GetMaterialParams for non-reactive material)
    double sigma_tot = reacoeff_;
    if (fluidAdjoint3Parameter_->whichtau_ == INPAR::FLUID::tau_codina)
      sigma_tot += 1.0/fluidAdjoint3Parameter_->dt_;

    // definition of constants as described above
    const double c1 = 1.0;
    const double c2 = 2.0;
    c3 = 4.0/mk;

    tau_(0) = 1.0/(sqrt(c1*dens_*sigma_tot
                      + c2*dens_*fluidvel_norm/strle
                      + c3*visc_/DSQR(strle)));
    tau_(1) = 1.0/(sqrt(c1*dens_*sigma_tot
                      + c2*dens_*fluidvel_norm/hk
                      + c3*visc_/DSQR(hk)));
    break;
  }
  default:
    dserror("unknown definition for tau_M\n %i  ", fluidAdjoint3Parameter_->whichtau_);
  }  // end switch (fluidAdjoint3Parameter_->whichtau_)


  //---------------------------------------------------------------------
  // second step: computation of tau_C with the following options:
  // A) definition according to Taylor et al. (1998)
  // B) definition according to Whiting (1999)/Whiting and Jansen (2001)
  // C) scaled version of definition according to Taylor et al. (1998)
  // D) definition according to Wall (1999)
  // E) definition according to Codina (2002)
  // F) definition according to Badia and Codina (2010)
  //    (only for Darcy or Darcy-Stokes/Brinkman flow)
  //---------------------------------------------------------------------
  // computation depending on which parameter definition is used
  switch (fluidAdjoint3Parameter_->whichtau_)
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

    tau_(2) = sqrt(Gnormu)/traceG;
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

    tau_(2) = 1.0/(tau_(0)*traceG);
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
    const double xi_tau_c = DMIN(reG,1.0);

    tau_(2) = xi_tau_c*sqrt(Gnormu)/traceG;
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
    const double xi_tau_c = DMIN(re12,1.0);

    tau_(2) = 0.5 * dens_ * fluidvel_norm * hk * xi_tau_c;
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

    tau_(2) = DSQR(hk)/(sqrt(c3)*tau_(1));
  }
  break;

  default: dserror("unknown definition for tau_C\n %i  ", fluidAdjoint3Parameter_->whichtau_);
  }  // end switch (fluidAdjoint3Parameter_->whichtau_)

  return;
}


/*----------------------------------------------------------------------*
 |  calculation of characteristic element length               vg 01/11 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::CalcCharEleLength(
    const double  vol,
    const double  fluidvel_norm,
    double&       strle,
    double&       hk
    )
{
  // cast dimension to a double varibale -> pow()
  const double dim = double (nsd_);

  //! direction of flow (normed velocity vector)
  LINALG::Matrix<nsd_,1> fluidvelino;

  //---------------------------------------------------------------------
  // various definitions for characteristic element length for tau_Mu
  //---------------------------------------------------------------------
  // a) streamlength due to Tezduyar et al. (1992) -> default
  // normed velocity vector
  if (fluidvel_norm>=1e-6) fluidvelino.Update(1.0/fluidvel_norm,fluidvelint_);
  else
  {
    fluidvelino.Clear();
    fluidvelino(0,0) = 1.0;
  }

  LINALG::Matrix<nen_,1> tmp;
  tmp.MultiplyTN(derxy_,fluidvelino);
  const double val = tmp.Norm1();
  strle = 2.0/val;

  // b) volume-equivalent diameter (warning: 3-D formula!)
  //strle = pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

  // c) cubic/square root of element volume/area
  //strle = pow(vol,1/dim);

  //---------------------------------------------------------------------
  // various definitions for characteristic element length for tau_Mp
  //---------------------------------------------------------------------
  // a) volume-equivalent diameter -> default for 3-D computations
  if (nsd_==3) hk = pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);

  // b) square root of element area -> default for 2-D computations,
  // may also alternatively be used for 3-D computations
  else if (nsd_==2) hk = pow(vol,1/dim);
  // check for potential 1-D computations
  else dserror("element length calculation not implemented for 1-D computation!");

  return;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::MassReactionGalPart(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_u,
    LINALG::Matrix<nsd_,nen_> &             velforce,
    const double &                          timefacfac,
    const double &                          timefacfacrhs
)
{
  double massreacfac = 0.0; // factor summing up coefficients of reactive term and mass-matrix
  if (fluidAdjoint3Parameter_->is_stationary_)
    massreacfac = reacoeff_*timefacfac/dens_;
  else
    massreacfac = fac_+reacoeff_*timefacfac/dens_; // fac -> mass matrix // reac*timefacfac/dens -> reactive

  for (int ui=0; ui<nen_; ++ui)
  {
    const int fui   = nsd_*ui;

    const double uifunct = massreacfac*funct_(ui);

    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi   = nsd_*vi;

      const double value = funct_(vi)*uifunct;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        estif_u(fvi+idim,fui+idim) += value;
      } // end for (idim)
    } //vi
  } // ui

  // rhs at new time step
  LINALG::Matrix<nsd_,1> scaled_vel(true);
  scaled_vel.Update(massreacfac,velint_);
  for (int vi=0;vi<nen_;++vi)
  {
    for (int jdim=0;jdim<nsd_;++jdim)
    {
      velforce(jdim,vi)-=funct_(vi)*scaled_vel(jdim);
    }
  }

  // rhs at old time step
  if (not fluidAdjoint3Parameter_->is_stationary_)
  {
    double massreacfacrhs = fac_+reacoeff_*timefacfacrhs/dens_; // fac -> mass matrix // reac*timefacfac/dens -> reactive
    scaled_vel.Update(massreacfacrhs,velint_old_);

    for (int vi=0;vi<nen_;++vi)
    {
      for (int jdim=0;jdim<nsd_;++jdim)
      {
        velforce(jdim,vi)-=funct_(vi)*scaled_vel(jdim);
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
  /*  reaction */
  /*
            /                \
           |                  |
           |    sigma*Du , v  |
           |                  |
            \                /
  */
//  if (fluidAdjoint3Parameter_->is_newton_)
//  {
//    for (int ui=0; ui<nen_; ++ui)
//    {
//      const int fui   = nsd_*ui;
//
//      for (int idim = 0; idim <nsd_; ++idim)
//      {
//        const int idim_nsd=idim*nsd_;
//
//        for (int vi=0; vi<nen_; ++vi)
//        {
//          const int fvi   = nsd_*vi;
//
//          const int fvi_p_idim = fvi+idim;
//
//          for (int jdim= 0; jdim<nsd_;++jdim)
//          {
//            estif_u(fvi_p_idim,fui+jdim) += funct_(vi)*lin_resM_Du(idim_nsd+jdim,ui);
//          } // end for (jdim)
//        } // end for (idim)
//      } //vi
//    } // ui
//  }
//  else
//  {
//    for (int ui=0; ui<nen_; ++ui)
//    {
//      const int fui   = nsd_*ui;
//
//      for (int vi=0; vi<nen_; ++vi)
//      {
//        const int fvi   = nsd_*vi;
//
//        for (int idim = 0; idim <nsd_; ++idim)
//        {
//          estif_u(fvi+idim,fui+idim) += funct_(vi)*lin_resM_Du(idim*nsd_+idim,ui);
//        } // end for (idim)
//      } //vi
//    } // ui
//  }
//
//  // inertia terms on the right hand side for instationary fluids
//  if (not fluidAdjoint3Parameter_->is_stationary_)
//  {
//    for (int idim = 0; idim <nsd_; ++idim)
//    {
//      resM_Du(idim)+=fac_*dens_*velint_(idim);
//    }
//  }  // end if (not stationary)
//
////  for (int idim = 0; idim <nsd_; ++idim)
////  {
////    resM_Du(idim)+=rhsfac*dens_*conv_old_(idim);
////  }  // end for(idim)
//
//  for (int idim = 0; idim <nsd_; ++idim)
//  {
//    resM_Du(idim) += rhsfac*reacoeff_*velint_(idim);
//  }
//
//  for (int vi=0; vi<nen_; ++vi)
//  {
//    for(int idim = 0; idim <nsd_; ++idim)
//    {
//      velforce(idim,vi)-=resM_Du(idim)*funct_(vi);
//    }
//  }
  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ConvectionGalPart(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_u,
    LINALG::Matrix<nsd_,nen_> &             velforce,
    const double &                          timefacfac,
    const double &                          timefacfacrhs
)
{
  double value = 0.0; // helper

  for (int vi=0; vi<nen_; ++vi)
  {
    const int fvi   = nsd_*vi;

    value = 0.0;

    for (int dim=0;dim<nsd_;++dim)
      value += derxy_(dim,vi)*fluidvelint_(dim);

    value *= timefacfac;

    for (int ui=0; ui<nen_; ++ui)
    {
      const int fui   = nsd_*ui;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        estif_u(fvi+idim,fui+idim) += value*funct_(ui);
      } // end for (idim)
    } //ui
  } // vi

  // rhs at new time step
  for (int vi=0; vi<nen_; ++vi)
  {
    value = 0.0;

    for (int dim=0;dim<nsd_;++dim)
      value += derxy_(dim,vi)*fluidvelint_(dim);

    value *= timefacfac;

    for (int jdim=0;jdim<nsd_;++jdim)
      velforce(jdim,vi) -= value*velint_(jdim);
  } // vi

  // rhs at old time step
  if (not fluidAdjoint3Parameter_->is_stationary_)
  {
    for (int vi=0; vi<nen_; ++vi)
    {
      value = 0.0;

      for (int dim=0;dim<nsd_;++dim)
        value += derxy_(dim,vi)*fluidvelint_old_(dim);

      value *= timefacfacrhs;

      for (int jdim=0;jdim<nsd_;++jdim)
        velforce(jdim,vi) -= value*velint_old_(jdim);
    } // vi
  }

  // 3) convective - part2
  for (int ui=0;ui<nen_;++ui)
  {
    value = timefacfac*funct_(ui);

    for (int idim=0;idim<nsd_;++idim)
    {
      const int fui = nsd_*ui+idim;

      for (int vi=0;vi<nen_;++vi)
      {
        const int fvi = nsd_*vi;

        for (int jdim=0;jdim<nsd_;++jdim)
        {
          estif_u(fvi+jdim,fui) += funct_(vi)*fluidvelxy_(jdim,idim)*value;
        }
      }
    }
  }

  // rhs at new time step
  for (int jdim=0;jdim<nsd_;++jdim)
  {
    value = 0.0; // product of fluid and adjoint velocity

    for (int dim=0;dim<nsd_;++dim)
      value+=fluidvelxy_(jdim,dim)*velint_(dim);

    value*=timefacfac;

    for (int vi=0;vi<nen_;++vi)
      velforce(jdim,vi) -= funct_(vi)*value;
  }

  // rhs at new and old time step
  if (not fluidAdjoint3Parameter_->is_stationary_)
  {
    for (int jdim=0;jdim<nsd_;++jdim)
    {
      value = 0.0; // product of fluid and adjoint velocity

      for (int dim=0;dim<nsd_;++dim)
        value+=fluidvelxy_old_(jdim,dim)*velint_old_(dim);

      value*=timefacfacrhs;

      for (int vi=0;vi<nen_;++vi)
        velforce(jdim,vi) -= funct_(vi)*value;
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ViscousGalPart(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &   estif_u,
    LINALG::Matrix<nsd_,nen_> &             velforce,
    const double &                          timefacfac,
    const double &                          timefacfacrhs
)
{
  /* viscosity term */
  /*
                   /                        \
                  |       /  \         / \   |
            2 mu  |  eps | Du | , eps | v |  |
                  |       \  /         \ /   |
                   \                        /
  */

  double value = 0.0; // helper
  const double viscdenstimefac = timefacfac*visc_/dens_;

  for (int ui=0; ui<nen_; ++ui)
  {
    const int fui   = nsd_*ui;

    for (int vi=0; vi<nen_; ++vi)
    {
      const int fvi   = nsd_*vi;

      value = 0.0;

      for (int dim=0;dim<nsd_;++dim)
        value += derxy_(dim,vi)*derxy_(dim,ui);

      value *= viscdenstimefac;

      for (int idim = 0; idim <nsd_; ++idim)
      {
        estif_u(fvi+idim,fui+idim) += value;
      } // end for (idim)
    } //vi
  } // ui

  for (int ui=0; ui<nen_; ++ui)
  {
    const int fui   = nsd_*ui; // shp fcn of u known, derivative not

    for (int jdim=0;jdim<nsd_;++jdim) // derivative of u known by jdim
    {
      value = viscdenstimefac*derxy_(jdim,ui);

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi   = nsd_*vi; // shp fcn of v known, derivative not

        for (int idim=0;idim<nsd_;++idim)
        {
          estif_u(fvi+jdim,fui+idim) += derxy_(idim,vi)*value;
        }
      } //vi
    }
  } // ui

  // viscosity at new and old time step (if instationary)
  LINALG::Matrix<nsd_,nsd_> viscstress(true);

  double viscdenstimefacrhs = timefacfacrhs*visc_/dens_; // for instationary problems
  for (int jdim = 0; jdim < nsd_; ++jdim)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      viscstress(idim,jdim)+=viscdenstimefac*(vderxy_(jdim,idim)+vderxy_(idim,jdim));

      if (not fluidAdjoint3Parameter_->is_stationary_)
        viscstress(idim,jdim)+=viscdenstimefacrhs*(vderxy_old_(jdim,idim)+vderxy_old_(idim,jdim));
    }
  }


  // computation of right-hand-side viscosity term
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      for (int jdim = 0; jdim < nsd_; ++jdim)
        velforce(idim,vi) -= viscstress(idim,jdim)*derxy_(jdim,vi);
    }
  }

  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::PressureGalPart(
    LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    const double &                            timefacfacpre,
    const double &                            timefacfacprerhs,
    const double &                            press,
    const double &                            press_old)
{
  double value = 0.0; // helper

  double pressfac = -timefacfacpre/dens_;
  for (int ui=0;ui<nen_;++ui)
  {
    value = pressfac*funct_(ui);

    for (int vi=0;vi<nen_;++vi)
    {
      const int fvi = vi*nsd_;

      for (int jdim=0;jdim<nsd_;++jdim)
        estif_p_v(fvi+jdim,ui) += derxy_(jdim,vi)*value;
    } // vi
  } // ui

  // rhs at new and old time step
  value = pressfac*press;
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int jdim=0;jdim<nsd_;++jdim) // derivative of u known by jdim
      velforce(jdim,vi) -= derxy_(jdim,vi)*value;
  } // vi

  // rhs at new and old time step
  if (not fluidAdjoint3Parameter_->is_stationary_)
  {
    value = -timefacfacprerhs/dens_*press_old;
    for (int vi=0; vi<nen_; ++vi)
    {
      for (int jdim=0;jdim<nsd_;++jdim) // derivative of u known by jdim
        velforce(jdim,vi) -= derxy_(jdim,vi)*value;
    } // vi
  }
//  for (int ui=0; ui<nen_; ++ui)
//  {
//    const double v = -timefacfacpre*funct_(ui);
//    for (int vi=0; vi<nen_; ++vi)
//    {
//      const int fvi = nsd_*vi;
      /* pressure term */
      /*
           /                \
          |                  |
          |  Dp , nabla o v  |
          |                  |
           \                /
      */
//      for (int idim = 0; idim <nsd_; ++idim)
//      {
//        estif_p_v(fvi + idim,ui) += v*derxy_(idim, vi);
//      }
//    }
//  }

  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContinuityGalPart(
    LINALG::Matrix<nen_, nen_*nsd_> &         estif_q_u,
    LINALG::Matrix<nen_,1> &                  preforce,
    const double &                            timefacfacdiv,
    const double &                            timefacfacdivrhs)
{
  double value = 0.0; // helper
  /* continuity term */
  /*
       /                \
      |                  |
      | nabla o Du  , q  |
      |                  |
       \                /
  */
  for (int vi=0;vi<nen_;++vi)
  {
    value = -timefacfacdiv*funct_(vi);

    for (int ui=0;ui<nen_;++ui)
    {
      const int fui = ui*nsd_;

      for (int idim=0;idim<nsd_;++idim)
        estif_q_u(vi,fui+idim) += value*derxy_(idim,ui);
    } // vi
  } // ui

  // rhs at new and old time step
  value = -timefacfacdiv*vdiv_;
  if (not fluidAdjoint3Parameter_->is_stationary_)
    value += -timefacfacdivrhs*vdiv_old_;

  for (int vi=0; vi<nen_; ++vi)
    preforce(vi) -= funct_(vi)*value;


  return;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::BodyForceGalPart(
    LINALG::Matrix<nsd_,nen_> &               velforce,
    const double &                            timefacfac,
    const double &                            timefacfacrhs)
{
  double value = 0.0;
  if (fluidAdjoint3Parameter_->dissipation_)
  {
    const double dissipation = fluidAdjoint3Parameter_->dissipation_fac_;

    for (int idim = 0; idim <nsd_; ++idim)
    {
      value = dissipation*reacoeff_*timefacfac*fluidvelint_(idim);

      if (fluidAdjoint3Parameter_->is_stationary_)
        value += dissipation*reacoeff_*timefacfacrhs*fluidvelint_old_(idim);

      for (int vi=0; vi<nen_; ++vi)
      {
        velforce(idim,vi)+=value*funct_(vi);
      }
    }  // end for(idim)
  }


  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::MomRes(
    LINALG::Matrix<nsd_*nsd_,nen_> &    GalMomResnU,
    LINALG::Matrix<nsd_,1> &            StrongResMomScaled,
    const double &                      timefacfac,
    const double &                      timefacfacrhs,
    const double &                      timefacfacpre,
    const double &                      timefacfacprerhs,
    const LINALG::Matrix<nsd_,nen_>&    eveln,
    const LINALG::Matrix<nsd_,nen_>&    evelnp,
    const LINALG::Matrix<nsd_,nen_>&    efluidveln,
    const LINALG::Matrix<nsd_,nen_>&    efluidvelnp
)
{
  /*
    instationary                          cross-stress, part 1
       +---+                             +-------------------+
       |   |                             |                   |

              dt*Theta   /       n+1       \        /      ~n+1       \
        Du +  --------  |   rho*u   o nabla | Du + |   rho*u   o nabla | Du +
                rho      \      (i)        /        \      (i)        /

                 /                \  n+1
              + |   rho*Du o nabla | u      +  sigma*Du
                 \                /   (i)
                |                        |     |       |
                +------------------------+     +-------+
                        Newton                  reaction
  */

  // mass matrix + reaction
  double massreacfac = 0.0; // factor summing up coefficients of reactive term and mass-matrix
  if (fluidAdjoint3Parameter_->is_stationary_)
    massreacfac = reacoeff_*timefacfac/dens_;
  else
    massreacfac = fac_+reacoeff_*timefacfac/dens_; // fac -> mass matrix // reac*timefacfac/dens -> reactive

  for (int ui=0; ui<nen_; ++ui)
  {
    const double uifunct = massreacfac*funct_(ui);

    for (int idim=0; idim<nsd_; ++idim)
      GalMomResnU(idim*nsd_+idim,ui) += uifunct;
  } // ui

  // convection
  for (int ui=0; ui<nen_; ++ui)
  {
    double value = 0.0;

    for (int dim=0;dim<nsd_;++dim)
      value -= timefacfac*fluidvelint_(dim)*derxy_(dim,ui);

    for (int idim = 0; idim <nsd_; ++idim)
      GalMomResnU(idim*nsd_+idim,ui) += value;
  } //ui

  for (int ui=0;ui<nen_;++ui)
  {
    const double uifunct = timefacfac*funct_(ui);

    for (int idim=0;idim<nsd_;++idim)
    {
      for (int jdim=0;jdim<nsd_;++jdim)
      {
        GalMomResnU(jdim+idim*nsd_,ui) += fluidvelxy_(jdim,idim)*uifunct;
      }
    }
  }

  // viscous
  LINALG::Matrix<nsd_,1> viscs(true);
  LINALG::Matrix<nsd_,1> viscs_old(true);
  if (is_higher_order_ele_)
  {
    // prework: evaluate div(eps(v))
    LINALG::Matrix<nsd_*nsd_,nen_> visc_shp(true);
    CalcDivEps(eveln,evelnp,viscs,viscs_old,visc_shp);

    // add viscous part
    GalMomResnU.Update(-2.0*timefacfac*visc_/dens_,visc_shp,1.0);
  }


  // evaluate bodyforce for strong residuum
  LINALG::Matrix<nsd_,1> bodyforce(true);
  LINALG::Matrix<nsd_,1> bodyforce_old(true);
  BodyForce(efluidveln,efluidvelnp,bodyforce,bodyforce_old);


  // residuum of momentum equation in strong form
  if (not fluidAdjoint3Parameter_->is_stationary_)
  {
    for (int idim=0;idim<nsd_;++idim)
    {
      StrongResMomScaled(idim) += velint_(idim)-velint_old_(idim) // mass term last iteration
          +(timefacfac* // velocity part of last iteration (at t^n) coming
            (dens_*(-conv1_(idim)+conv2_(idim))-2*visc_*viscs(idim)
            +reacoeff_*velint_(idim)-bodyforce(idim))
          +timefacfacpre*gradp_(idim) // pressure part of last iteration (at t^n)
          +timefacfacrhs* // last time step (= t^n+1) coming
            (dens_*(-conv1_old_(idim)+conv2_old_(idim))-2*visc_*viscs_old(idim)
            +reacoeff_*velint_old_(idim)-bodyforce_old(idim))
          +timefacfacprerhs*gradp_old_(idim))
          /dens_; // scale with density
    }
  }
  else
  {
    for (int idim=0;idim<nsd_;++idim)
    {
      StrongResMomScaled(idim) = dens_*(-conv1_(idim)+conv2_(idim))-2*visc_*viscs(idim)
                      +reacoeff_*velint_(idim)+gradp_(idim)-bodyforce(idim);
    }
  }

  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::CalcDivEps(
    const LINALG::Matrix<nsd_,nen_>&      eveln,
    const LINALG::Matrix<nsd_,nen_>&      evelnp,
    LINALG::Matrix<nsd_,1>&               viscs,
    LINALG::Matrix<nsd_,1>&               viscs_old,
    LINALG::Matrix<nsd_*nsd_,nen_>&       visc_shp
)
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

  if (nsd_==3)
  {
    for (int inode=0; inode<nen_; ++inode)
    {
      double sum = (derxy2_(0,inode)+derxy2_(1,inode)+derxy2_(2,inode));
      visc_shp(0,inode) = 0.5 * (sum + derxy2_(0,inode));
      visc_shp(1,inode) = 0.5 *  derxy2_(3,inode);
      visc_shp(2,inode) = 0.5 *  derxy2_(4,inode);
      visc_shp(3,inode) = 0.5 *  derxy2_(3,inode);
      visc_shp(4,inode) = 0.5 * (sum + derxy2_(1,inode));
      visc_shp(5,inode) = 0.5 *  derxy2_(5,inode);
      visc_shp(6,inode) = 0.5 *  derxy2_(4,inode);
      visc_shp(7,inode) = 0.5 *  derxy2_(5,inode);
      visc_shp(8,inode) = 0.5 * (sum + derxy2_(2,inode));
    }
  }
  else if (nsd_==2)
  {
    for (int inode=0; inode<nen_; ++inode)
    {
      double sum = (derxy2_(0,inode)+derxy2_(1,inode));
      visc_shp(0,inode) = 0.5 * (sum + derxy2_(0,inode));
      visc_shp(1,inode) = 0.5 * derxy2_(2,inode);
      visc_shp(2,inode) = 0.5 * derxy2_(2,inode);
      visc_shp(3,inode) = 0.5 * (sum + derxy2_(1,inode));
    }
  }
  else dserror("Epsilon(N) is not implemented for the 1D case");

  for (int inode=0; inode<nen_; ++inode)
  {
    for (int idim=0; idim<nsd_; ++idim)
    {
      const int nsd_idim = idim*nsd_;

      for (int jdim=0; jdim<nsd_; ++jdim)
        viscs(idim) += visc_shp(nsd_idim+jdim,inode)*eveln(jdim,inode);
    }
  }

  if (not fluidAdjoint3Parameter_->is_stationary_)
  {
    for (int inode=0; inode<nen_; ++inode)
    {
      for (int idim=0; idim<nsd_; ++idim)
      {
        const int nsd_idim = idim*nsd_;

        for (int jdim=0; jdim<nsd_; ++jdim)
          viscs_old(idim) += visc_shp(nsd_idim+jdim,inode)*evelnp(jdim,inode);
      }
    }
  }

  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::BodyForce(
           const LINALG::Matrix<nsd_,nen_>&           efluidveln,
           const LINALG::Matrix<nsd_,nen_>&           efluidvelnp,
           LINALG::Matrix<nsd_,1>&                    bodyforce,
           LINALG::Matrix<nsd_,1>&                    bodyforce_old
)
{
  bodyforce.Clear();
  bodyforce_old.Clear();

  if (fluidAdjoint3Parameter_->dissipation_)
  {
    const double dissipation = fluidAdjoint3Parameter_->dissipation_fac_;

    /* ------------------------------------------------------------------------ *
     * 1) evaluate bodyforce at new time step                                   *
     * ------------------------------------------------------------------------ */

    // dissipation term due to reaction
    bodyforce.Update(2*dissipation*reacoeff_,fluidvelint_);

    // dissipation term due to viscosity
    if (is_higher_order_ele_)
    {
      LINALG::Matrix<nsd_,numderiv2_> fluidvelxy2(true);
      fluidvelxy2.MultiplyNT(efluidveln,derxy2_);

      LINALG::Matrix<nsd_,1> laplaceU(true);
      for (int idim=0;idim<nsd_;++idim)
      {
        for (int jdim=0;jdim<nsd_;++jdim)
          laplaceU(idim) += fluidvelxy2(idim,jdim);
      }

      bodyforce.Update(-dissipation*visc_,laplaceU,1.0);
    }


    /* ------------------------------------------------------------------------ *
     * 2) evaluate bodyforce at old time step in instationary case              *
     * ------------------------------------------------------------------------ */
    if (not fluidAdjoint3Parameter_->is_stationary_)
    {
      bodyforce_old.Update(2*dissipation*reacoeff_,fluidvelint_old_);

      // dissipation term due to viscosity
      if (is_higher_order_ele_)
      {
        LINALG::Matrix<nsd_,numderiv2_> fluidvelxy2_old(true);
        fluidvelxy2_old.MultiplyNT(efluidvelnp,derxy2_);

        LINALG::Matrix<nsd_,1> laplaceU_old(true);
        for (int idim=0;idim<nsd_;++idim)
        {
          for (int jdim=0;jdim<nsd_;++jdim)
            laplaceU_old(idim) += fluidvelxy2_old(idim,jdim);
        }

        bodyforce_old.Update(-dissipation*visc_,laplaceU_old,1.0);
      }
    }
  }
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::PSPG(
    LINALG::Matrix<nen_, nen_*nsd_> &         estif_q_u,
    LINALG::Matrix<nen_,nen_> &               ppmat,
    LINALG::Matrix<nen_,1> &                  preforce,
    LINALG::Matrix<nsd_*nsd_,nen_> &          GalMomResnU,
    LINALG::Matrix<nsd_,1> &                  StrongResMomScaled,
    const double &                            timefacfac,
    const double &                            timefacfacrhs,
    const double &                            timefacfacpre,
    const double &                            timefacfacprerhs)
{
  const double tau=tau_(1);

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
          estif_q_u(vi,fui_p_jdim) += tau*derxy_(idim,vi)*GalMomResnU(nsd_idim+jdim,ui);
        } // jdim
      } // vi
    } // ui
  } //idim

  for (int ui=0; ui<nen_; ++ui)
  {
    for (int idim = 0; idim <nsd_; ++idim)
    {
      const double v=timefacfacpre*derxy_(idim,ui)*tau;

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
        ppmat(vi,ui)+=v*derxy_(idim,vi);
      } // vi
    } // end for(idim)
  }  // ui

  // rhs for new and old time step
  for (int idim = 0; idim <nsd_; ++idim)
  {
    const double resmom_scaled = -tau*StrongResMomScaled(idim);

    for (int vi=0; vi<nen_; ++vi)
    {
      // pressure stabilisation
      preforce(vi) += derxy_(idim, vi)*resmom_scaled;
    }
  } // end for(idim)
  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::SUPG(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
    LINALG::Matrix<nen_*nsd_,nen_> &          estif_p_v,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    LINALG::Matrix<nsd_*nsd_,nen_> &          GalMomResnU,
    LINALG::Matrix<nsd_,1> &                  StrongResMomScaled,
    const double &                            timefacfac,
    const double &                            timefacfacrhs,
    const double &                            timefacfacpre,
    const double &                            timefacfacprerhs)
{
  /*
                    /                                \
                   |  ~n+af    /     n+af       \     |
                 - |  u     , | rho*u    o nabla | v  |
                   |           \     (i)        /     |
                    \                                /
   */

  double supgfac=-dens_*tau_(0);

  LINALG::Matrix<nen_,1> supg_test(true);
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int dim=0;dim<nsd_;++dim)
    {
      supg_test(vi)+=supgfac*derxy_(dim,vi)*fluidvelint_(dim);
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

          estif_u(fvi_p_idim,fui_p_jdim) += supg_test(vi)*GalMomResnU(nsd_idim_p_jdim,ui);
        } // jdim
      } // vi
    } // ui
  } //idim

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
    const double v = timefacfacpre*supg_test(vi);

    for (int idim = 0; idim <nsd_; ++idim)
    {
      const int fvi   = nsd_*vi + idim;

      for (int ui=0; ui<nen_; ++ui)
      {
        estif_p_v(fvi,ui) += v*derxy_(idim, ui);
      }
    }
  }  // end for(idim)


  // rhs for new and old time step
  for (int idim = 0; idim <nsd_; ++idim)
  {
    for (int vi=0; vi<nen_; ++vi)
      velforce(idim,vi) += supg_test(vi)*StrongResMomScaled(idim);
  }  // end for(idim)
  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContRes(
    double &                                  StrongResContScaled,
    const double &                            timefacfacdiv,
    const double &                            timefacfacdivrhs
)
{
  double contforce = 0.0;
  double contforce_old = 0.0;

  ContForce(contforce,contforce_old);

  StrongResContScaled = timefacfacdiv*(vdiv_-contforce);

  if (not fluidAdjoint3Parameter_->is_stationary_)
    StrongResContScaled += timefacfacdivrhs*(vdiv_old_-contforce_old);
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContForce(
    double &                            contForce,
    double &                            contForce_old
)
{
  contForce = 0.0;
  contForce_old = 0.0;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ContStab(
    LINALG::Matrix<nen_*nsd_,nen_*nsd_> &     estif_u,
    LINALG::Matrix<nsd_,nen_> &               velforce,
    const double &                            timefacfacdiv,
    const double &                            timefacfacdivrhs
)
{
  double cstabfac = timefacfacdiv*tau_(2);
  double value = 0.0;

  /* continuity stabilisation on left hand side */
  /*
              /                        \
             |                          |
        tauC | nabla o Du  , nabla o v  |
             |                          |
              \                        /
  */

  for (int ui=0; ui<nen_; ++ui)
  {
    const int fui = nsd_*ui;

    for (int idim = 0; idim <nsd_; ++idim)
    {
      const int fui_p_idim = fui+idim;

      value = cstabfac*derxy_(idim,ui);

      for (int vi=0; vi<nen_; ++vi)
      {
        const int fvi = nsd_*vi;

        for(int jdim=0;jdim<nsd_;++jdim)
        {
          estif_u(fvi+jdim,fui_p_idim) += value*derxy_(jdim, vi) ;
        }
      }
    } // end for(idim)
  }

  double StrongResContScaled = 0.0;

  ContRes(StrongResContScaled,
      timefacfacdiv,
      timefacfacdivrhs);

  // computation of rhs viscosity term at new time step
  for (int vi=0; vi<nen_; ++vi)
  {
    for (int idim = 0; idim < nsd_; ++idim)
    {
      /* viscosity term on right-hand side */
      velforce(idim,vi)+= derxy_(idim,vi)*StrongResContScaled;
    }
  }

  return;
}



template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3Impl<distype>::ExtractValuesFromGlobalVector(
    const DRT::Discretization&   discretization, ///< discretization
    const vector<int>&           lm,             ///<
    LINALG::Matrix<nsd_,nen_> *  matrixtofill,   ///< vector field
    LINALG::Matrix<nen_,1> *     vectortofill,   ///< scalar field
    const std::string            state          ///< state of the global vector
)
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(state);
  if(matrix_state == null)
    dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  for (int inode=0; inode<nen_; ++inode)  // number of nodes
  {
    // fill a vector field via a pointer
    if (matrixtofill != NULL)
    {
      for(int idim=0; idim<nsd_; ++idim) // number of dimensions
      {
        (*matrixtofill)(idim,inode) = mymatrix[idim+(inode*numdofpernode_)];
      }  // end for(idim)
    }
    // fill a scalar field via a pointer
    if (vectortofill != NULL)
      (*vectortofill)(inode,0) = mymatrix[nsd_+(inode*numdofpernode_)];
  }
}



#endif
