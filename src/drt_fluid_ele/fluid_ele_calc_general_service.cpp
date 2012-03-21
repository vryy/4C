/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_general_service.cpp

\brief general service routines for calculation of fluid element

<pre>
Maintainer: Volker Gravemeier & Andreas Ehrl
            {vgravem,ehrl}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15245/15252
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_calc.H"
#include "fluid_ele.H"
#include "fluid_ele_parameter.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"

#include "../drt_geometry/position_array.H"

#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/newtonianfluid.H"

#include "../drt_nurbs_discret/drt_nurbs_utils.H"

/*----------------------------------------------------------------------*
 * Action type: Integrate shape function
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::IntegrateShapeFunction(
    DRT::ELEMENTS::Fluid3*    ele,
    DRT::Discretization&      discretization,
    vector<int>&              lm            ,
    Epetra_SerialDenseVector& elevec1       )
{
  // --------------------------------------------------
  // construct views
  LINALG::Matrix<numdofpernode_*nen_,    1> vector(elevec1.A(),true);

  // get Gaussrule
  //const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  //----------------------------------------------------------------
  // Now do the nurbs specific stuff (for isogeometric elements)
  //----------------------------------------------------------------
  if(isNurbs_)
  {
    // access knots and weights for this element
    bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization,ele,myknots_,weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if(zero_size)
      return(0);
  } // Nurbs specific stuff

  if (ele->IsAle())
  {
    LINALG::Matrix<nsd_,nen_>       edispnp(true);
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"dispnp");

    // get new node positions for isale
     xyze_ += edispnp;
  }

//------------------------------------------------------------------
//                       INTEGRATION LOOP
//------------------------------------------------------------------

  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad,ele->Id());

    for (int ui=0; ui<nen_; ++ui) // loop rows  (test functions)
    {
      // integrated shape function is written into the pressure dof
      int fuippp=numdofpernode_*ui+nsd_;
      vector(fuippp)+=fac_*funct_(ui);
    }
  }

  return 0;
}


/*----------------------------------------------------------------------*
 * Action type: Compute Error                              shahmiri 01/12
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::ComputeError(
    DRT::ELEMENTS::Fluid3*          ele,
    ParameterList&                  params,
    Teuchos::RCP<MAT::Material>&    mat,
    DRT::Discretization&            discretization,
    vector<int>&                    lm,
    Epetra_SerialDenseVector&       elevec1
    )
{
  // integrations points and weights
  // more GP than usual due to (possible) cos/exp fcts in analytical solutions
  // degree 5
  const DRT::UTILS::GaussIntegration intpoints(distype, 5);
  return ComputeError( ele, params, mat,
                       discretization, lm,
                       elevec1, intpoints);
}
/*----------------------------------------------------------------------*
 * Action type: Compute Error                              shahmiri 01/12
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::ComputeError(
    DRT::ELEMENTS::Fluid3*          ele,
    ParameterList&                  params,
    Teuchos::RCP<MAT::Material>&    mat,
    DRT::Discretization&            discretization,
    vector<int>&                    lm,
    Epetra_SerialDenseVector&       elevec1,
    const DRT::UTILS::GaussIntegration & intpoints
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
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelaf, &epreaf,"u and p at time n+1 (converged)");

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  //----------------------------------------------------------------
  // Now do the nurbs specific stuff (for isogeometric elements)
  //----------------------------------------------------------------
  if(isNurbs_)
  {
    // access knots and weights for this element
    bool zero_size = DRT::NURBS::GetMyNurbsKnotsAndWeights(discretization,ele,myknots_,weights_);

    // if we have a zero sized element due to a interpolated point -> exit here
    if(zero_size)
      return(0);
  } // Nurbs specific stuff

  if (ele->IsAle())
  {
    LINALG::Matrix<nsd_,nen_>       edispnp(true);
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"dispnp");

    // get new node positions for isale
     xyze_ += edispnp;
  }

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

         const double t = f3Parameter_->time_;

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
      }
      break;
    case INPAR::FLUID::shear_flow:
      {
        const double maxvel = 1.0;
        const double hight = 1.0;

        // y=0 is located in the middle of the domain
        if (nsd_ == 2)
        {
          p = 1.0;
          u(0) = xyzint(1)*maxvel + hight/2*maxvel;
          u(1) = 0.0;
        }
        if (nsd_ == 3)
        {
          p = 0.0;
          u(0) = xyzint(1)*maxvel + hight/2*maxvel;
          u(1) = 0.0;
          u(2) = 0.0;
        }
      }
      break;
    case INPAR::FLUID::gravitation:
      {
        const double gravity = 10.0;
        const double hight = 1.0;

        // 2D: rectangle 1.0x1.0
        // 3D: cube 1.0x1.0x1.0
        // y=0 is located in the middle of the domain
        if (nsd_ == 2)
        {
          p = -xyzint(1)*gravity + hight/2*gravity;
          u(0) = 0.0;
          u(1) = 0.0;
        }
        if (nsd_ == 3)
        {
          p = -xyzint(1)*gravity + hight/2*gravity;
          u(0) = 0.0;
          u(1) = 0.0;
          u(2) = 0.0;
        }
      }
      break;
    case INPAR::FLUID::channel2D:
      {
        const double maxvel=1.25;
        const double hight = 1.0;
        const double visc = 1.0;
        const double pressure_gradient = 10.0;

        // u_max = 1.25
        // y=0 is located in the middle of the channel
        if (nsd_ == 2)
        {
          p = 1.0;
          //p = -10*xyzint(0)+20;
          u(0) = maxvel -((hight*hight)/(2.0*visc)*pressure_gradient*(xyzint(1)/hight)*(xyzint(1)/hight));
          u(1) = 0.0;
        }
        else
          dserror("3D analytical solution is not implemented yet");
      }
      break;
    case INPAR::FLUID::jeffery_hamel_flow:
    {
      //LINALG::Matrix<3,1> physpos(true);
      //GEO::elementToCurrentCoordinates(distype, xyzint, xsi_, physpos);

      double position[2];
      position[0] = xyzint(0);
      position[1] = xyzint(1);

      if (1.0 < position[0] and position[0] < 2.0 and 0.0 < position[1] and position[1] < position[0])
      {
        const double u_exact_x = DRT::Problem::Instance()->Funct(0).Evaluate(0,position,0.0,NULL);
        const double u_exact_y = DRT::Problem::Instance()->Funct(0).Evaluate(1,position,0.0,NULL);
        u(0) = u_exact_x;
        u(1) = u_exact_y;
      }

    }
      break;
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


/*!
 * \brief fill elment matrix and vectors with the global values
 */
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalc<distype>::ExtractValuesFromGlobalVector( const DRT::Discretization&   discretization, ///< discretization
                                    const vector<int>&           lm,             ///<
                                    FLD::RotationallySymmetricPeriodicBC<distype> & rotsymmpbc, ///<
                                    LINALG::Matrix<nsd_,nen_> *  matrixtofill,   ///< vector field
                                    LINALG::Matrix<nen_,1> *     vectortofill,   ///< scalar field
                                    const std::string            state)          ///< state of the global vector
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(state);
  if(matrix_state == null)
    dserror("Cannot get state vector %s", state.c_str());

  // extract local values of the global vectors
  std::vector<double> mymatrix(lm.size());
  DRT::UTILS::ExtractMyValues(*matrix_state,mymatrix,lm);

  // rotate the vector field in the case of rotationally symmetric boundary conditions
  if(matrixtofill != NULL)
    rotsymmpbc.RotateMyValuesIfNecessary(mymatrix);

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


/*--------------------------------------------------------------------------------
 * additional output for turbulent channel flow                    rasthofer 12/10
 * -> dissipation
 *--------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::CalcDissipation(
  Fluid3*                    ele,
  ParameterList&             params,
  DRT::Discretization&       discretization,
  vector<int>&               lm,
  RefCountPtr<MAT::Material> mat)
{
  // TODO: Ursula
  dserror("This function is currently not used! -> Check it!");
  // -> has to be adapted to changes in Sysmat()
  // -> implemented for scale similarity model
  // -> not adapted to multifractal subgrid-scale modeling
  // -> (residual-based) cross- and Reynolds-stress terms not included, yet
  // -> has to be merged with corresponding function of np-gen-alpha (Gammis Code)

//  // create matrix objects for nodal values
//  // and extract velocities, pressure and accelerations
//  // from the global distributed vectors
//  LINALG::Matrix<my::nen_,1> epre;
//  LINALG::Matrix<my::nsd_,my::nen_> evel;
//  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &evel, &epre,"vel");
//  LINALG::Matrix<my::nsd_,my::nen_> eacc;
//  my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &eacc, NULL,"acc");
//  LINALG::Matrix<my::nsd_,my::nen_> fsevel(true);
//  if (my::f3Parameter_->fssgv_ != INPAR::FLUID::no_fssgv)
//  {
//    my::ExtractValuesFromGlobalVector(discretization,lm, *my::rotsymmpbc_, &fsevel, NULL,"fsvel");
//  }
//  LINALG::Matrix<my::nsd_,my::nen_> evel_hat;
//  LINALG::Matrix<my::nsd_*my::nsd_,my::nen_> ereynoldsstress_hat;
//  if (my::f3Parameter_->turb_mod_action_ == INPAR::FLUID::scale_similarity_basic)
//  {
//    RCP<Epetra_MultiVector> filtered_vel = params.get<RCP<Epetra_MultiVector> >("filtered vel");
//    RCP<Epetra_MultiVector> filtered_reystre = params.get<RCP<Epetra_MultiVector> >("filtered reystr");
//    for (int nn=0;nn<my::nen_;++nn)
//    {
//      int lid = (ele->Nodes()[nn])->LID();
//
//      for (int dimi=0;dimi<3;++dimi)
//      {
//        evel_hat(dimi,nn) = (*((*filtered_vel)(dimi)))[lid];
//        for (int dimj=0;dimj<3;++dimj)
//        {
//          int index=3*dimi+dimj;
//          ereynoldsstress_hat(index,nn) = (*((*filtered_reystre)(index)))[lid];
//        }
//      }
//    }
//  }
//
//
//  // the coordinates of the element layers in the channel
//  // planecoords are named nodeplanes in turbulence_statistics_channel!
//  RefCountPtr<vector<double> > planecoords  = params.get<RefCountPtr<vector<double> > >("planecoords_",Teuchos::null);
//  if(planecoords==Teuchos::null)
//    dserror("planecoords is null, but need channel_flow_of_height_2\n");
//
//  //this will be the y-coordinate of a point in the element interior
//  double center = 0.0;
//  // get node coordinates of element
//  for(int inode=0;inode<my::nen_;inode++)
//  {
//    my::xyze_(0,inode)=ele->Nodes()[inode]->X()[0];
//    my::xyze_(1,inode)=ele->Nodes()[inode]->X()[1];
//    my::xyze_(2,inode)=ele->Nodes()[inode]->X()[2];
//
//    center+=my::xyze_(1,inode);
//  }
//  center/=my::nen_;
//
//
//  // ---------------------------------------------------------------------
//  // calculate volume
//  // ---------------------------------------------------------------------
//  // evaluate shape functions and derivatives at element center
//  my::EvalShapeFuncAndDerivsAtEleCenter(ele->Id());
//  // element area or volume
//  const double vol = fac_;
//
//  // get velocity at integration point
//  my::velint_.Multiply(evel,my::funct_);
//  // convective term
//  convvelint_.Update(velint_);
//
//  if (f3Parameter_->mat_gp_ or f3Parameter_->tau_gp_)
//   dserror ("Evaluation of material or stabilization parameters at gauss point not supported,yet!");
//  // ---------------------------------------------------------------------
//  // get material
//  // ---------------------------------------------------------------------
//  if (mat->MaterialType() == INPAR::MAT::m_fluid)
//  {
//    const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());
//
//    // get constant viscosity
//    visc_ = actmat->Viscosity();
//    // get constant density
//    densaf_ = actmat->Density();
//    densam_ = densaf_;
//    densn_  = densaf_;
//  }
//  else dserror("Only material m_fluid supported");
//  densaf_ = 1.0;
//  if (f3Parameter_->physicaltype_ != INPAR::FLUID::incompressible)
//    dserror("CalcDissipation() only for incompressible flows!");
//
//
//  // ---------------------------------------------------------------------
//  // calculate turbulent viscosity at element center
//  // ---------------------------------------------------------------------
//  double Cs_delta_sq = ele->CsDeltaSq();
//  double visceff = visc_;
//  if (f3Parameter_->turb_mod_action_ == INPAR::FLUID::smagorinsky or f3Parameter_->turb_mod_action_ == INPAR::FLUID::dynamic_smagorinsky)
//  {
//    CalcSubgrVisc(evel,vol,f3Parameter_->Cs_,Cs_delta_sq,f3Parameter_->l_tau_);
//    // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
//    visceff += sgvisc_;
//  }
//  else if (f3Parameter_->fssgv_ != INPAR::FLUID::no_fssgv)
//    CalcFineScaleSubgrVisc(evel,fsevel,vol,f3Parameter_->Cs_);
//
//  // ---------------------------------------------------------------------
//  // calculate stabilization parameters at element center
//  // ---------------------------------------------------------------------
//  // Stabilization parameter
//  CalcStabParameter(vol);
//  const double tau_M       = tau_(0);
//  const double tau_Mp      = tau_(1);
//  const double tau_C       = tau_(2);
//
//  // ---------------------------------------------------------------------
//  // get bodyforce
//  // ---------------------------------------------------------------------
//  LINALG::Matrix<nsd_,nen_> ebofo(true);
//  LINALG::Matrix<nsd_,nen_> epgrad(true);
//  LINALG::Matrix<nen_,1>    escabofo(true);
//  BodyForce(ele, f3Parameter_, ebofo, epgrad, escabofo);
//
//  // working arrays for the quantities we want to compute
//  double eps_visc        = 0.0;
//  double eps_smag        = 0.0;
//  double eps_avm3        = 0.0;
//  double eps_scsim       = 0.0;
//  double eps_scsimfs     = 0.0;
//  double eps_scsimbs     = 0.0;
//  double eps_supg        = 0.0;
//  double eps_cstab       = 0.0;
//  double eps_pspg        = 0.0;
//
//  // gaussian points
//  //const DRT::UTILS::IntPointsAndWeights<nsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);
//  //DRT::UTILS::GaussIntegration intpoints( distype );
//  //------------------------------------------------------------------
//  //                       INTEGRATION LOOP
//  //------------------------------------------------------------------
//  //for (int iquad=0;iquad<intpoints.IP().nquad;++iquad)
//  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
//  {
//    // evaluate shape functions and derivatives at integration point
//    //EvalShapeFuncAndDerivsAtIntPoint(intpoints,iquad,ele->Id());
//    EvalShapeFuncAndDerivsAtIntPoint(iquad,ele->Id());
//
//
//    // get velocity at integration point
//    velint_.Multiply(evel,funct_);
//    // get velocity derivatives at integration point
//    vderxy_.MultiplyNT(evel,derxy_);
//
//    if (f3Parameter_->fssgv_ != INPAR::FLUID::no_fssgv)
//      fsvderxy_.MultiplyNT(fsevel,derxy_);
//    // get pressure gradient at integration point
//    gradp_.Multiply(derxy_,epre);
//    // convective term
//    convvelint_.Update(velint_);
//    conv_old_.Multiply(vderxy_,convvelint_);
//    // get bodyforce at integration point
//    bodyforce_.Multiply(ebofo,funct_);
//    // prescribed pressure gradient
//    prescribedpgrad_.Multiply(epgrad,funct_);
//    // get acceleration at integration point
//    accint_.Multiply(eacc,funct_);
//
//    if(f3Parameter_->turb_mod_action_ == INPAR::FLUID::scale_similarity_basic)
//    {
//      reystressinthat_.Clear();
//      velinthat_.Clear();
//      // get filtered velocity at integration point
//      velinthat_.Multiply(evel_hat,funct_);
//      // get filtered reynoldsstress at integration point
//      for (int dimi=0;dimi<nsd_;dimi++)
//      {
//        for (int dimj=0;dimj<nsd_;dimj++)
//        {
//          for (int inode=0;inode<nen_;inode++)
//          {
//            reystressinthat_(dimi,dimj) += funct_(inode) * ereynoldsstress_hat(3*dimi+dimj,inode);
//          }
//        }
//      }
//    }
//    // calculate residual of momentum equation at integration point
    /*
                               /                  \
      r    (x) = acc    (x) + | vel    (x) o nabla | vel    (x) +
       M                       \                  /

                 + nabla p    - f             (not higher order, i.e. 2 * visceff * nabla o eps(vel))
    */
//    for (int rr=0;rr<nsd_;rr++)
//    {
//      momres_old_(rr,0) = densaf_ * (accint_(rr,0) + conv_old_(rr,0) + gradp_(rr,0) - bodyforce_(rr,0)) - prescribedpgrad_(rr,0);
//    }
//    // get second derivative of the viscous term:
//    // div(epsilon(u))
//    if (is_higher_order_ele_)
//    {
//      CalcDivEps(evel);
//      for (int rr=0;rr<nsd_;rr++)
//      {
//        momres_old_(rr,0) -= 2*visceff*visc_old_(rr,0);
//      }
//    }
//    else
//    {
//      viscs2_.Clear();
//      visc_old_.Clear();
//    }
//
//    // calculate residual of continuity equation integration point
//    vdiv_ = 0.0;
//    for (int rr=0;rr<nsd_;rr++)
//    {
//      vdiv_ += vderxy_(rr, rr);
//    }
//
//    LINALG::Matrix<nsd_,nsd_> two_epsilon;
//    for(int rr=0;rr<nsd_;++rr)
//    {
//      for(int mm=0;mm<nsd_;++mm)
//      {
//        two_epsilon(rr,mm) = vderxy_(rr,mm) + vderxy_(mm,rr);
//      }
//    }
//
//    // contribution of this gausspoint to viscous energy
//    // dissipation (Galerkin)
    /*
                     /                                \
                    |       / n+1 \         / n+1 \   |
          2* visc * |  eps | u     | , eps | u     |  |
                    |       \     /         \     /   |
                     \                                /
    */
//    for(int rr=0;rr<nsd_;++rr)
//    {
//      for(int mm=0;mm<nsd_;++mm)
//      {
//        eps_visc += 0.5*visc_*fac_*two_epsilon(rr,mm)*two_epsilon(mm,rr);
//      }
//    }
//
//    // contribution of this gausspoint to viscous energy
//    // dissipation (Smagorinsky)
    /*
                         /                                \
                        |       / n+1 \         / n+1 \   |
          2* visc    *  |  eps | u     | , eps | u     |  |
                 turb   |       \     /         \     /   |
                         \                                /
    */
//    if(f3Parameter_->turb_mod_action_ == INPAR::FLUID::dynamic_smagorinsky or f3Parameter_->turb_mod_action_ == INPAR::FLUID::smagorinsky)
//    {
//      for(int rr=0;rr<nsd_;++rr)
//      {
//        for(int mm=0;mm<nsd_;++mm)
//        {
//          eps_smag += 0.5*sgvisc_*fac_*two_epsilon(rr,mm)*two_epsilon(mm,rr);
//        }
//      }
//    }
//
//    // contribution of this gausspoint to viscous energy
//    // dissipation (AVM3)
    /*
                         /                                \
                        |       /  n+1 \         / n+1 \   |
          2* visc    *  |  eps | du     | , eps | u     |  |
                 turb   |       \      /         \     /   |
                         \                                /
    */
//    if (f3Parameter_->fssgv_ != INPAR::FLUID::no_fssgv)
//    {
//      LINALG::Matrix<nsd_,nsd_> fstwo_epsilon;
//      for(int rr=0;rr<nsd_;++rr)
//      {
//        for(int mm=0;mm<nsd_;++mm)
//        {
//          fstwo_epsilon(rr,mm) = fsvderxy_(rr,mm) + fsvderxy_(mm,rr);
//        }
//      }
//      for(int rr=0;rr<nsd_;++rr)
//      {
//        for(int mm=0;mm<nsd_;++mm)
//        {
////          eps_avm3 += 0.5*fssgvisc_*fac_*fstwo_epsilon(rr,mm)*two_epsilon(mm,rr);
//          eps_avm3 += 0.5*fssgvisc_*fac_*fstwo_epsilon(rr,mm)*fstwo_epsilon(mm,rr);
//        }
//      }
//    }
//
//    // contribution of this gausspoint to viscous energy
//    // dissipation (Scale Similarity)
    /*
             /                                \
            |   ssm  /^n+1 \         / n+1 \   |
            |  tau  | u     | , eps | u     |  |
            |        \     /         \     /   |
             \                                /
    */
//    if(f3Parameter_->turb_mod_action_ == INPAR::FLUID::scale_similarity_basic)
//    {
//      LINALG::Matrix<nsd_,nsd_> tau_scale_sim;
//      for(int rr=0;rr<nsd_;++rr)
//      {
//        for(int mm=0;mm<nsd_;++mm)
//        {
//          tau_scale_sim(rr,mm) = reystressinthat_(rr,mm) - velinthat_(rr) * velinthat_(mm);
//        }
//      }
//
//      //old version
//        double Production = 0.0;
//
//        for (int dimi=0;dimi<nsd_;dimi++)
//        {
//          for (int dimj=0;dimj<nsd_;dimj++)
//          {
//            Production += - tau_scale_sim(dimi,dimj)*0.5*two_epsilon(dimi,dimj);
//          }
//        }
//
//      // dissipation due to scale similarity model
//      for(int rr=0;rr<nsd_;++rr)
//      {
//        for(int mm=0;mm<nsd_;++mm)
//        {
//          eps_scsim += -0.5*fac_*densaf_*f3Parameter_->Cl_*tau_scale_sim(rr,mm)*two_epsilon(mm,rr);
//        }
//      }
//      if (Production >= 0.0)
//      {
//        // forwardscatter
//        for(int rr=0;rr<nsd_;++rr)
//        {
//          for(int mm=0;mm<nsd_;++mm)
//          {
//            eps_scsimfs += -0.5*fac_*densaf_*f3Parameter_->Cl_*tau_scale_sim(rr,mm)*two_epsilon(mm,rr);
//          }
//        }
//      }
//      else
//      {
//        // backscatter
//        for(int rr=0;rr<nsd_;++rr)
//        {
//          for(int mm=0;mm<nsd_;++mm)
//          {
//            eps_scsimbs += -0.5*fac_*densaf_*f3Parameter_->Cl_*tau_scale_sim(rr,mm)*two_epsilon(mm,rr);
//          }
//        }
//      }
//    }
//
//    // contribution of this gausspoint to energy
//    // dissipation by supg-stabilization
//    if (f3Parameter_->supg_ == INPAR::FLUID::convective_stab_supg)
//    {
//      for (int rr=0;rr<nsd_;rr++)
//      {
//        eps_supg += densaf_ * fac_ * conv_old_(rr,0) * tau_M * momres_old_(rr,0);
//      }
//    }
//
//    // contribution of this gausspoint to energy
//    // dissipation by continuity-stabilization
//    if (f3Parameter_->cstab_ == INPAR::FLUID::continuity_stab_yes)
//    {
//      eps_cstab += fac_ * vdiv_ * tau_C * vdiv_;
//    }
//
//    // contribution of this gausspoint to energy
//    // dissipation by pspg-stabilization
//    if (f3Parameter_->pspg_ == INPAR::FLUID::pstab_use_pspg)
//    {
//      for (int rr=0;rr<nsd_;rr++)
//      {
//        eps_pspg += fac_ * gradp_(rr,0) * tau_Mp * momres_old_(rr,0);
//      }
//    }
//
//    velint_.Clear();
//    vderxy_.Clear();
//    fsvderxy_.Clear();
//    gradp_.Clear();
//    convvelint_.Clear();
//    conv_old_.Clear();
//    bodyforce_.Clear();
//    prescribedpgrad_.Clear();
//    accint_.Clear();
//    velinthat_.Clear();
//    reystressinthat_.Clear();
//    momres_old_.Clear();
//    viscs2_.Clear();
//    visc_old_.Clear();
//    vdiv_ = 0.0;
//
//  }
//
//
//  eps_visc /= vol;
//  eps_smag /= vol;
//  eps_avm3 /= vol;
//  eps_scsim /= vol;
//  eps_scsimfs /= vol;
//  eps_scsimbs /= vol;
//  eps_supg /= vol;
//  eps_cstab /= vol;
//  eps_pspg /= vol;
//
//  RefCountPtr<vector<double> > incrvol           = params.get<RefCountPtr<vector<double> > >("incrvol"          );
//
//  RefCountPtr<vector<double> > incr_eps_visc      = params.get<RefCountPtr<vector<double> > >("incr_eps_visc"    );
//  RefCountPtr<vector<double> > incr_eps_eddyvisc  = params.get<RefCountPtr<vector<double> > >("incr_eps_eddyvisc");
//  RefCountPtr<vector<double> > incr_eps_avm3      = params.get<RefCountPtr<vector<double> > >("incr_eps_avm3"    );
//  RefCountPtr<vector<double> > incr_eps_scsim     = params.get<RefCountPtr<vector<double> > >("incr_eps_scsim"   );
//  RefCountPtr<vector<double> > incr_eps_scsimfs   = params.get<RefCountPtr<vector<double> > >("incr_eps_scsimfs" );
//  RefCountPtr<vector<double> > incr_eps_scsimbs   = params.get<RefCountPtr<vector<double> > >("incr_eps_scsimbs" );
//  RefCountPtr<vector<double> > incr_eps_supg      = params.get<RefCountPtr<vector<double> > >("incr_eps_supg"    );
//  RefCountPtr<vector<double> > incr_eps_cstab     = params.get<RefCountPtr<vector<double> > >("incr_eps_cstab"   );
//  RefCountPtr<vector<double> > incr_eps_pspg      = params.get<RefCountPtr<vector<double> > >("incr_eps_pspg"    );
//
//  bool found = false;
//
//  int nlayer = 0;
//  for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
//  {
//    if(center<(*planecoords)[nlayer+1])
//    {
//      found = true;
//      break;
//    }
//    nlayer++;
//  }
//  if (found ==false)
//  {
//    dserror("could not determine element layer");
//  }
//
//  // collect layer volume
//  (*incrvol      )[nlayer] += vol;
//
//  (*incr_eps_visc    )[nlayer] += eps_visc       ;
//  (*incr_eps_eddyvisc)[nlayer] += eps_smag       ;
//  (*incr_eps_avm3    )[nlayer] += eps_avm3       ;
//  (*incr_eps_scsim   )[nlayer] += eps_scsim      ;
//  (*incr_eps_scsimfs )[nlayer] += eps_scsimfs    ;
//  (*incr_eps_scsimbs )[nlayer] += eps_scsimbs    ;
//  (*incr_eps_supg    )[nlayer] += eps_supg       ;
//  (*incr_eps_cstab   )[nlayer] += eps_cstab      ;
//  (*incr_eps_pspg    )[nlayer] += eps_pspg       ;

  return 0;
}


/*!
      \brief do finite difference check for given element ID
             --> for debugging purposes only

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
void DRT::ELEMENTS::FluidEleCalc<distype>::FDcheck(
  int                                                 eid,
  const LINALG::Matrix<nsd_,nen_>&                    evelaf,
  const LINALG::Matrix<nsd_,nen_>&                    eveln,
  const LINALG::Matrix<nsd_,nen_>&                    fsevelaf,
  const LINALG::Matrix<nen_,1>&                       epreaf,
  const LINALG::Matrix<nsd_,nen_>&                    eaccam,
  const LINALG::Matrix<nen_,1>&                       escaaf,
  const LINALG::Matrix<nen_,1>&                       escaam,
  const LINALG::Matrix<nen_,1>&                       escadtam,
  const LINALG::Matrix<nsd_,nen_>&                    emhist,
  const LINALG::Matrix<nsd_,nen_>&                    edispnp,
  const LINALG::Matrix<nsd_,nen_>&                    egridv,
  const LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  estif,
  const LINALG::Matrix<(nsd_+1)*nen_,(nsd_+1)*nen_>&  emesh,
  const LINALG::Matrix<(nsd_+1)*nen_,    1>&          eforce,
  const double                                        thermpressaf,
  const double                                        thermpressam,
  const double                                        thermpressdtaf,
  const double                                        thermpressdtam,
  const Teuchos::RCP<const MAT::Material>             material,
  const double                                        timefac,
  const double&                                       Cs,
  const double&                                       Cs_delta_sq,
  const double&                                       l_tau)
{
  // magnitude of dof perturbation
  const double epsilon=1e-14;

  // make a copy of all input parameters potentially modified by Sysmat
  // call --- they are not intended to be modified
//  double copy_Cs         =Cs;
//  double copy_Cs_delta_sq=Cs_delta_sq;
//  double copy_l_tau      =l_tau;

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
  printf("| FINITE DIFFERENCE CHECK FOR ELEMENT %5d |\n",eid);
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

        if (f3Parameter_->is_genalpha_)
        {
          checkepreaf(nn)+=f3Parameter_->alphaF_*epsilon;
        }
        else
        {
          checkepreaf(nn)+=epsilon;
        }
      }
      else
      {
        printf("velocity dof %d (%d)\n",rr,nn);

        if (f3Parameter_->is_genalpha_)
        {
          checkevelaf(rr,nn)+=f3Parameter_->alphaF_*epsilon;
          checkeaccam(rr,nn)+=f3Parameter_->alphaM_/(f3Parameter_->gamma_*f3Parameter_->dt_)*epsilon;
        }
        else
        {
          checkevelaf(rr,nn)+=epsilon;
        }
      }

      // TODO: Andi
      // calculate the right hand side for the perturbed vector
//      Sysmat2D3D(checkevelaf,
//                 eveln,
//                 fsevelaf,
//                 checkepreaf,
//                 checkeaccam,
//                 escaaf,
//                 escaam,
//                 escadtam,
//                 emhist,
//                 edispnp,
//                 egridv,
//                 checkmat1,
//                 checkmat2,
//                 checkvec1,
//                 thermpressaf,
//                 thermpressam,
//                 thermpressdtaf,
//                 thermpressdtam,
//                 copy_material,
//                 timefac,
//                 copy_Cs,
//                 copy_Cs_delta_sq,
//                 copy_l_tau);

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
        if (f3Parameter_->is_genalpha_)
        {
          val   =-(eforce(mm)   /(epsilon))*(f3Parameter_->gamma_*f3Parameter_->dt_)/(f3Parameter_->alphaM_);
          lin   =-(eforce(mm)   /(epsilon))*(f3Parameter_->gamma_*f3Parameter_->dt_)/(f3Parameter_->alphaM_)+estif(mm,dof);
          nonlin=-(checkvec1(mm)/(epsilon))*(f3Parameter_->gamma_*f3Parameter_->dt_)/(f3Parameter_->alphaM_);
        }
        else
        {
          val   =-eforce(mm)/epsilon;
          lin   =-eforce(mm)/epsilon+estif(mm,dof);
          nonlin=-checkvec1(mm)/epsilon;
        }

        double norm=abs(lin);
        if(norm<1e-12)
        {
          norm=1e-12;
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


// tenplate classes
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex8>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex20>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex27>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet4>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet10>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::wedge6>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::pyramid5>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad4>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad8>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad9>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tri3>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tri6>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::nurbs9>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::nurbs27>;
