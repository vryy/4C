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
#include "fluid_ele_utils.H"

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
    DRT::ELEMENTS::Fluid*    ele,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm            ,
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


/*---------------------------------------------------------------------*
 | Action type: calc_divop                                             |
 | calculate integrated divergence operator              mayr.mt 04/12 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::CalcDivOp(
    DRT::ELEMENTS::Fluid*     ele,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm            ,
    Epetra_SerialDenseVector& elevec1       )
{
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  if (ele->IsAle()) // Do ALE specific updates if necessary
  {
    LINALG::Matrix<nsd_,nen_> edispnp(true);
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"dispnp");

    // get new node positions of ALE mesh
     xyze_ += edispnp;
  }

  // integration loop
  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad,ele->Id());

    for (int nodes = 0; nodes < nen_; nodes++) // loop over nodes
    {
      for (int dim = 0; dim < nsd_; dim++) // loop over spatial dimensions
      {
        elevec1((nsd_+1) * nodes + dim) +=  derxy_(dim,nodes) * fac_;
      }
    }
  } // end of integration loop

  return 0;
}

/*----------------------------------------------------------------------*
 * Action type: Compute Error                              shahmiri 01/12
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::ComputeError(
    DRT::ELEMENTS::Fluid*           ele,
    Teuchos::ParameterList&         params,
    Teuchos::RCP<MAT::Material>&    mat,
    DRT::Discretization&            discretization,
    std::vector<int>&               lm,
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
    DRT::ELEMENTS::Fluid*           ele,
    Teuchos::ParameterList&         params,
    Teuchos::RCP<MAT::Material>&    mat,
    DRT::Discretization&            discretization,
    std::vector<int>&               lm,
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

    const double t = fldpara_->Time();

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

      // function evaluation requires a 3D position vector!!
      double position[3];
      position[0] = xyzint(0);
      position[1] = xyzint(1);
      position[2] = 0.0;

      if (1.0 < position[0] and position[0] < 2.0 and 0.0 < position[1] and position[1] < position[0])
      {
        const double u_exact_x = DRT::Problem::Instance()->Funct(0).Evaluate(0,position,t,NULL);
        const double u_exact_y = DRT::Problem::Instance()->Funct(0).Evaluate(1,position,t,NULL);
        u(0) = u_exact_x;
        u(1) = u_exact_y;
      }

    }
    break;
    case INPAR::FLUID::byfunct1:
    {
      const int func_no = 1;


      // function evaluation requires a 3D position vector!!
      double position[3];

      if (nsd_ == 2)
      {

        position[0] = xyzint(0);
        position[1] = xyzint(1);
        position[2] = 0.0;
      }
      else if(nsd_ == 3)
      {
        position[0] = xyzint(0);
        position[1] = xyzint(1);
        position[2] = xyzint(2);
      }
      else dserror("invalid nsd %d", nsd_);

      if(nsd_ == 2)
      {
        const double u_exact_x = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(0,position,t,NULL);
        const double u_exact_y = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(1,position,t,NULL);
        const double p_exact   = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(2,position,t,NULL);

        u(0) = u_exact_x;
        u(1) = u_exact_y;
        p    = p_exact;
      }
      else if(nsd_==3)
      {
        const double u_exact_x = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(0,position,t,NULL);
        const double u_exact_y = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(1,position,t,NULL);
        const double u_exact_z = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(2,position,t,NULL);
        const double p_exact   = DRT::Problem::Instance()->Funct(func_no-1).Evaluate(3,position,t,NULL);

        u(0) = u_exact_x;
        u(1) = u_exact_y;
        u(2) = u_exact_z;
        p    = p_exact;
      }
      else dserror("invalid dimension");

    }
    break;
    default:
      dserror("analytical solution is not defined");
      break;
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
                                    const std::vector<int>&      lm,             ///<
                                    FLD::RotationallySymmetricPeriodicBC<distype> & rotsymmpbc, ///<
                                    LINALG::Matrix<nsd_,nen_> *  matrixtofill,   ///< vector field
                                    LINALG::Matrix<nen_,1> *     vectortofill,   ///< scalar field
                                    const std::string            state)          ///< state of the global vector
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(state);
  if(matrix_state == Teuchos::null)
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
  Fluid*                     ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  std::vector<int>&          lm,
  RCP<MAT::Material> mat)
{
  //----------------------------------------------------------------------
  // get all nodal values
  // ---------------------------------------------------------------------

  // call routine for calculation of body force in element nodes,
  // with pressure gradient prescribed as body force included for turbulent
  // channel flow and with scatra body force included for variable-density flow
  // (evaluation at time n+alpha_F for generalized-alpha scheme,
  //  and at time n+1 otherwise)
  LINALG::Matrix<nsd_,nen_> ebofoaf(true);
  LINALG::Matrix<nsd_,nen_> eprescpgaf(true);
  LINALG::Matrix<nen_,1>    escabofoaf(true);
  BodyForce(ele,fldpara_,ebofoaf,eprescpgaf,escabofoaf);

  // if not available, the arrays for the subscale quantities have to be
  // resized and initialised to zero
  if (fldpara_->Tds()==INPAR::FLUID::subscales_time_dependent)
   dserror("Time-dependent subgrid scales not supported");

  // get all general state vectors: velocity/pressure, scalar,
  // acceleration/scalar time derivative and history
  // velocity/pressure and scalar values are at time n+alpha_F/n+alpha_M
  // for generalized-alpha scheme and at time n+1/n for all other schemes
  // acceleration/scalar time derivative values are at time n+alpha_M for
  // generalized-alpha scheme and at time n+1 for all other schemes
  // fill the local element vector/matrix with the global values
  // af_genalpha: velocity/pressure at time n+alpha_F
  // np_genalpha: velocity at time n+alpha_F, pressure at time n+1
  // ost:         velocity/pressure at time n+1
  LINALG::Matrix<nsd_,nen_> evelaf(true);
  LINALG::Matrix<nen_,1>    epreaf(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelaf, &epreaf,"velaf");

  // np_genalpha: additional vector for velocity at time n+1
  LINALG::Matrix<nsd_,nen_> evelnp(true);
  if (fldpara_->IsGenalphaNP())
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelnp, NULL,"velnp");

  LINALG::Matrix<nen_,1> escaaf(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, NULL, &escaaf,"scaaf");

  LINALG::Matrix<nsd_,nen_> emhist(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &emhist, NULL,"hist");

  LINALG::Matrix<nsd_,nen_> eaccam(true);
  LINALG::Matrix<nen_,1>    escadtam(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &eaccam, &escadtam,"accam");

  LINALG::Matrix<nsd_,nen_> eveln(true);
  LINALG::Matrix<nen_,1>    escaam(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &eveln, &escaam,"scaam");

  if (fldpara_->IsGenalpha()) eveln.Clear();
  else                            eaccam.Clear();

  // get additional state vectors for ALE case: grid displacement and vel.
  LINALG::Matrix<nsd_, nen_> edispnp(true);
  LINALG::Matrix<nsd_, nen_> egridv(true);

  if (ele->IsAle())
  {
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"dispnp");
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &egridv, NULL,"gridv");
  }

  // get additional state vector for AVM3 case and multifractal subgrid scales:
  // fine-scale velocity values are at time n+alpha_F for generalized-alpha
  // scheme and at time n+1 for all other schemes
  LINALG::Matrix<nsd_,nen_> fsevelaf(true);
  LINALG::Matrix<nen_,1>    fsescaaf(true);
  if (fldpara_->Fssgv() != INPAR::FLUID::no_fssgv
   or fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales)
  {
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &fsevelaf, NULL,"fsvelaf");
    if(fldpara_->PhysicalType() == INPAR::FLUID::loma and fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales)
     ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, NULL, &fsescaaf,"fsscaaf");
  }

  // get node coordinates and number of elements per node
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // for scale similarity model:
  // get filtered veolcities and reynoldsstresses
  LINALG::Matrix<nsd_,nen_> evel_hat(true);
  LINALG::Matrix<nsd_*nsd_,nen_> ereynoldsstress_hat(true);
  if (fldpara_->TurbModAction() == INPAR::FLUID::scale_similarity
      or fldpara_->TurbModAction() == INPAR::FLUID::scale_similarity_basic)
  {
    RCP<Epetra_MultiVector> filtered_vel = params.get<RCP<Epetra_MultiVector> >("Filtered velocity");
    RCP<Epetra_MultiVector> fs_vel = params.get<RCP<Epetra_MultiVector> >("Fine scale velocity");
    RCP<Epetra_MultiVector> filtered_reystre = params.get<RCP<Epetra_MultiVector> >("Filtered reynoldsstress");

    for (int nn=0;nn<nen_;++nn)
    {
      int lid = (ele->Nodes()[nn])->LID();

      for (int dimi=0;dimi<3;++dimi)
      {
        evel_hat(dimi,nn) = (*((*filtered_vel)(dimi)))[lid];
        fsevelaf(dimi,nn) = (*((*fs_vel)(dimi)))[lid];

        for (int dimj=0;dimj<3;++dimj)
        {
          int index=3*dimi+dimj;

          ereynoldsstress_hat(index,nn) = (*((*filtered_reystre)(index)))[lid];

        }
      }
    }
  }

  // flag for higher order elements
  is_higher_order_ele_ = IsHigherOrder<distype>::ishigherorder;
  // overrule higher_order_ele if input-parameter is set
  // this might be interesting for fast (but slightly
  // less accurate) computations
  if (fldpara_->IsInconsistent() == true) is_higher_order_ele_ = false;

  // set thermodynamic pressure at n+1/n+alpha_F and n+alpha_M/n and
  // its time derivative at n+alpha_M/n+1
  const double thermpressaf   = params.get<double>("thermpress at n+alpha_F/n+1");
  const double thermpressam   = params.get<double>("thermpress at n+alpha_M/n");
  const double thermpressdtaf = params.get<double>("thermpressderiv at n+alpha_F/n+1");
  const double thermpressdtam = params.get<double>("thermpressderiv at n+alpha_M/n+1");

  // ---------------------------------------------------------------------
  // set parameters for classical turbulence models
  // ---------------------------------------------------------------------
  Teuchos::ParameterList& turbmodelparams = params.sublist("TURBULENCE MODEL");

  double Ci_delta_sq = 0.0;
  double Cs_delta_sq = 0.0;
  visceff_ = 0.0;

  // remember the layer of averaging for the dynamic Smagorinsky model
  int  smaglayer=0;

  double CsDeltaSq = 0.0;
  double CiDeltaSq = 0.0;
  if (fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky)
  {
    RCP<Epetra_Vector> ele_CsDeltaSq = params.sublist("TURBULENCE MODEL").get<RCP<Epetra_Vector> >("col_Cs_delta_sq");
    RCP<Epetra_Vector> ele_CiDeltaSq = params.sublist("TURBULENCE MODEL").get<RCP<Epetra_Vector> >("col_Ci_delta_sq");
    const int id = ele->LID();
    CsDeltaSq = (*ele_CsDeltaSq)[id];
    CiDeltaSq = (*ele_CiDeltaSq)[id];
  }
  GetTurbulenceParams(turbmodelparams,
                      Cs_delta_sq,
                      Ci_delta_sq,
                      smaglayer,
                      CsDeltaSq,
                      CiDeltaSq);


  //----------------------------------------------------------------------
  // prepare mean values
  // ---------------------------------------------------------------------

  // the coordinates of the element layers in the channel
  // planecoords are named nodeplanes in turbulence_statistics_channel!
  RCP<std::vector<double> > planecoords  = params.get<RCP<std::vector<double> > >("planecoords_",Teuchos::null);
  if(planecoords==Teuchos::null)
    dserror("planecoords is null, but need channel_flow_of_height_2\n");

  //this will be the y-coordinate of a point in the element interior
  double center = 0.0;
  // get node coordinates of element
  for(int inode=0;inode<nen_;inode++)
    center+=xyze_(1,inode);

  center/=nen_;

  // working arrays for the quantities we want to compute
  LINALG::Matrix<nsd_,1>  mean_res        ;
  LINALG::Matrix<nsd_,1>  mean_sacc       ;
  LINALG::Matrix<nsd_,1>  mean_svelaf     ;
  LINALG::Matrix<nsd_,1>  mean_res_sq     ;
  LINALG::Matrix<nsd_,1>  mean_sacc_sq    ;
  LINALG::Matrix<nsd_,1>  mean_svelaf_sq  ;
  LINALG::Matrix<nsd_,1>  mean_tauinvsvel ;

  LINALG::Matrix<2*nsd_,1>  mean_crossstress;
  LINALG::Matrix<2*nsd_,1>  mean_reystress  ;

  double vol             = 0.0;

  double h               = 0.0;
  double h_bazilevs      = 0.0;
  double strle           = 0.0;
  double gradle          = 0.0;
  double averaged_tauC   = 0.0;
  double averaged_tauM   = 0.0;

  double abs_res         = 0.0;
  double abs_svel        = 0.0;

  double mean_resC       = 0.0;
  double mean_resC_sq    = 0.0;
  double mean_sprenp     = 0.0;
  double mean_sprenp_sq  = 0.0;

  double eps_visc        = 0.0;
  double eps_conv        = 0.0;
  double eps_smag        = 0.0;
  double eps_avm3        = 0.0;
  double eps_mfs         = 0.0;
  double eps_mfscross    = 0.0;
  double eps_mfsrey      = 0.0;
  double eps_supg        = 0.0;
  double eps_cross       = 0.0;
  double eps_rey         = 0.0;
  double eps_cstab       = 0.0;
  double eps_pspg        = 0.0;

  mean_res        .Clear();
  mean_sacc       .Clear();
  mean_svelaf     .Clear();
  mean_res_sq     .Clear();
  mean_sacc_sq    .Clear();
  mean_svelaf_sq  .Clear();
  mean_tauinvsvel .Clear();
  mean_crossstress.Clear();
  mean_reystress  .Clear();


  // ---------------------------------------------------------------------
  // calculate volume and evaluate material, tau ... at element center
  // ---------------------------------------------------------------------

  // evaluate shape functions and derivatives at element center
  EvalShapeFuncAndDerivsAtEleCenter(ele->Id());

  // set element area or volume
  vol = fac_;

  //------------------------------------------------------------------------
  // potential evaluation of material parameters, subgrid viscosity
  // and/or stabilization parameters at element center
  //------------------------------------------------------------------------
  // get material parameters at element center
  if (not fldpara_->MatGp() or not fldpara_->TauGp())
  {
    GetMaterialParams(mat,evelaf,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam);

    // calculate all-scale or fine-scale subgrid viscosity at element center
    visceff_ = visc_;
    if (fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky)
    {
      CalcSubgrVisc(evelaf,vol,fldpara_->Cs_,Cs_delta_sq,Ci_delta_sq,fldpara_->l_tau_);
      // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
      visceff_ += sgvisc_;
    }
    else if (fldpara_->Fssgv() != INPAR::FLUID::no_fssgv)
      CalcFineScaleSubgrVisc(evelaf,fsevelaf,vol,fldpara_->Cs_);
  }

  // potential evaluation of multifractal subgrid-scales at element center
  // coefficient B of fine-scale velocity
  LINALG::Matrix<nsd_,1> B_mfs(true);
  // coefficient D of fine-scale scalar (loma only)
  double D_mfs = 0.0;
  if (fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales)
  {
    if (not fldpara_->BGp())
    {
      // make sure to get material parameters at element center
      if (fldpara_->MatGp())
        //GetMaterialParams(material,evelaf,escaaf,escaam,thermpressaf,thermpressam,thermpressdtam);
        GetMaterialParams(mat,evelaf,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam);

      // provide necessary velocities and gradients at element center
      velint_.Multiply(evelaf,funct_);
      fsvelint_.Multiply(fsevelaf,funct_);
      vderxy_.MultiplyNT(evelaf,derxy_);
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
  if (not fldpara_->TauGp())
  {
    // get convective velocity at element center for evaluation of
    // stabilization parameter
    velint_.Multiply(evelaf,funct_);
    convvelint_.Update(velint_);
    if (ele->IsAle()) convvelint_.Multiply(-1.0,egridv,funct_,1.0);

    // calculate stabilization parameters at element center
    CalcStabParameter(vol);
  }


  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
  {
    //---------------------------------------------------------------
    // evaluate shape functions and derivatives at integration point
    //---------------------------------------------------------------
    EvalShapeFuncAndDerivsAtIntPoint(iquad,ele->Id());

    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    velint_.Multiply(evelaf,funct_);

    // get velocity derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    vderxy_.MultiplyNT(evelaf,derxy_);

    // get fine-scale velocity and its derivatives at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    if (fldpara_->Fssgv() != INPAR::FLUID::no_fssgv)
    {
      fsvderxy_.MultiplyNT(fsevelaf,derxy_);
    }
    else
    {
      fsvderxy_.Clear();
    }
    if (fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales)
    {
      fsvelint_.Multiply(fsevelaf,funct_);
      fsvderxy_.MultiplyNT(fsevelaf,derxy_);
    }
    else
    {
      fsvelint_.Clear();
    }

    // get convective velocity at integration point
    // (ALE case handled implicitly here using the (potential
    //  mesh-movement-dependent) convective velocity, avoiding
    //  various ALE terms used to be calculated before)
    convvelint_.Update(velint_);
    if (ele->IsAle())
    {
      gridvelint_.Multiply(egridv,funct_);
      convvelint_.Update(-1.0,gridvelint_,1.0);
    }

    // get pressure gradient at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    gradp_.Multiply(derxy_,epreaf);

    // get bodyforce at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    bodyforce_.Multiply(ebofoaf,funct_);
    // get prescribed pressure gradient acting as body force
    // (required for turbulent channel flow)
    prescribedpgrad_.Multiply(eprescpgaf,funct_);

    // get momentum history data at integration point
    // (only required for one-step-theta and BDF2 time-integration schemes)
    histmom_.Multiply(emhist,funct_);

//    if(fldpara_->TurbModAction() == INPAR::FLUID::scale_similarity_basic)
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

    // evaluation of various partial operators at integration point
    // compute convective term from previous iteration and convective operator
    conv_old_.Multiply(vderxy_,convvelint_);

    // compute viscous term from previous iteration and viscous operator
    if (is_higher_order_ele_) CalcDivEps(evelaf);
    else
      visc_old_.Clear();

    // compute divergence of velocity from previous iteration
    vdiv_ = 0.0;
    if (not fldpara_->IsGenalphaNP())
    {
      for (int idim = 0; idim <nsd_; ++idim)
      {
        vdiv_ += vderxy_(idim, idim);
      }
    }
    else
    {
      for (int idim = 0; idim <nsd_; ++idim)
      {
        //get vdiv at time n+1 for np_genalpha,
        LINALG::Matrix<nsd_,nsd_> vderxy(true);
        vderxy.MultiplyNT(evelnp,derxy_);
        vdiv_ += vderxy(idim, idim);
      }
    }

    // get material parameters at integration point
    if (fldpara_->MatGp())
    {
      GetMaterialParams(mat,evelaf,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam);

      // calculate all-scale or fine-scale subgrid viscosity at integration point
      visceff_ = visc_;
      if (fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky)
      {
        CalcSubgrVisc(evelaf,vol,fldpara_->Cs_,Cs_delta_sq,Ci_delta_sq,fldpara_->l_tau_);
        // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
        visceff_ += sgvisc_;
      }
      else if (fldpara_->Fssgv() != INPAR::FLUID::no_fssgv)
        CalcFineScaleSubgrVisc(evelaf,fsevelaf,vol,fldpara_->Cs_);
    }

    // potential evaluation of coefficient of multifractal subgrid-scales at integration point
    if (fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales)
    {
      if (fldpara_->BGp())
      {
        // make sure to get material parameters at gauss point
        if (not fldpara_->MatGp())
          GetMaterialParams(mat,evelaf,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam);

        // calculate parameters of multifractal subgrid-scales
        PrepareMultifractalSubgrScales(B_mfs, D_mfs, evelaf, fsevelaf, vol);
      }

      // calculate fine-scale velocity, its derivative and divergence for multifractal subgrid-scale modeling
      for (int idim=0; idim<nsd_; idim++)
        mffsvelint_(idim,0) = fsvelint_(idim,0) * B_mfs(idim,0);
    }
    else
    {
      mffsvelint_.Clear();
      mffsvderxy_.Clear();
      mffsvdiv_ = 0.0;
    }

    // calculate stabilization parameter at integration point
    if (fldpara_->TauGp())
      CalcStabParameter(vol);


    // compute residual of continuity equation
    // residual contains velocity divergence only for incompressible flow
    conres_old_ = vdiv_;

    // following computations only required for variable-density flow at low Mach number
    if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
    {
      // compute additional Galerkin terms on right-hand side of continuity equation
      // -> different for generalized-alpha and other time-integration schemes
      ComputeGalRHSContEq(eveln,escaaf,escaam,escadtam,ele->IsAle());

      // add to residual of continuity equation
      conres_old_ -= rhscon_;

      // compute subgrid-scale part of scalar
      // -> different for generalized-alpha and other time-integration schemes
      ComputeSubgridScaleScalar(escaaf,escaam);

      // update material parameters including subgrid-scale part of scalar
      if (fldpara_->UpdateMat())
      {
        if (fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales)
          UpdateMaterialParams(mat,evelaf,escaaf,escaam,thermpressaf,thermpressam,mfssgscaint_);
        else
          UpdateMaterialParams(mat,evelaf,escaaf,escaam,thermpressaf,thermpressam,sgscaint_);
        visceff_ = visc_;
        if (fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky)
        visceff_ += sgvisc_;
      }
    }

    // evaluate momentum residual once for all stabilization right hand sides
    if (fldpara_->IsGenalpha())
    {
      // get acceleration at time n+alpha_M at integration point
      accint_.Multiply(eaccam,funct_);

      for (int rr=0;rr<nsd_;++rr)
      {
        momres_old_(rr) = densam_*accint_(rr)+densaf_*conv_old_(rr)+gradp_(rr)
                       -2*visceff_*visc_old_(rr)-densaf_*bodyforce_(rr)-prescribedpgrad_(rr);
      }
    }
    else
    {
      rhsmom_.Update((densn_/fldpara_->Dt()/fldpara_->Theta()),histmom_,densaf_,bodyforce_);
      // and pressure gradient prescribed as body force
      // caution: not density weighted
      rhsmom_.Update(1.0,prescribedpgrad_,1.0);
      // compute instationary momentum residual:
      // momres_old = u_(n+1)/dt + theta ( ... ) - histmom_/dt - theta*bodyforce_
      for (int rr=0;rr<nsd_;++rr)
      {
        momres_old_(rr) = ((densaf_*velint_(rr)/fldpara_->Dt()
                         +fldpara_->Theta()*(densaf_*conv_old_(rr)+gradp_(rr)
                         -2*visceff_*visc_old_(rr)))/fldpara_->Theta())-rhsmom_(rr);
      }
    }


    //---------------------------------------------------------------
    // element average dissipation and production rates
    //---------------------------------------------------------------

    //---------------------------------------------------------------
    // residual-based subgrid-scale modeling terms
    //---------------------------------------------------------------

    // dissipation by supg-stabilization
    if (fldpara_->SUPG() == INPAR::FLUID::convective_stab_supg)
    {
      for (int rr=0;rr<nsd_;rr++)
      {
        eps_supg += densaf_ * fac_ * tau_(0) * momres_old_(rr,0) * conv_old_(rr,0);
      }
    }

    // dissipation by cross-stress-stabilization
    if (fldpara_->Cross() != INPAR::FLUID::cross_stress_stab_none)
    {
      for (int rr=0;rr<nsd_;rr++)
      {
        eps_cross += densaf_ * fac_ * tau_(0) * velint_(rr,0) * ( momres_old_(0,0) * vderxy_ (rr,0)
                                                                + momres_old_(1,0) * vderxy_ (rr,1)
                                                                + momres_old_(2,0) * vderxy_ (rr,2));
      }
    }

    // dissipation by reynolds-stress-stabilization
    if (fldpara_->Reynolds() != INPAR::FLUID::reynolds_stress_stab_none)
    {
      for (int rr=0;rr<nsd_;rr++)
      {
        eps_rey -= densaf_ * fac_ * tau_(0) * tau_(0) * momres_old_(rr,0) * ( momres_old_(0,0) * vderxy_ (rr,0)
                                                                            + momres_old_(1,0) * vderxy_ (rr,1)
                                                                            + momres_old_(2,0) * vderxy_ (rr,2));
      }
    }

    // dissipation by pspg-stabilization
    if (fldpara_->PSPG() == INPAR::FLUID::pstab_use_pspg)
    {
      for (int rr=0;rr<nsd_;rr++)
      {
        eps_pspg += fac_ * gradp_(rr,0) * tau_(1) * momres_old_(rr,0);
      }
    }

    // dissipation by continuity-stabilization
    if (fldpara_->CStab() == INPAR::FLUID::continuity_stab_yes)
    {
      eps_cstab += fac_ * vdiv_ * tau_(2) * conres_old_;
    }

    //---------------------------------------------------------------
    // multifractal subgrid-scale modeling terms
    //---------------------------------------------------------------

    // dissipation multifractal subgrid-scales
    if(fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales)
    {
      for (int rr=0;rr<nsd_;rr++)
      {
        eps_mfs -= densaf_ * fac_ * ( mffsvelint_(rr,0) * conv_old_(rr,0)
                                    + velint_(rr,0) * ( mffsvelint_(0,0) * vderxy_ (rr,0)
                                                      + mffsvelint_(1,0) * vderxy_ (rr,1)
                                                      + mffsvelint_(2,0) * vderxy_ (rr,2))
                                    + mffsvelint_(rr,0) * ( mffsvelint_(0,0) * vderxy_ (rr,0)
                                                          + mffsvelint_(1,0) * vderxy_ (rr,1)
                                                          + mffsvelint_(2,0) * vderxy_ (rr,2)));;

        eps_mfscross -= densaf_ * fac_ * ( mffsvelint_(rr,0) * conv_old_(rr,0)
                                         + velint_(rr,0) * ( mffsvelint_(0,0) * vderxy_ (rr,0)
                                                           + mffsvelint_(1,0) * vderxy_ (rr,1)
                                                           + mffsvelint_(2,0) * vderxy_ (rr,2)));

        eps_mfsrey -= densaf_ * fac_ * mffsvelint_(rr,0) * ( mffsvelint_(0,0) * vderxy_ (rr,0)
                                                           + mffsvelint_(1,0) * vderxy_ (rr,1)
                                                           + mffsvelint_(2,0) * vderxy_ (rr,2));
      }
    }

    //---------------------------------------------------------------
    // small-scale subgrid-viscosity subgrid-scale modeling terms
    //---------------------------------------------------------------

    // dissipation AVM3
    /*
                         /                                \
                        |       /  n+1 \         / n+1 \   |
          2* visc    *  |  eps | du     | , eps | u     |  |
                 turb   |       \      /         \     /   |
                         \                                /
    */
    if (fldpara_->Fssgv() != INPAR::FLUID::no_fssgv)
    {
      LINALG::Matrix<nsd_,nsd_> fstwo_epsilon;
      for(int rr=0;rr<nsd_;++rr)
      {
        for(int mm=0;mm<nsd_;++mm)
        {
          fstwo_epsilon(rr,mm) = fsvderxy_(rr,mm) + fsvderxy_(mm,rr);
        }
      }
      for(int rr=0;rr<nsd_;++rr)
      {
        for(int mm=0;mm<nsd_;++mm)
        {
//          eps_avm3 += 0.5*fssgvisc_*fac_*fstwo_epsilon(rr,mm)*two_epsilon(rr,mm);
          eps_avm3 += 0.5*fssgvisc_*fac_*fstwo_epsilon(rr,mm)*fstwo_epsilon(rr,mm);
        }
      }
    }

    //---------------------------------------------------------------
    // Smagorinsky model
    //---------------------------------------------------------------

    // dissipation (Smagorinsky)
    /*
                         /                                \
                        |       / n+1 \         / n+1 \   |
          2* visc    *  |  eps | u     | , eps | u     |  |
                 turb   |       \     /         \     /   |
                         \                                /
    */
    LINALG::Matrix<nsd_,nsd_> two_epsilon;
    for(int rr=0;rr<nsd_;++rr)
    {
      for(int mm=0;mm<nsd_;++mm)
      {
        two_epsilon(rr,mm) = vderxy_(rr,mm) + vderxy_(mm,rr);
      }
    }
    if(fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky
      or fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky)
    {
      for(int rr=0;rr<nsd_;++rr)
      {
        for(int mm=0;mm<nsd_;++mm)
        {
          eps_smag += 0.5*sgvisc_*fac_*two_epsilon(rr,mm)*two_epsilon(rr,mm);
        }
      }
      if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
        eps_smag -= (2.0/3.0)*fac_*(sgvisc_*vdiv_+q_sq_)*vdiv_;
    }
    

    //---------------------------------------------------------------
    // scale-similarity model
    //---------------------------------------------------------------

    // dissipation Scale Similarity
    /*
             /                                \
            |   ssm  /^n+1 \         / n+1 \   |
            |  tau  | u     | , eps | u     |  |
            |        \     /         \     /   |
             \                                /
    */
//    if(fldpara_->TurbModAction() == INPAR::FLUID::scale_similarity_basic)
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
//          eps_scsim += -0.5*fac_*densaf_*fldpara_->Cl()*tau_scale_sim(rr,mm)*two_epsilon(mm,rr);
//        }
//      }
//      if (Production >= 0.0)
//      {
//        // forwardscatter
//        for(int rr=0;rr<nsd_;++rr)
//        {
//          for(int mm=0;mm<nsd_;++mm)
//          {
//            eps_scsimfs += -0.5*fac_*densaf_*fldpara_->Cl()*tau_scale_sim(rr,mm)*two_epsilon(mm,rr);
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
//            eps_scsimbs += -0.5*fac_*densaf_*fldpara_->Cl()*tau_scale_sim(rr,mm)*two_epsilon(mm,rr);
//          }
//        }
//      }
//    }

    //---------------------------------------------------------------
    // standard Galerkin terms
    //---------------------------------------------------------------

    // convective (Galerkin)
    /*
                 /                          \
                |   n+1   / n+1 \   /  n+1\  |
                |  u    , | u   | o | u   |  |
                |         \     /   \     /  |
                 \                          /
    */
    for (int rr=0;rr<nsd_;rr++)
    {
      eps_conv -= densaf_ * fac_ * velint_(rr,0) * ( velint_(0,0) * vderxy_ (rr,0)
                                                   + velint_(1,0) * vderxy_ (rr,1)
                                                   + velint_(2,0) * vderxy_ (rr,2));
    }

    // dissipation (Galerkin)
    /*
                     /                                \
                    |       / n+1 \         / n+1 \   |
          2* visc * |  eps | u     | , eps | u     |  |
                    |       \     /         \     /   |
                     \                                /
    */
    for(int rr=0;rr<nsd_;++rr)
    {
      for(int mm=0;mm<nsd_;++mm)
      {
        eps_visc += 0.5*visc_*fac_*two_epsilon(rr,mm)*two_epsilon(rr,mm);
      }
    }
    if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
      eps_visc -= (2.0/3.0)*visc_*fac_*vdiv_*vdiv_;


    //---------------------------------------------------------------
    // reference length for stabilization parameters
    //---------------------------------------------------------------
    // volume based element size
    double hk = std::pow((6.*vol/M_PI),(1.0/3.0))/sqrt(3.0);
    h += fac_*hk;

    // streamlength based element size
    {
      const double vel_norm=velint_.Norm2();

      // this copy of velintaf_ will be used to store the normed velocity
      LINALG::Matrix<3,1> normed_velint;

      // normed velocity at element center (we use the copy for safety reasons!)
      if (vel_norm >= 1e-6)
      {
        for (int rr=0;rr<3;++rr) /* loop element nodes */
        {
          normed_velint(rr)=velint_(rr)/vel_norm;
        }
      }
      else
      {
        normed_velint(0) = 1.;
        for (int rr=1;rr<3;++rr) /* loop element nodes */
        {
          normed_velint(rr)=0.0;
        }
      }

      // get streamlength
      double val = 0.0;
      for (int rr=0;rr<nen_;++rr) /* loop element nodes */
      {
        val += fabs( normed_velint(0)*derxy_(0,rr)
                    +normed_velint(1)*derxy_(1,rr)
                    +normed_velint(2)*derxy_(2,rr));
      } /* end of loop over element nodes */
      strle += 2.0/val*fac_;
    }

    // element size in main gradient direction
    {
      // this copy of velintaf_ will be used to store the normed velocity
      LINALG::Matrix<3,1> normed_velgrad;

      for (int rr=0;rr<3;++rr)
      {
        normed_velgrad(rr)=sqrt(vderxy_(0,rr)*vderxy_(0,rr)
                                +
                                vderxy_(1,rr)*vderxy_(1,rr)
                                +
                                vderxy_(2,rr)*vderxy_(2,rr));
      }
      double norm=normed_velgrad.Norm2();

      // normed gradient
      if (norm>1e-6)
      {
        for (int rr=0;rr<3;++rr)
        {
          normed_velgrad(rr)/=norm;
        }
      }
      else
      {
        normed_velgrad(0) = 1.;
        for (int rr=1;rr<3;++rr)
        {
          normed_velgrad(rr)=0.0;
        }
      }

      // get length in this direction
      double val = 0.0;
      for (int rr=0;rr<nen_;++rr) /* loop element nodes */
      {
        val += fabs( normed_velgrad(0)*derxy_(0,rr)
                    +normed_velgrad(1)*derxy_(1,rr)
                    +normed_velgrad(2)*derxy_(2,rr));
      } /* end of loop over element nodes */
      gradle += 2.0/val*fac_;
    }

    {
      /*          +-           -+   +-           -+   +-           -+
                  |             |   |             |   |             |
                  |  dr    dr   |   |  ds    ds   |   |  dt    dt   |
            G   = |  --- * ---  | + |  --- * ---  | + |  --- * ---  |
             ij   |  dx    dx   |   |  dx    dx   |   |  dx    dx   |
                  |    i     j  |   |    i     j  |   |    i     j  |
                  +-           -+   +-           -+   +-           -+
      */
      LINALG::Matrix<3,3> G;

      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          G(nn,rr) = xji_(nn,0)*xji_(rr,0);
          for (int mm=1;mm<3;++mm)
          {
            G(nn,rr) += xji_(nn,mm)*xji_(rr,mm);
          }
        }
      }

      /*          +----
                   \
          G : G =   +   G   * G
          -   -    /     ij    ij
          -   -   +----
                   i,j
      */
      double normG = 0;
      for (int nn=0;nn<3;++nn)
      {
        for (int rr=0;rr<3;++rr)
        {
          normG+=G(nn,rr)*G(nn,rr);
        }
      }

      h_bazilevs+=1./sqrt(sqrt(normG))*fac_;
     }


    //---------------------------------------------------------------
    // element averages of residual and subgrid scales
    //---------------------------------------------------------------
    for(int rr=0;rr<3;++rr)
    {
      mean_res    (rr) += momres_old_(rr)*fac_;
      mean_res_sq (rr) += momres_old_(rr)*momres_old_(rr)*fac_;
    }
    abs_res    += sqrt(momres_old_(0)*momres_old_(0)+momres_old_(1)*momres_old_(1)+momres_old_(2)*momres_old_(2))*fac_;

    for(int rr=0;rr<3;++rr)
    {
      const double aux = tau_(0)*momres_old_(rr);

      mean_svelaf   (rr) -= aux*fac_;
      mean_svelaf_sq(rr) += aux*aux*fac_;
    }

    abs_svel +=sqrt( momres_old_(0)*momres_old_(0)
                   + momres_old_(1)*momres_old_(1)
                   + momres_old_(2)*momres_old_(2))*tau_(0)*fac_;

    for(int rr=0;rr<3;++rr)
    {
      mean_tauinvsvel(rr)+=mean_svelaf(rr)/tau_(0);
    }


    {
      const double aux = tau_(2)*conres_old_;

      mean_sprenp     -= aux*fac_;
      mean_sprenp_sq  += aux*aux*fac_;
    }


    //---------------------------------------------------------------
    // element averages of cross stresses and cross stresses
    //---------------------------------------------------------------
    if(fldpara_->Cross() != INPAR::FLUID::cross_stress_stab_none)
    {
      mean_crossstress(0)+=fac_*tau_(0)*(momres_old_(0)*velint_(0)+velint_(0)*momres_old_(0));
      mean_crossstress(1)+=fac_*tau_(0)*(momres_old_(1)*velint_(1)+velint_(1)*momres_old_(1));
      mean_crossstress(2)+=fac_*tau_(0)*(momres_old_(2)*velint_(2)+velint_(2)*momres_old_(2));
      mean_crossstress(3)+=fac_*tau_(0)*(momres_old_(0)*velint_(1)+velint_(0)*momres_old_(1));
      mean_crossstress(4)+=fac_*tau_(0)*(momres_old_(1)*velint_(2)+velint_(1)*momres_old_(2));
      mean_crossstress(5)+=fac_*tau_(0)*(momres_old_(2)*velint_(0)+velint_(2)*momres_old_(0));
    }

    if(fldpara_->Reynolds() != INPAR::FLUID::reynolds_stress_stab_none)
    {
      mean_reystress(0)  -=fac_*tau_(0)*tau_(0)*(momres_old_(0)*momres_old_(0)+momres_old_(0)*momres_old_(0));
      mean_reystress(1)  -=fac_*tau_(0)*tau_(0)*(momres_old_(1)*momres_old_(1)+momres_old_(1)*momres_old_(1));
      mean_reystress(2)  -=fac_*tau_(0)*tau_(0)*(momres_old_(2)*momres_old_(2)+momres_old_(2)*momres_old_(2));
      mean_reystress(3)  -=fac_*tau_(0)*tau_(0)*(momres_old_(0)*momres_old_(1)+momres_old_(1)*momres_old_(0));
      mean_reystress(4)  -=fac_*tau_(0)*tau_(0)*(momres_old_(1)*momres_old_(2)+momres_old_(2)*momres_old_(1));
      mean_reystress(5)  -=fac_*tau_(0)*tau_(0)*(momres_old_(2)*momres_old_(0)+momres_old_(0)*momres_old_(2));
      }


    //---------------------------------------------------------------
    // element averages of tau_Mu and tau_C
    //---------------------------------------------------------------
    averaged_tauM+=tau_(0)*fac_;
    averaged_tauC+=tau_(2)*fac_;

    mean_resC    += conres_old_*fac_;
    mean_resC_sq += conres_old_*conres_old_*fac_;
  }// end integration loop


  for(int rr=0;rr<3;++rr)
  {
    mean_res        (rr)/= vol;
    mean_res_sq     (rr)/= vol;
    mean_sacc       (rr)/= vol;
    mean_sacc_sq    (rr)/= vol;
    mean_svelaf     (rr)/= vol;
    mean_svelaf_sq  (rr)/= vol;
    mean_tauinvsvel (rr)/= vol;
  }


  for(int rr=0;rr<6;++rr)
  {
    mean_crossstress(rr)/=vol;
    mean_reystress  (rr)/=vol;
  }

  abs_res         /= vol;
  abs_svel        /= vol;

  mean_resC       /= vol;
  mean_resC_sq    /= vol;
  mean_sprenp     /= vol;
  mean_sprenp_sq  /= vol;

  h               /= vol;
  h_bazilevs      /= vol;
  strle           /= vol;
  gradle          /= vol;

  averaged_tauC   /= vol;
  averaged_tauM   /= vol;

  eps_visc /= vol;
  eps_conv /= vol;
  eps_smag /= vol;
  eps_avm3 /= vol;
  eps_mfs /= vol;
  eps_mfscross /= vol;
  eps_mfsrey /= vol;
  eps_supg /= vol;
  eps_cross /= vol;
  eps_rey /= vol;
  eps_cstab /= vol;
  eps_pspg /= vol;

  RCP<std::vector<double> > incrvol           = params.get<RCP<std::vector<double> > >("incrvol"          );

  RCP<std::vector<double> > incr_eps_visc      = params.get<RCP<std::vector<double> > >("incr_eps_visc"    );
  RCP<std::vector<double> > incr_eps_conv      = params.get<RCP<std::vector<double> > >("incr_eps_conv"    );
  RCP<std::vector<double> > incr_eps_smag      = params.get<RCP<std::vector<double> > >("incr_eps_eddyvisc");
  RCP<std::vector<double> > incr_eps_avm3      = params.get<RCP<std::vector<double> > >("incr_eps_avm3"    );
  RCP<std::vector<double> > incr_eps_mfs       = params.get<RCP<std::vector<double> > >("incr_eps_mfs"     );
  RCP<std::vector<double> > incr_eps_mfscross  = params.get<RCP<std::vector<double> > >("incr_eps_mfscross");
  RCP<std::vector<double> > incr_eps_mfsrey    = params.get<RCP<std::vector<double> > >("incr_eps_mfsrey"  );
  RCP<std::vector<double> > incr_eps_supg      = params.get<RCP<std::vector<double> > >("incr_eps_supg"    );
  RCP<std::vector<double> > incr_eps_cross     = params.get<RCP<std::vector<double> > >("incr_eps_cross"   );
  RCP<std::vector<double> > incr_eps_rey       = params.get<RCP<std::vector<double> > >("incr_eps_rey"     );
  RCP<std::vector<double> > incr_eps_cstab     = params.get<RCP<std::vector<double> > >("incr_eps_cstab"   );
  RCP<std::vector<double> > incr_eps_pspg      = params.get<RCP<std::vector<double> > >("incr_eps_pspg"    );

  RCP<std::vector<double> > incrhk            = params.get<RCP<std::vector<double> > >("incrhk"           );
  RCP<std::vector<double> > incrhbazilevs     = params.get<RCP<std::vector<double> > >("incrhbazilevs"    );
  RCP<std::vector<double> > incrstrle         = params.get<RCP<std::vector<double> > >("incrstrle"        );
  RCP<std::vector<double> > incrgradle        = params.get<RCP<std::vector<double> > >("incrgradle"       );

  RCP<std::vector<double> > incrres           = params.get<RCP<std::vector<double> > >("incrres"          );
  RCP<std::vector<double> > incrres_sq        = params.get<RCP<std::vector<double> > >("incrres_sq"       );
  RCP<std::vector<double> > incrabsres        = params.get<RCP<std::vector<double> > >("incrabsres"       );
  RCP<std::vector<double> > incrtauinvsvel    = params.get<RCP<std::vector<double> > >("incrtauinvsvel"   );

  RCP<std::vector<double> > incrsvelaf        = params.get<RCP<std::vector<double> > >("incrsvelaf"       );
  RCP<std::vector<double> > incrsvelaf_sq     = params.get<RCP<std::vector<double> > >("incrsvelaf_sq"    );
  RCP<std::vector<double> > incrabssvelaf     = params.get<RCP<std::vector<double> > >("incrabssvelaf"    );

  RCP<std::vector<double> > incrresC          = params.get<RCP<std::vector<double> > >("incrresC"         );
  RCP<std::vector<double> > incrresC_sq       = params.get<RCP<std::vector<double> > >("incrresC_sq"      );
  RCP<std::vector<double> > spressnp          = params.get<RCP<std::vector<double> > >("incrspressnp"     );
  RCP<std::vector<double> > spressnp_sq       = params.get<RCP<std::vector<double> > >("incrspressnp_sq"  );

  RCP<std::vector<double> > incrtauC          = params.get<RCP<std::vector<double> > >("incrtauC"         );
  RCP<std::vector<double> > incrtauM          = params.get<RCP<std::vector<double> > >("incrtauM"         );

  RCP<std::vector<double> > incrcrossstress   = params.get<RCP<std::vector<double> > >("incrcrossstress"  );
  RCP<std::vector<double> > incrreystress     = params.get<RCP<std::vector<double> > >("incrreystress"    );

  bool found = false;

  int nlayer = 0;
  for (nlayer=0;nlayer<(int)(*planecoords).size()-1;)
  {
    if(center<(*planecoords)[nlayer+1])
    {
      found = true;
      break;
    }
    nlayer++;
  }
  if (found ==false)
  {
    dserror("could not determine element layer");
  }

  // collect layer volume
  (*incrvol      )[nlayer] += vol;

  // element length in stabilisation parameter
  (*incrhk       )[nlayer] += h;

  // element length in viscous regime defined by the Bazilevs parameter
  (*incrhbazilevs)[nlayer] += h_bazilevs;

  // stream length
  (*incrstrle    )[nlayer] += strle;

  // gradient based element length
  (*incrgradle   )[nlayer] += gradle;

  // averages of stabilisation parameters
  (*incrtauC     )[nlayer] += averaged_tauC;
  (*incrtauM     )[nlayer] += averaged_tauM;

  // averages of momentum residuals, subscale velocity and accelerations
  for(int mm=0;mm<3;++mm)
  {
    (*incrres       )[3*nlayer+mm] += mean_res       (mm);
    (*incrres_sq    )[3*nlayer+mm] += mean_res_sq    (mm);

    (*incrsvelaf    )[3*nlayer+mm] += mean_svelaf    (mm);
    (*incrsvelaf_sq )[3*nlayer+mm] += mean_svelaf_sq (mm);

    (*incrtauinvsvel)[3*nlayer+mm] += mean_tauinvsvel(mm);
  }

  (*incrabsres       )[nlayer] += abs_res;
  (*incrabssvelaf    )[nlayer] += abs_svel;

  // averages of subscale pressure and continuity residuals
  (*incrresC         )[nlayer] += mean_resC      ;
  (*incrresC_sq      )[nlayer] += mean_resC_sq   ;

  (*spressnp         )[nlayer] += mean_sprenp    ;
  (*spressnp_sq      )[nlayer] += mean_sprenp_sq ;


  (*incr_eps_visc    )[nlayer] += eps_visc       ;
  (*incr_eps_conv    )[nlayer] += eps_conv       ;
  (*incr_eps_smag    )[nlayer] += eps_smag       ;
  (*incr_eps_avm3    )[nlayer] += eps_avm3       ;
  (*incr_eps_mfs     )[nlayer] += eps_mfs        ;
  (*incr_eps_mfscross)[nlayer] += eps_mfscross   ;
  (*incr_eps_mfsrey  )[nlayer] += eps_mfsrey     ;
  (*incr_eps_supg    )[nlayer] += eps_supg       ;
  (*incr_eps_cross   )[nlayer] += eps_cross      ;
  (*incr_eps_rey     )[nlayer] += eps_rey        ;
  (*incr_eps_cstab   )[nlayer] += eps_cstab      ;
  (*incr_eps_pspg    )[nlayer] += eps_pspg       ;

  // averages of subgrid stress tensors
  for(int mm=0;mm<6;++mm)
  {
    (*incrcrossstress)[6*nlayer+mm] += mean_crossstress(mm);
    (*incrreystress  )[6*nlayer+mm] += mean_reystress  (mm);
  }

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

        if (fldpara_->IsGenalpha())
        {
          checkepreaf(nn)+=fldpara_->AlphaF()*epsilon;
        }
        else
        {
          checkepreaf(nn)+=epsilon;
        }
      }
      else
      {
        printf("velocity dof %d (%d)\n",rr,nn);

        if (fldpara_->IsGenalpha())
        {
          checkevelaf(rr,nn)+=fldpara_->AlphaF()*epsilon;
          checkeaccam(rr,nn)+=fldpara_->AlphaM()/(fldpara_->Gamma()*fldpara_->Dt())*epsilon;
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
        if (fldpara_->IsGenalpha())
        {
          val   =-(eforce(mm)   /(epsilon))*(fldpara_->Gamma()*fldpara_->Dt())/(fldpara_->AlphaM());
          lin   =-(eforce(mm)   /(epsilon))*(fldpara_->Gamma()*fldpara_->Dt())/(fldpara_->AlphaM())+estif(mm,dof);
          nonlin=-(checkvec1(mm)/(epsilon))*(fldpara_->Gamma()*fldpara_->Dt())/(fldpara_->AlphaM());
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


// template classes
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
