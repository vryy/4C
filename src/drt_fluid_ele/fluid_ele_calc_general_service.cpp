/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_general_service.cpp

\brief general service routines for calculation of fluid element

<pre>
Maintainer: Ursula Rasthofer & Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236/-245
</pre>
*/
/*----------------------------------------------------------------------*/

#include "fluid_ele_factory.H"
#include "fluid_ele_calc.H"
#include "fluid_ele.H"
#include "fluid_ele_parameter.H"
#include "fluid_ele_parameter_timint.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "fluid_ele_action.H"

#include "../drt_fluid/fluid_rotsym_periodicbc.H"

#include "../drt_geometry/position_array.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_elementtype.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils.H"

#include "../drt_mat/newtonianfluid.H"

#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "Sacado.hpp"

/*----------------------------------------------------------------------*
 * Evaluate supporting methods of the element
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::EvaluateService(
    DRT::ELEMENTS::Fluid*     ele,
    Teuchos::ParameterList&   params,
    Teuchos::RCP<MAT::Material> & mat,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::Action act = DRT::INPUT::get<FLD::Action>(params,"action");

  switch(act)
  {
    case FLD::calc_div_u:
    {
      // compute divergence of velocity field at the element
      return ComputeDivU(ele, discretization, lm, elevec1);
    }
    break;
    case FLD::calc_fluid_error:
    {
      // compute error for a known analytical solution
      return ComputeError(ele, params, mat, discretization, lm, elevec1);
    }
    break;
    case FLD::calc_dissipation:
    {
      if (nsd_ == 3)
      {
        if (ele->Owner() == discretization.Comm().MyPID()) // don't store values of ghosted elements
        {
          return CalcDissipation(
              ele,
              params,
              discretization,
              lm,
              mat);
        }
      }
      else dserror("%i D elements does not support calculation of dissipation", nsd_);
    }
    break;
    case FLD::integrate_shape:
    {
      // integrate shape function for this element
      // (results assembled into element vector)
      return IntegrateShapeFunction(ele, discretization, lm, elevec1);
    }
    break;
    case FLD::calc_divop:
    {
      // calculate the integrated divergence operator
      return CalcDivOp(ele, discretization, lm, elevec1);
    }
    break;
    case FLD::calc_mat_deriv_u_and_rot_u:
    {
      // calculate material derivative at specified element coordinates
      return CalcMatDerivAndRotU(
          ele,
          params,
          discretization,
          lm,
          elevec1,
          elevec2,
          elevec3);
    }
    break;
    case FLD::void_fraction_gaussian_integration:
    {
      // calculate void fraction of the element for cavitation problems
      return ComputeVoidFraction(ele, params, discretization, lm, elevec1);
    }
    break;
    case FLD::calc_mass_matrix:
    {
      // compute element mass matrix
      return CalcMassMatrix(ele, discretization, lm, mat, elemat1);
    }
    break;
    case FLD::calc_turbulence_statistics:
    {
      if (nsd_ == 3)
      {
        if(ele->Owner() == discretization.Comm().MyPID())
        {
          //it is quite expensive to calculate second order derivatives, since a matrix has to be inverted...
          //so let's save time here
          is_higher_order_ele_=false;
          return CalcChannelStatistics(ele,params,discretization,lm,mat);
        }
      } // end if (nsd == 3)
      else dserror("action 'calc_turbulence_statistics' is a 3D specific action");
      return 0;
    }
    break;
    case FLD::calc_dt_via_cfl:
    {
      return CalcTimeStep(ele, discretization, lm, elevec1);
    }
    break;
    default:
      dserror("Unknown type of action '%i' for Fluid EvaluateService()", act);
    break;
  } // end of switch(act)

  return 0;
}


/*----------------------------------------------------------------------*
 * Action type: Integrate shape function
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::IntegrateShapeFunction(
    DRT::ELEMENTS::Fluid*     ele,
    DRT::Discretization&      discretization,
    const std::vector<int>&   lm,
    Epetra_SerialDenseVector& elevec1)
{
  // integrations points and weights
  return IntegrateShapeFunction( ele, discretization, lm, elevec1, intpoints_);
}


/*----------------------------------------------------------------------*
 * Action type: Integrate shape function
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::IntegrateShapeFunction(
    DRT::ELEMENTS::Fluid*     ele,
    DRT::Discretization&      discretization,
    const std::vector<int>&   lm            ,
    Epetra_SerialDenseVector& elevec1,
    const DRT::UTILS::GaussIntegration & intpoints)
{
  // --------------------------------------------------
  // construct views
  LINALG::Matrix<numdofpernode_*nen_,    1> vector(elevec1.A(),true);

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);
  // set element id
  eid_ = ele->Id();

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
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

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

  // set element id
  eid_ = ele->Id();

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
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

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


/*---------------------------------------------------------------------*
 | Action type: calc_mat_derivative_u                                  |
 | calculate material derivative of velocity               ghamm 01/13 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::CalcMatDerivAndRotU(
  DRT::ELEMENTS::Fluid*     ele,
  Teuchos::ParameterList&   params,
  DRT::Discretization&      discretization,
  std::vector<int>&         lm,
  Epetra_SerialDenseVector& elevec1,
  Epetra_SerialDenseVector& elevec2,
  Epetra_SerialDenseVector& elevec3)
{
  // fill the local element vector/matrix with the global values
  LINALG::Matrix<nsd_,nen_> eveln(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &eveln, NULL,"veln");
  LINALG::Matrix<nsd_,nen_> evelnm(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelnm, NULL,"velnm");

  // coordinates of the current integration point
  xsi_.Update(1.0, params.get<LINALG::Matrix<nsd_,1> >("elecoords"));

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

  // get velocities u_n and u_nm at integration point
  velint_.Multiply(eveln,funct_);
  LINALG::Matrix<nsd_,1> velintnm(true);
  velintnm.Multiply(evelnm,funct_);

  for (int isd=0;isd<nsd_;isd++)
  {
    elevec1[isd] = velint_(isd);
  }

  // get gradient of velocity at integration point
  vderxy_.MultiplyNT(eveln,deriv_);

  // calculate (u_n * nabla) u_n
  conv_old_.Multiply(vderxy_,velint_);

  double invdt = 1.0 / params.get<double>("timestep");

  // calculate (u_n - u_nm)/dt + (u_n * nabla) u_n
  conv_old_.Update(invdt, velint_, -invdt, velintnm, 1.0);

  for (int isd=0;isd<nsd_;isd++)
  {
    elevec2[isd] = conv_old_(isd);
  }

  // velocity gradient stored in vderxy_
  /*
     +-            -+
     | du   du   du |
     | --   --   -- |
     | dx   dy   dz |
     |              |
     | dv   dv   dv |
     | --   --   -- |
     | dx   dy   dz |
     |              |
     | dw   dw   dw |
     | --   --   -- |
     | dx   dy   dz |
     +-            -+
  */

  // rotation of fluid
  if(nsd_ == 3)
  {
    elevec3[0] = vderxy_(2,1) - vderxy_(1,2);
    elevec3[1] = vderxy_(0,2) - vderxy_(2,0);
    elevec3[2] = vderxy_(1,0) - vderxy_(0,1);
  }
  else
    dserror("not implemented");

  return 0;
}


/*---------------------------------------------------------------------*
 | Action type: void_fraction_gaussian_integration                     |
 | calculate void fraction for this element                ghamm 01/13 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::ComputeVoidFraction(
  DRT::ELEMENTS::Fluid*     ele,
  Teuchos::ParameterList&   params,
  DRT::Discretization&      discretization,
  std::vector<int>&         lm,
  Epetra_SerialDenseVector& elevec1)
{
  const double influence = params.get<double>("influence");
  LINALG::Matrix<3,1> pos = params.get<LINALG::Matrix<3,1> >("particle_pos");

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);
  // set element id
  eid_ = ele->Id();

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
  DRT::UTILS::GaussRule3D gaussrule = DRT::UTILS::intrule3D_undefined;
  if (nsd_ == 3)
  {
    const int gp_per_dir = params.get<int>("gp_per_dir");
    switch(gp_per_dir)
    {
    case 1:
      gaussrule = DRT::UTILS::intrule_hex_1point;
    break;
    case 2:
      gaussrule = DRT::UTILS::intrule_hex_8point;
    break;
    case 3:
      gaussrule = DRT::UTILS::intrule_hex_27point;
    break;
    case 4:
      gaussrule = DRT::UTILS::intrule_hex_64point;
    break;
    case 5:
      gaussrule = DRT::UTILS::intrule_hex_125point;
    break;
    case 6:
      gaussrule = DRT::UTILS::intrule_hex_216point;
    break;
    case 7:
      gaussrule = DRT::UTILS::intrule_hex_343point;
    break;
    case 8:
      gaussrule = DRT::UTILS::intrule_hex_512point;
    break;
    case 9:
      gaussrule = DRT::UTILS::intrule_hex_729point;
    break;
    }
  }
  else
    dserror("gauss rule not implemented");

  const DRT::UTILS::IntegrationPoints3D intpoints(gaussrule);

  //----------------------------------------------------------------------
  // loop over all gauss points of the actual element
  //----------------------------------------------------------------------
  for (int gp = 0; gp < intpoints.nquad; gp++)
  {
    // compute fac_
    const double e0 = intpoints.qxg[gp][0];
    const double e1 = intpoints.qxg[gp][1];
    const double e2 = intpoints.qxg[gp][2];

    // get shape functions and derivatives of the element
    DRT::UTILS::shape_function_3D(funct_,e0,e1,e2,distype);
    DRT::UTILS::shape_function_3D_deriv1(deriv_,e0,e1,e2,distype);


    xjm_.MultiplyNT(deriv_,xyze_);
    det_ = xji_.Invert(xjm_);

    if (det_ < 1E-16)
      dserror("GLOBAL ELEMENT NO.%i\nZERO OR NEGATIVE JACOBIAN DETERMINANT: %f", eid_, det_);

    // compute integration factor
    fac_ = intpoints.qwgt[gp]*det_;

    LINALG::Matrix<nsd_,1> x_gp(true);
    x_gp.Multiply(xyze_,funct_);


    // evaluate fourth order clipped polynomial
    double integrand = 1.0;
    for(int d=0; d<nsd_; ++d)
    {
      if( abs(x_gp(d) - pos(d)) < influence )
        integrand *= 15.0/16.0 * ( (pow((x_gp(d)-pos(d)), 4.0)/pow(influence, 5.0))
                                  - 2.0 * (pow((x_gp(d)-pos(d)), 2.0)/pow(influence, 3.0))
                                  + 1.0 / influence );
      else
      {
        integrand = 0.0;
        break;
      }
    }

    for (int i=0; i<nen_; ++i)
    {
      elevec1[0] += integrand*fac_*funct_(i);
    }

  }

  return 0;
}


/*----------------------------------------------------------------------*
 * Action type: Compute Div u                                 ehrl 12/12|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::ComputeDivU(
    DRT::ELEMENTS::Fluid*           ele,
    DRT::Discretization&            discretization,
    std::vector<int>&               lm,
    Epetra_SerialDenseVector&       elevec1)
{
  double area = 0.0;
  double divu = 0.0;

  //----------------------------------------------------------------------------
  //   Extract velocity/pressure from global vectors
  //----------------------------------------------------------------------------

  // fill the local element vector/matrix with the global values
  // af_genalpha: velocity/pressure at time n+alpha_F
  // np_genalpha: velocity at time n+alpha_F, pressure at time n+1, velocity for continuity equ. at time n+1
  // ost:         velocity/pressure at time n+1
  LINALG::Matrix<nsd_,nen_> evelaf(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelaf, NULL,"velaf");

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // set element id
  eid_ = ele->Id();

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

  {
    // evaluate shape functions and derivatives at element center
    EvalShapeFuncAndDerivsAtEleCenter();

//------------------------------------------------------------------
//                       INTEGRATION LOOP
//------------------------------------------------------------------

  // option to evaluate div u at Gauss point

  // loop over Gauss points if div u needs to be evaluated at the Gauss points
  /*
  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());
  */

    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    vderxy_.MultiplyNT(evelaf,derxy_);

    vdiv_= 0.0;
    for (int idim = 0; idim <nsd_; ++idim)
    {
      vdiv_ += vderxy_(idim, idim);
    }

    divu += vdiv_*fac_;
    area += fac_;
  }
  elevec1[0] = divu/area;
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
  const DRT::UTILS::GaussIntegration intpoints(distype, ele->Degree()*2+3);
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
  LINALG::Matrix<nsd_,nsd_> dervel(true);

  // error
  LINALG::Matrix<nsd_,1> deltavel(true);
  double         deltap=0.0;
  LINALG::Matrix<nsd_,nsd_> deltadervel(true);
  LINALG::Matrix<nsd_,nsd_> dervelint(true);

  const INPAR::FLUID::CalcError calcerr = DRT::INPUT::get<INPAR::FLUID::CalcError>(params,"calculate error");

  //----------------------------------------------------------------------------
  //   Extract velocity/pressure from global vectors
  //----------------------------------------------------------------------------

  // fill the local element vector/matrix with the global values
  // af_genalpha: velocity/pressure at time n+alpha_F
  // np_genalpha: velocity at time n+alpha_F, pressure at time n+1
  // ost:         velocity/pressure at time n+1
  LINALG::Matrix<nsd_,nen_> evelaf(true);
  LINALG::Matrix<nen_,1>    epreaf(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelaf, &epreaf,"velaf");

  // np_genalpha: additional vector for velocity at time n+1
  LINALG::Matrix<nsd_,nen_> evelnp(true);
  LINALG::Matrix<nen_,1>    eprenp(true);
  if (fldparatimint_->IsGenalphaNP())
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelnp, &eprenp,"velnp");

  //----------------------------------------------------------------------------
  //                         ELEMENT GEOMETRY
  //----------------------------------------------------------------------------

  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);
  // set element id
  eid_ = ele->Id();

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
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    // get velocity at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    velint_.Multiply(evelaf,funct_);

    // get pressure at integration point
    // (value at n+alpha_F for generalized-alpha scheme,
    //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
    double preint(true);
    if(fldparatimint_->IsGenalphaNP())
      preint= funct_.Dot(eprenp);
    else
      preint = funct_.Dot(epreaf);

    // H1 -error norm
    // compute first derivative of the velocity
    dervelint.MultiplyNT(evelaf,derxy_);

    // get coordinates at integration point
    LINALG::Matrix<nsd_,1> xyzint(true);
    xyzint.Multiply(xyze_,funct_);

    //  the error is evaluated at the specific time of the used time integration scheme
    //  n+alpha_F for generalized-alpha scheme
    //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)

    EvaluateAnalyticSolutionPoint(xyzint, fldparatimint_->Time(), calcerr, mat, u, p, dervel);

    if (calcerr == INPAR::FLUID::topoptchannel &&
        !(xyzint(1)>-0.2-1.0e-014 && xyzint(1)<0.2+1.0e-014))
      preint = 0.0;

    // compute difference between analytical solution and numerical solution
    deltap    = preint - p;
    deltavel.Update(1.0, velint_, -1.0, u);

    // H1 -error norm
    // compute error for first velocity derivative
    for(int i=0;i<nsd_;++i)
      for(int j=0;j<nsd_;++j)
        deltadervel(i,j)= dervelint(i,j) - dervel(i,j);

    // 0: delta velocity L2-error norm
    // 1: delta p L2-error norm
    // 2: delta velocity H1-error norm
    // 3: analytical velocity L2 norm
    // 4: analytical p L2 norm
    // 5: analytical velocity H1 norm

    // the error for the L2 and H1 norms are evaluated at the Gauss point
    for (int isd=0;isd<nsd_;isd++)
    {
      // integrate delta velocity for L2-error norm
      elevec1[0] += deltavel(isd)*deltavel(isd)*fac_;
      // integrate delta velocity for H1-error norm
      elevec1[2] += deltavel(isd)*deltavel(isd)*fac_;
      // integrate analytical velocity for L2 norm
      elevec1[3] += u(isd)*u(isd)*fac_;
      // integrate analytical velocity for H1 norm
      elevec1[5] += u(isd)*u(isd)*fac_;
    }
    // integrate delta p for L2-error norm
    elevec1[1] += deltap*deltap*fac_;
    //integrate analytical p for L2 norm
    elevec1[4] += p*p*fac_;

    // integrate delta velocity derivative for H1-error norm
    elevec1[2] += deltadervel.Dot(deltadervel)*fac_;
    // integrate analytical velocity for H1 norm
    elevec1[5] += dervel.Dot(dervel)*fac_;
  }

  return 0;
}


template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalc<distype>::EvaluateAnalyticSolutionPoint (
      const LINALG::Matrix<nsd_,1>      &xyzint,
      const double                       t,
      const INPAR::FLUID::CalcError      calcerr,
      const Teuchos::RCP<MAT::Material> &mat,
      LINALG::Matrix<nsd_,1>            &u,
      double                            &p,
      LINALG::Matrix<nsd_,nsd_>         &dervel
      )
{
  // Compute analytical solution
  switch(calcerr)
  {
  case INPAR::FLUID::beltrami_flow:
  {
    if (nsd_ == 3)
    {
      double visc = 1.;
      // get viscosity
      if (mat->MaterialType() == INPAR::MAT::m_fluid)
      {
        const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());

        // get constant kinematic viscosity
        visc = actmat->Viscosity()/actmat->Density();
      }
      else dserror("Material is not Newtonian Fluid");

      const double a      = M_PI/4.0;
      const double d      = M_PI/2.0;

      // compute analytical pressure
      p = -a*a/2.0 *
          ( std::exp(2.0*a*xyzint(0))
              + std::exp(2.0*a*xyzint(1))
              + std::exp(2.0*a*xyzint(2))
              + 2.0 * std::sin(a*xyzint(0) + d*xyzint(1)) * std::cos(a*xyzint(2) + d*xyzint(0)) * std::exp(a*(xyzint(1)+xyzint(2)))
              + 2.0 * std::sin(a*xyzint(1) + d*xyzint(2)) * std::cos(a*xyzint(0) + d*xyzint(1)) * std::exp(a*(xyzint(2)+xyzint(0)))
              + 2.0 * std::sin(a*xyzint(2) + d*xyzint(0)) * std::cos(a*xyzint(1) + d*xyzint(2)) * std::exp(a*(xyzint(0)+xyzint(1)))
          )* std::exp(-2.0*visc*d*d*t);

      // H1 -error norm
      // sacado data type replaces "double"
      typedef Sacado::Fad::DFad<double> FAD;  // for first derivs

      FAD x = xyzint(0);
      x.diff(0,3);  // independent variable 0 out of a total of 3

      FAD y = xyzint(1);
      y.diff(1,3);  // independent variable 1 out of a total of 3

      FAD z = xyzint(2);
      z.diff(2,3);  // independent variable 2 out of a total of 3

      // compute the function itself AND its derivatives w.r.t. ALL indep. variables
      FAD uu = -a * ( std::exp(a*x) * std::sin(a*y + d*z) +
          std::exp(a*z) * std::cos(a*x + d*y) ) * std::exp(-visc*d*d*t);
      FAD vv = -a * ( std::exp(a*y) * std::sin(a*z + d*x) +
          std::exp(a*x) * std::cos(a*y + d*z) ) * std::exp(-visc*d*d*t);
      FAD ww = -a * ( std::exp(a*z) * std::sin(a*x + d*y) +
          std::exp(a*y) * std::cos(a*z + d*x) ) * std::exp(-visc*d*d*t);

      u(0) = uu.val();
      u(1) = vv.val();
      u(2) = ww.val();

      dervel(0,0)=uu.dx(0);
      dervel(0,1)=uu.dx(1);
      dervel(0,2)=uu.dx(2);
      dervel(1,0)=vv.dx(0);
      dervel(1,1)=vv.dx(1);
      dervel(1,2)=vv.dx(2);
      dervel(2,0)=ww.dx(0);
      dervel(2,1)=ww.dx(1);
      dervel(2,2)=ww.dx(2);
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
    const double maxvel=1.0;
    const double height = 1.0;
    const double visc = 1.0;
    const double pressure_gradient = 10.0;

    // u_max = 1.25
    // y=0 is located in the middle of the channel
    if (nsd_ == 2)
    {
      p = 1.0;
      //p = -10*xyzint(0)+20;
      u(0) = maxvel -((height*height)/(2.0*visc)*pressure_gradient*(xyzint(1)/height)*(xyzint(1)/height));
      u(1) = 0.0;
    }
    else
      dserror("3D analytical solution is not implemented yet");
  }
  break;
  case INPAR::FLUID::topoptchannel:
  {
    const double visc = static_cast<const MAT::NewtonianFluid*>(mat.get())->Viscosity();
    // Y=xyzint(1); y=0 is located in the middle of the channel

    if (xyzint(1)>-0.2-1.0e-014 && xyzint(1)<0.2+1.0e-014)
    {
      u(0) = 1-25*xyzint(1)*xyzint(1);
      u(1) = 0.0;
      p = (xyzint(0)-0.5)*(10 - 50*visc); // 10 bof-o, 50 visc-part

      dervel(0,0)=0.0;
      dervel(0,1)=-50*xyzint(1);
      dervel(1,0)=0.0;
      dervel(1,1)=0.0;
    }
    else
    {
      u(0) = 0.0;
      u(1) = 0.0;
      //p = preint; //pressure error outside of channel not factored in
      p=0.0;

      dervel(0,0)=0.0;
      dervel(0,1)=0.0;
      dervel(1,0)=0.0;
      dervel(1,1)=0.0;
    }
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
  case INPAR::FLUID::fsi_fluid_pusher:
  {
    // get pointer to material in order to access density
    if (mat->MaterialType() == INPAR::MAT::m_fluid)
    {
      const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());

      // compute analytical velocities and pressure for
      // Note: consider offset of coordinate system
      // d(t) = -t^2
      {
        u(0) = -2.0 * t;
        u(1) = 0.0;
        u(2) = 0.0;
        p = actmat->Density() * 2.0 * (xyzint(0) + 1.5);
      }

//        // d(t) = -t^3
//        {
//          u(0) = -3.0 * t * t;
//          u(1) = 0.0;
//          u(2) = 0.0;
//          p = actmat->Density() * 6.0 * t * (xyzint(0) + 1.5);
//        }

//        // d(t) = -t^4
//        {
//          u(0) = -4.0 * t * t * t;
//          u(1) = 0.0;
//          u(2) = 0.0;
//          p = actmat->Density() * 12.0 * t * t * (xyzint(0) + 1.5);
//        }

//        // d(t) = -t^5
//        {
//          u(0) = -5.0 * t * t * t *t;
//          u(1) = 0.0;
//          u(2) = 0.0;
//          p = actmat->Density() * 20.0 * t * t * t * (xyzint(0) + 1.5);
//        }
    }
    else dserror("Material is not a Newtonian Fluid");
  }
  break;
  default:
    dserror("analytical solution is not defined");
    break;
  }

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
  Teuchos::RCP<MAT::Material> mat)
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
  BodyForce(ele,ebofoaf,eprescpgaf,escabofoaf);

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
  LINALG::Matrix<nen_,1>    eprenp(true);
  if (fldparatimint_->IsGenalphaNP())
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelnp, &eprenp,"velnp");

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

  if (fldparatimint_->IsGenalpha()) eveln.Clear();
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

  // set element id
  eid_ = ele->Id();

  // for scale similarity model:
  // get filtered veolcities and reynoldsstresses
  LINALG::Matrix<nsd_,nen_> evel_hat(true);
  LINALG::Matrix<nsd_*nsd_,nen_> ereynoldsstress_hat(true);
  if (fldpara_->TurbModAction() == INPAR::FLUID::scale_similarity
      or fldpara_->TurbModAction() == INPAR::FLUID::scale_similarity_basic)
  {
    Teuchos::RCP<Epetra_MultiVector> filtered_vel = params.get<Teuchos::RCP<Epetra_MultiVector> >("Filtered velocity");
    Teuchos::RCP<Epetra_MultiVector> fs_vel = params.get<Teuchos::RCP<Epetra_MultiVector> >("Fine scale velocity");
    Teuchos::RCP<Epetra_MultiVector> filtered_reystre = params.get<Teuchos::RCP<Epetra_MultiVector> >("Filtered reynoldsstress");

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
    Teuchos::RCP<Epetra_Vector> ele_CsDeltaSq = params.sublist("TURBULENCE MODEL").get<Teuchos::RCP<Epetra_Vector> >("col_Cs_delta_sq");
    Teuchos::RCP<Epetra_Vector> ele_CiDeltaSq = params.sublist("TURBULENCE MODEL").get<Teuchos::RCP<Epetra_Vector> >("col_Ci_delta_sq");
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
  Teuchos::RCP<std::vector<double> > planecoords  = params.get<Teuchos::RCP<std::vector<double> > >("planecoords_",Teuchos::null);
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
  double eps_graddiv       = 0.0;
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
  EvalShapeFuncAndDerivsAtEleCenter();

  // set element area or volume
  vol = fac_;

  //------------------------------------------------------------------------
  // potential evaluation of material parameters, subgrid viscosity
  // and/or stabilization parameters at element center
  //------------------------------------------------------------------------
  // get material parameters at element center
  if (not fldpara_->MatGp() or not fldpara_->TauGp())
  {

    GetMaterialParams(mat,evelaf,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam,vol);

    // calculate all-scale or fine-scale subgrid viscosity at element center
    visceff_ = visc_;
    if (fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky or fldpara_->TurbModAction() == INPAR::FLUID::vreman)
    {
      CalcSubgrVisc(evelaf,vol,Cs_delta_sq,Ci_delta_sq);
      // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
      visceff_ += sgvisc_;
    }
    else if (fldpara_->Fssgv() != INPAR::FLUID::no_fssgv)
      CalcFineScaleSubgrVisc(evelaf,fsevelaf,vol);
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
        //GetMaterialParams(material,evelaf,escaaf,escaam,thermpressaf,thermpressam,thermpressdtam,vol);
        GetMaterialParams(mat,evelaf,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam,vol);

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
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

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
    // (value at n+alpha_F for generalized-alpha scheme,
    //  value at n+alpha_F for generalized-alpha-NP schemen, n+1 otherwise)
    if(fldparatimint_->IsGenalphaNP())
      gradp_.Multiply(derxy_,eprenp);
    else
      gradp_.Multiply(derxy_,epreaf);

    // get bodyforce at integration point
    // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
    bodyforce_.Multiply(ebofoaf,funct_);
    // get prescribed pressure gradient acting as body force
    // (required for turbulent channel flow)
    generalbodyforce_.Multiply(eprescpgaf,funct_);

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
    if (not fldparatimint_->IsGenalphaNP())
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
      GetMaterialParams(mat,evelaf,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam,vol);

      // calculate all-scale or fine-scale subgrid viscosity at integration point
      visceff_ = visc_;
      if (fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky or fldpara_->TurbModAction() == INPAR::FLUID::vreman)
      {
        CalcSubgrVisc(evelaf,vol,Cs_delta_sq,Ci_delta_sq);
        // effective viscosity = physical viscosity + (all-scale) subgrid viscosity
        visceff_ += sgvisc_;
      }
      else if (fldpara_->Fssgv() != INPAR::FLUID::no_fssgv)
        CalcFineScaleSubgrVisc(evelaf,fsevelaf,vol);
    }

    // potential evaluation of coefficient of multifractal subgrid-scales at integration point
    if (fldpara_->TurbModAction() == INPAR::FLUID::multifractal_subgrid_scales)
    {
      if (fldpara_->BGp())
      {
        // make sure to get material parameters at gauss point
        if (not fldpara_->MatGp())
          GetMaterialParams(mat,evelaf,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam,vol);

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
      if (not fldparatimint_->IsGenalpha()) dserror("Does ComputeGalRHSContEq() for ost really the right thing?");
      // remark: I think the term theta*u^n+1*nabla T^n+1 is missing.
      //         Moreover, the resulting conres_old_ should be multiplied by theta (see monres_old_).

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
        if (fldpara_->TurbModAction() == INPAR::FLUID::smagorinsky or fldpara_->TurbModAction() == INPAR::FLUID::dynamic_smagorinsky or fldpara_->TurbModAction() == INPAR::FLUID::vreman)
        visceff_ += sgvisc_;
      }
    }

    // evaluate momentum residual once for all stabilization right hand sides
    if (fldparatimint_->IsGenalpha())
    {
      // get acceleration at time n+alpha_M at integration point
      accint_.Multiply(eaccam,funct_);

      for (int rr=0;rr<nsd_;++rr)
      {
        momres_old_(rr) = densam_*accint_(rr)+densaf_*conv_old_(rr)+gradp_(rr)
                       -2*visceff_*visc_old_(rr)-densaf_*bodyforce_(rr)-generalbodyforce_(rr);
      }
    }
    else
    {
      rhsmom_.Update((densn_/fldparatimint_->Dt()),histmom_,densaf_*fldparatimint_->Theta(),bodyforce_);
      // and pressure gradient prescribed as body force
      // caution: not density weighted
      rhsmom_.Update(fldparatimint_->Theta(),generalbodyforce_,1.0);
      // compute instationary momentum residual:
      // momres_old = u_(n+1)/dt + theta ( ... ) - histmom_/dt - theta*bodyforce_
      for (int rr=0;rr<nsd_;++rr)
      {
        momres_old_(rr) = (densaf_*velint_(rr)/fldparatimint_->Dt()
                         +fldparatimint_->Theta()*(densaf_*conv_old_(rr)+gradp_(rr)
                         -2*visceff_*visc_old_(rr)))-rhsmom_(rr);
      }
    }


    //---------------------------------------------------------------
    // element average dissipation and production rates
    //---------------------------------------------------------------

    //---------------------------------------------------------------
    // residual-based subgrid-scale modeling terms
    //---------------------------------------------------------------

    // dissipation by supg-stabilization
    if (fldpara_->SUPG())
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
    if (fldpara_->PSPG())
    {
      for (int rr=0;rr<nsd_;rr++)
      {
        eps_pspg += fac_ * gradp_(rr,0) * tau_(1) * momres_old_(rr,0);
      }
    }

    // dissipation by continuity-stabilization
    if (fldpara_->CStab())
    {
      eps_graddiv += fac_ * vdiv_ * tau_(2) * conres_old_;
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
      if (fldpara_->PhysicalType() == INPAR::FLUID::loma)
      {
        dserror("Read warning before usage!");
        // Warning: Here, we should use the deviatoric part of the strain-rate tensor.
        //          However, I think this is not done in the element Sysmat-routine.
        //          Hence, I skipped it here.
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
  eps_graddiv /= vol;
  eps_pspg /= vol;

  Teuchos::RCP<std::vector<double> > incrvol           = params.get<Teuchos::RCP<std::vector<double> > >("incrvol"          );

  Teuchos::RCP<std::vector<double> > incr_eps_visc      = params.get<Teuchos::RCP<std::vector<double> > >("incr_eps_visc"    );
  Teuchos::RCP<std::vector<double> > incr_eps_conv      = params.get<Teuchos::RCP<std::vector<double> > >("incr_eps_conv"    );
  Teuchos::RCP<std::vector<double> > incr_eps_smag      = params.get<Teuchos::RCP<std::vector<double> > >("incr_eps_eddyvisc");
  Teuchos::RCP<std::vector<double> > incr_eps_avm3      = params.get<Teuchos::RCP<std::vector<double> > >("incr_eps_avm3"    );
  Teuchos::RCP<std::vector<double> > incr_eps_mfs       = params.get<Teuchos::RCP<std::vector<double> > >("incr_eps_mfs"     );
  Teuchos::RCP<std::vector<double> > incr_eps_mfscross  = params.get<Teuchos::RCP<std::vector<double> > >("incr_eps_mfscross");
  Teuchos::RCP<std::vector<double> > incr_eps_mfsrey    = params.get<Teuchos::RCP<std::vector<double> > >("incr_eps_mfsrey"  );
  Teuchos::RCP<std::vector<double> > incr_eps_supg      = params.get<Teuchos::RCP<std::vector<double> > >("incr_eps_supg"    );
  Teuchos::RCP<std::vector<double> > incr_eps_cross     = params.get<Teuchos::RCP<std::vector<double> > >("incr_eps_cross"   );
  Teuchos::RCP<std::vector<double> > incr_eps_rey       = params.get<Teuchos::RCP<std::vector<double> > >("incr_eps_rey"     );
  Teuchos::RCP<std::vector<double> > incr_eps_graddiv     = params.get<Teuchos::RCP<std::vector<double> > >("incr_eps_graddiv"   );
  Teuchos::RCP<std::vector<double> > incr_eps_pspg      = params.get<Teuchos::RCP<std::vector<double> > >("incr_eps_pspg"    );

  Teuchos::RCP<std::vector<double> > incrhk            = params.get<Teuchos::RCP<std::vector<double> > >("incrhk"           );
  Teuchos::RCP<std::vector<double> > incrhbazilevs     = params.get<Teuchos::RCP<std::vector<double> > >("incrhbazilevs"    );
  Teuchos::RCP<std::vector<double> > incrstrle         = params.get<Teuchos::RCP<std::vector<double> > >("incrstrle"        );
  Teuchos::RCP<std::vector<double> > incrgradle        = params.get<Teuchos::RCP<std::vector<double> > >("incrgradle"       );

  Teuchos::RCP<std::vector<double> > incrres           = params.get<Teuchos::RCP<std::vector<double> > >("incrres"          );
  Teuchos::RCP<std::vector<double> > incrres_sq        = params.get<Teuchos::RCP<std::vector<double> > >("incrres_sq"       );
  Teuchos::RCP<std::vector<double> > incrabsres        = params.get<Teuchos::RCP<std::vector<double> > >("incrabsres"       );
  Teuchos::RCP<std::vector<double> > incrtauinvsvel    = params.get<Teuchos::RCP<std::vector<double> > >("incrtauinvsvel"   );

  Teuchos::RCP<std::vector<double> > incrsvelaf        = params.get<Teuchos::RCP<std::vector<double> > >("incrsvelaf"       );
  Teuchos::RCP<std::vector<double> > incrsvelaf_sq     = params.get<Teuchos::RCP<std::vector<double> > >("incrsvelaf_sq"    );
  Teuchos::RCP<std::vector<double> > incrabssvelaf     = params.get<Teuchos::RCP<std::vector<double> > >("incrabssvelaf"    );

  Teuchos::RCP<std::vector<double> > incrresC          = params.get<Teuchos::RCP<std::vector<double> > >("incrresC"         );
  Teuchos::RCP<std::vector<double> > incrresC_sq       = params.get<Teuchos::RCP<std::vector<double> > >("incrresC_sq"      );
  Teuchos::RCP<std::vector<double> > spressnp          = params.get<Teuchos::RCP<std::vector<double> > >("incrspressnp"     );
  Teuchos::RCP<std::vector<double> > spressnp_sq       = params.get<Teuchos::RCP<std::vector<double> > >("incrspressnp_sq"  );

  Teuchos::RCP<std::vector<double> > incrtauC          = params.get<Teuchos::RCP<std::vector<double> > >("incrtauC"         );
  Teuchos::RCP<std::vector<double> > incrtauM          = params.get<Teuchos::RCP<std::vector<double> > >("incrtauM"         );

  Teuchos::RCP<std::vector<double> > incrcrossstress   = params.get<Teuchos::RCP<std::vector<double> > >("incrcrossstress"  );
  Teuchos::RCP<std::vector<double> > incrreystress     = params.get<Teuchos::RCP<std::vector<double> > >("incrreystress"    );

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
  (*incr_eps_graddiv   )[nlayer] += eps_graddiv      ;
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
      \param graddiv            (i) boolean flag for stabilisation (pass-through)
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

  if(fldparatimint_->IsGenalphaNP())
    dserror("FD check not available for NP genalpha!!");

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
  printf("| FINITE DIFFERENCE CHECK FOR ELEMENT %5d |\n",eid_);
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

        if (fldparatimint_->IsGenalpha())
        {
          checkepreaf(nn)+=fldparatimint_->AlphaF()*epsilon;
        }
        else
        {
          checkepreaf(nn)+=epsilon;
        }
      }
      else
      {
        printf("velocity dof %d (%d)\n",rr,nn);

        if (fldparatimint_->IsGenalpha())
        {
          checkevelaf(rr,nn)+=fldparatimint_->AlphaF()*epsilon;
          checkeaccam(rr,nn)+=fldparatimint_->AlphaM()/(fldparatimint_->Gamma()*fldparatimint_->Dt())*epsilon;
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
        if (fldparatimint_->IsGenalpha())
        {
          val   =-(eforce(mm)   /(epsilon))*(fldparatimint_->Gamma()*fldparatimint_->Dt())/(fldparatimint_->AlphaM());
          lin   =-(eforce(mm)   /(epsilon))*(fldparatimint_->Gamma()*fldparatimint_->Dt())/(fldparatimint_->AlphaM())+estif(mm,dof);
          nonlin=-(checkvec1(mm)/(epsilon))*(fldparatimint_->Gamma()*fldparatimint_->Dt())/(fldparatimint_->AlphaM());
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


/*-------------------------------------------------------------------------------*
 |find elements of inflow section                                rasthofer 10/12 |
 *-------------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidEleCalc<distype>::InflowElement(DRT::Element* ele)
{
   is_inflow_ele_ = false;

  std::vector<DRT::Condition*> myinflowcond;

  // check whether all nodes have a unique inflow condition
  DRT::UTILS::FindElementConditions(ele, "TurbulentInflowSection", myinflowcond);
  if (myinflowcond.size()>1)
    dserror("More than one inflow condition on one node!");

  if (myinflowcond.size()==1)
    is_inflow_ele_ = true;

  return;
}

/*-----------------------------------------------------------------------------*
 | Calculate element mass matrix                               mayr.mt 05/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::CalcMassMatrix(
    DRT::ELEMENTS::Fluid*                ele,
//    Teuchos::ParameterList&              params,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm,
    Teuchos::RCP<MAT::Material> &        mat,
    Epetra_SerialDenseMatrix&            elemat1_epetra
    )
{
  // set element id
  eid_ = ele->Id();

  // ---------------------------------------------------------------------------
  // Prepare material parameters
  // ---------------------------------------------------------------------------
  // Since we need only the density, we use a lot of dummy values.

  // create dummy matrices
  LINALG::Matrix<nsd_, nen_> mat1(true);
  LINALG::Matrix<nen_, 1> mat2(true);

  GetMaterialParams(mat, mat1, mat2, mat2, mat2, 0.0, 0.0, 0.0, 0.0, 0.0);

  // ---------------------------------------------------------------------------
  // Geometry
  // ---------------------------------------------------------------------------
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // Do ALE specific updates if necessary
  if (ele->IsAle())
  {
    LINALG::Matrix<nsd_,nen_> edispnp(true);
    ExtractValuesFromGlobalVector(discretization, lm, *rotsymmpbc_, &edispnp, NULL, "dispnp");

    // get new node positions of ALE mesh
     xyze_ += edispnp;
  }

  // definition of matrices
  LINALG::Matrix<nen_*nsd_,nen_*nsd_> estif_u(true);

  // ---------------------------------------------------------------------------
  // Integration loop
  // ---------------------------------------------------------------------------
  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    for (int ui=0; ui<nen_; ++ui)
    {
      for (int vi=0; vi<nen_; ++vi)
      {
        for (int jdim= 0; jdim<nsd_;++jdim)
        {
          for (int idim = 0; idim<nsd_; ++idim)
          {
            estif_u(nsd_*vi+idim,nsd_*ui+jdim) += funct_(vi)*funct_(ui)*fac_*densaf_;
          } // end for (idim)
        } // end for (jdim)
      } // end for (vi)
    } // end for (ui)
  } // end of integration loop

  // ---------------------------------------------------------------------------
  // Add velocity-velocity part to matrix
  // ---------------------------------------------------------------------------
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
          elemat1_epetra(numdof_vi+idim, numdof_ui_jdim) += estif_u(nsd_vi+idim, nsd_ui_jdim);
        } // end for (idim)
      } // end for (vi)
    } // end for (jdim)
  } // end for (ui)

  return 0;
}

/*-----------------------------------------------------------------------------*
 | Calculate channel statistics                                     bk 05/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::CalcChannelStatistics(
    DRT::ELEMENTS::Fluid*                ele,
    Teuchos::ParameterList&              params,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm,
    Teuchos::RCP<MAT::Material> &        mat
    )
{
  //moved here from fluid_ele_evaluate_utils,  f3_calc_means()

  /*!
    \brief calculate spatial mean values for channel flow
    (requires wall parallel layers of elements)

                                                             gammi 07/07

    <pre>

    this method assumes that each 2 dimensional integration element
    in the homogeneous plane is parallel to the wall!!!

    The necessary element integration is done in here. The element
    is cut into two (HEX8) or three (quadratic elements) planes (plus
    additional planes for visualisation purposes, defined by planes
    vector), the spatial functions (velocity, pressure etc.) are
    integrated over this plane and this element contribution is added
    to a processor local vector (see formulas below for a exact
    description of the output).
    It is assumed that the sampling planes are distributed equidistant
    in the element. The result is normalized by the area afterwards


                       ^ normdirect       integration plane
                       |                /
                       |               /
                       |
                 +-----|-------------+
                /|     |            /|
               / |     |           / |
              /  |     |          /  |
             /   |     |         /   |
            /    +-----|--------/----+ ---- additional integration
           /    /|     |       /    /|      plane (for quadratic elements)
          /    / |     |      /    / |
         +-------------------+    /  |
         |   /   |     *-----|---+------------>
         |  /    +----/------|--/----+         inplanedirect[1]
         | /    /    /       | /    /
         |/    /    /        |/    /   \
         +---------+---------+    /     \
         |   /    /          |   /       integration plane
         |  /    /           |  /
         | /    /            | /
         |/    /             |/
         +----/--------------+
             /
            /   inplanedirect[0]


    Example for a mean value evaluation:

           1.0       /                     1.0      /
    _               |                              |            detJ
    u = -------- *  | u(x,y,z) dx dz =  -------- * | u(r,s,t) * ---- dr ds
        +---        |                   +---       |             h
         \         / A                   \        /  [-1:1]^2     y
         / area                          / area
        +---                            +---                    +--+
                                                               Jacobi-
                                                             determinant
                                                              computed
                                                           from 3d mapping
                                                         h  is the (constant)
                                                          y
                                                        height of the element

    The method computes:
                        _             _             _             _
               numele * u  , numele * v  , numele * w  , numele * p
                        ___           ___           ___           ___
                         ^2            ^2            ^2            ^2
    and        numele * u  , numele * v  , numele * w  , numele * p

                        _ _           _ _           _ _
    as well as  numele * u*v, numele * u*w, numele * v*w


    as well as numele and the element area.
    All results are communicated via the parameter list!

    </pre>

   */

  // --------------------------------------------------
  // extract velocity and pressure from global
  // distributed vectors
  // --------------------------------------------------
  // velocity and pressure values (n+1)
  Teuchos::RCP<const Epetra_Vector> velnp
  = discretization.GetState("u and p (n+1,converged)");
  if (velnp==Teuchos::null) dserror("Cannot get state vector 'velnp'");

  // extract local values from the global vectors
  std::vector<double> mysol  (lm.size());
  DRT::UTILS::ExtractMyValues(*velnp,mysol,lm);
  // get view of solution and subgrid-viscosity vector
  LINALG::Matrix<4*nen_,1> sol(&(mysol[0]),true);

  // the plane normal tells you in which plane the integration takes place
  const int normdirect = params.get<int>("normal direction to homogeneous plane");

  // the vector planes contains the coordinates of the homogeneous planes (in
  // wall normal direction)
  Teuchos::RCP<std::vector<double> > planes = params.get<Teuchos::RCP<std::vector<double> > >("coordinate vector for hom. planes");

  // get the pointers to the solution vectors
  Teuchos::RCP<std::vector<double> > sumarea= params.get<Teuchos::RCP<std::vector<double> > >("element layer area");

  Teuchos::RCP<std::vector<double> > sumu   = params.get<Teuchos::RCP<std::vector<double> > >("mean velocity u");
  Teuchos::RCP<std::vector<double> > sumv   = params.get<Teuchos::RCP<std::vector<double> > >("mean velocity v");
  Teuchos::RCP<std::vector<double> > sumw   = params.get<Teuchos::RCP<std::vector<double> > >("mean velocity w");
  Teuchos::RCP<std::vector<double> > sump   = params.get<Teuchos::RCP<std::vector<double> > >("mean pressure p");

  Teuchos::RCP<std::vector<double> > sumsqu = params.get<Teuchos::RCP<std::vector<double> > >("mean value u^2");
  Teuchos::RCP<std::vector<double> > sumsqv = params.get<Teuchos::RCP<std::vector<double> > >("mean value v^2");
  Teuchos::RCP<std::vector<double> > sumsqw = params.get<Teuchos::RCP<std::vector<double> > >("mean value w^2");
  Teuchos::RCP<std::vector<double> > sumuv  = params.get<Teuchos::RCP<std::vector<double> > >("mean value uv");
  Teuchos::RCP<std::vector<double> > sumuw  = params.get<Teuchos::RCP<std::vector<double> > >("mean value uw");
  Teuchos::RCP<std::vector<double> > sumvw  = params.get<Teuchos::RCP<std::vector<double> > >("mean value vw");
  Teuchos::RCP<std::vector<double> > sumsqp = params.get<Teuchos::RCP<std::vector<double> > >("mean value p^2");

  // get node coordinates of element
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // Do ALE specific updates if necessary
  if (ele->IsAle())
  {
    LINALG::Matrix<nsd_,nen_> edispnp(true);
    ExtractValuesFromGlobalVector(discretization, lm, *rotsymmpbc_, &edispnp, NULL, "dispnp");

    // get new node positions of ALE mesh
     xyze_ += edispnp;

     for (int inode=0; inode<nen_; inode++)
     {
       if(abs(edispnp(normdirect,inode))>1e-6)
       {
         dserror("no sampling possible if homogeneous planes are not conserved\n");
       }
     }

  }

  if(distype == DRT::Element::hex8
     ||
     distype == DRT::Element::hex27
     ||
     distype == DRT::Element::hex20)
  {
    double min = xyze_(normdirect,0);
    double max = xyze_(normdirect,0);

    // set maximum and minimum value in wall normal direction
    for(int inode=0;inode<nen_;inode++)
    {
      if(min > xyze_(normdirect,inode))
      {
        min=xyze_(normdirect,inode);
      }
      if(max < xyze_(normdirect,inode))
      {
         max=xyze_(normdirect,inode);
      }
    }

    // determine the ids of the homogeneous planes intersecting this element
    std::set<int> planesinele;
    for(unsigned nplane=0;nplane<planes->size();++nplane)
    {
    // get all available wall normal coordinates
      for(int nn=0;nn<nsd_;++nn)
      {
        if (min-2e-9 < (*planes)[nplane] && max+2e-9 > (*planes)[nplane])
        {
          planesinele.insert(nplane);
         }
      }
    }

    // remove lowest layer from planesinele to avoid double calculations. This is not done
    // for the first level (index 0) --- if deleted, shift the first integration point in
    // wall normal direction
    // the shift depends on the number of sampling planes in the element
    double shift=0;

    // set the number of planes which cut the element
    const int numplanesinele = planesinele.size();

    if(*planesinele.begin() != 0)
    {
      // this is not an element of the lowest element layer
      planesinele.erase(planesinele.begin());

      shift=2.0/(static_cast<double>(numplanesinele-1));
    }
    else
    {
      // this is an element of the lowest element layer. Increase the counter
      // in order to compute the total number of elements in one layer
      int* count = params.get<int*>("count processed elements");

      (*count)++;
    }

    // determine the orientation of the rst system compared to the xyz system
    int elenormdirect=-1;
    bool upsidedown =false;
    // the only thing of interest is how normdirect is oriented in the
    // element coordinate system
    if(xyze_(normdirect,4)-xyze_(normdirect,0)>2e-9)
    {
      // t aligned
      elenormdirect =2;
    }
    else if (xyze_(normdirect,3)-xyze_(normdirect,0)>2e-9)
    {
      // s aligned
      elenormdirect =1;
    }
    else if (xyze_(normdirect,1)-xyze_(normdirect,0)>2e-9)
    {
      // r aligned
      elenormdirect =0;
    }
    else if(xyze_(normdirect,4)-xyze_(normdirect,0)<-2e-9)
    {
      // -t aligned
      elenormdirect =2;
      upsidedown =true;
    }
    else if (xyze_(normdirect,3)-xyze_(normdirect,0)<-2e-9)
    {
      // -s aligned
      elenormdirect =1;
      upsidedown =true;
    }
    else if (xyze_(normdirect,1)-xyze_(normdirect,0)<-2e-9)
    {
      // -r aligned
      elenormdirect =0;
      upsidedown =true;
    }
    else
    {
      dserror("cannot determine orientation of plane normal in local coordinate system of element");
    }
    std::vector<int> inplanedirect;
    {
      std::set <int> inplanedirectset;
      for(int i=0;i<3;++i)
      {
         inplanedirectset.insert(i);
      }
      inplanedirectset.erase(elenormdirect);

      for(std::set<int>::iterator id = inplanedirectset.begin();id!=inplanedirectset.end() ;++id)
      {
        inplanedirect.push_back(*id);
      }
    }

    // get the quad9 gaussrule for the in plane integration
    DRT::UTILS::GaussIntegration intpoints( DRT::Element::quad9 );
//    const DRT::UTILS::IntegrationPoints2D  intpoints(DRT::UTILS::intrule_quad_9point);

    // a hex8 element has two levels, the hex20 and hex27 element have three layers to sample
    // (now we allow even more)
    double layershift=0;

    // loop all levels in element
    for(std::set<int>::const_iterator id = planesinele.begin();id!=planesinele.end() ;++id)
    {
      // reset temporary values
      double area=0;

      double ubar=0;
      double vbar=0;
      double wbar=0;
      double pbar=0;

      double usqbar=0;
      double vsqbar=0;
      double wsqbar=0;
      double uvbar =0;
      double uwbar =0;
      double vwbar =0;
      double psqbar=0;

      // get the integration point in wall normal direction
      double e[3];

      e[elenormdirect]=-1.0+shift+layershift;
      if(upsidedown)
      {
        e[elenormdirect]*=-1;
      }

      // start loop over integration points in layer
      for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
      {
        // get the other gauss point coordinates
        for(int i=0;i<2;++i)
        {
          e[inplanedirect[i]]=iquad.Point()[i];
        }

        const double* econst=e;

        // evaluate shape functions and derivatives at integration point
        EvalShapeFuncAndDerivsAtIntPoint(econst,iquad.Weight());

        // we assume that every plane parallel to the wall is preserved
        // hence we can compute the jacobian determinant of the 2d cutting
        // element by replacing max-min by one on the diagonal of the
        // jacobi matrix (the two non-diagonal elements are zero)
        if(xjm_(elenormdirect,normdirect)<0)
        {
          xjm_(elenormdirect,normdirect)=-1.0;
        }
        else
        {
          xjm_(elenormdirect,normdirect)= 1.0;
        }

        //get local copy of determinant/ adjusted to plane integration
        const double det =
          xjm_(0,0)*xjm_(1,1)*xjm_(2,2)
          +
          xjm_(0,1)*xjm_(1,2)*xjm_(2,0)
          +
          xjm_(0,2)*xjm_(1,0)*xjm_(2,1)
          -
          xjm_(0,2)*xjm_(1,1)*xjm_(2,0)
          -
          xjm_(0,0)*xjm_(1,2)*xjm_(2,1)
          -
          xjm_(0,1)*xjm_(1,0)*xjm_(2,2);

        // check for degenerated elements
        if (det <= 0.0)
        {
          dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
        }

#ifdef DEBUG
    // check whether this gausspoint is really inside the desired plane
    {
      double x[3];
      x[0]=0;
      x[1]=0;
      x[2]=0;
      for(int inode=0;inode<nen_;inode++)
      {
        x[0]+=funct_(inode)*xyze_(0,inode);
        x[1]+=funct_(inode)*xyze_(1,inode);
        x[2]+=funct_(inode)*xyze_(2,inode);
      }

      if(abs(x[normdirect]-(*planes)[*id])>2e-9)
      {
        dserror("Mixing up element cut planes during integration");
      }
    }
#endif

        //interpolated values at gausspoints
        double ugp=0;
        double vgp=0;
        double wgp=0;
        double pgp=0;

        // the computation of this jacobian determinant from the 3d
        // mapping is based on the assumption that we do not deform
        // our elements in wall normal direction!
        const double fac=det*iquad.Weight();

        // increase area of cutting plane in element
        area += fac;

        for(int inode=0;inode<nen_;inode++)
        {
          int finode=inode*4;

          ugp  += funct_(inode)*sol(finode++);
          vgp  += funct_(inode)*sol(finode++);
          wgp  += funct_(inode)*sol(finode++);
          pgp  += funct_(inode)*sol(finode  );
        }

        // add contribution to integral

        double dubar  = ugp*fac;
        double dvbar  = vgp*fac;
        double dwbar  = wgp*fac;
        double dpbar  = pgp*fac;

        ubar   += dubar;
        vbar   += dvbar;
        wbar   += dwbar;
        pbar   += dpbar;

        usqbar += ugp*dubar;
        vsqbar += vgp*dvbar;
        wsqbar += wgp*dwbar;
        uvbar  += ugp*dvbar;
        uwbar  += ugp*dwbar;
        vwbar  += vgp*dwbar;
        psqbar += pgp*dpbar;
      } // end loop integration points

      // add increments from this layer to processor local vectors
      (*sumarea)[*id] += area;

      (*sumu   )[*id] += ubar;
      (*sumv   )[*id] += vbar;
      (*sumw   )[*id] += wbar;
      (*sump   )[*id] += pbar;

      (*sumsqu )[*id] += usqbar;
      (*sumsqv )[*id] += vsqbar;
      (*sumsqw )[*id] += wsqbar;
      (*sumuv  )[*id] += uvbar;
      (*sumuw  )[*id] += uwbar;
      (*sumvw  )[*id] += vwbar;
      (*sumsqp )[*id] += psqbar;

      // jump to the next layer in the element.
      // in case of an hex8 element, the two coordinates are -1 and 1(+2)
      // for quadratic elements with three sample planes, we have -1,0(+1),1(+2)

      layershift+=2.0/(static_cast<double>(numplanesinele-1));
    }
  }
  else if(distype == DRT::Element::nurbs8 || distype == DRT::Element::nurbs27)
  {
    // get size of planecoords
    int size = planes->size();

    DRT::NURBS::NurbsDiscretization* nurbsdis
      =
      dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(discretization));

    if(nurbsdis == NULL)
    {
      dserror("we need a nurbs discretisation for nurbs elements\n");
    }

    // get nurbs dis' element numbers
    std::vector<int> nele_x_mele_x_lele(nurbsdis->Return_nele_x_mele_x_lele(0));

    // use size of planes and mele to determine number of layers
    int numsublayers=(size-1)/nele_x_mele_x_lele[1];

    // get the knotvector itself
    Teuchos::RCP<DRT::NURBS::Knotvector> knots=nurbsdis->GetKnotVector();

    DRT::Node**   nodes = ele->Nodes();

    // get gid, location in the patch
    int gid = ele->Id();

    std::vector<int> ele_cart_id(3);

    int npatch = -1;

    knots->ConvertEleGidToKnotIds(gid,npatch,ele_cart_id);
    if(npatch!=0)
    {
      dserror("expected single patch nurbs problem for calculating means");
    }

    bool zero_size = false;
    zero_size = knots->GetEleKnots(myknots_,gid);

    // if we have a zero sized element due to a interpolated
    // point --- exit here
    if(zero_size)
    {
      return 0;
    }

    for (int inode=0; inode<nen_; ++inode)
    {
      DRT::NURBS::ControlPoint* cp
    =
    dynamic_cast<DRT::NURBS::ControlPoint* > (nodes[inode]);

      weights_(inode) = cp->W();
    }

    // there's one additional plane for the last element layer
    int endlayer=0;
    if(ele_cart_id[1]!=nele_x_mele_x_lele[1]-1)
    {
      endlayer=numsublayers;
    }
    else
    {
      endlayer=numsublayers+1;
    }



    //!!see below for more information in green!!
    //this is dangerous because we don't check anywhere, if the wall normal points in y direction.
    //please make this more general!
    //we don't have a test case for this routine either
    dserror("Warning: Nurbs channel statistics work only if the element wall normal points in y direction.");


    // loop layers in element
    for(int rr=0;rr<endlayer;++rr)
    {
      // set gauss point coordinates
      double gp[3];
      gp[1]=-1.0+rr*2.0/((double)numsublayers);

      // get the quad9 gaussrule for the in plane integration
      DRT::UTILS::GaussIntegration intpoints( DRT::Element::quad9 );

      // reset temporary values
      double area=0;

      double ubar=0;
      double vbar=0;
      double wbar=0;
      double pbar=0;

      double usqbar=0;
      double vsqbar=0;
      double wsqbar=0;
      double uvbar =0;
      double uwbar =0;
      double vwbar =0;
      double psqbar=0;


      // start loop over integration points in layer
      for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints.begin(); iquad!=intpoints.end(); ++iquad )
      {

        // get the other gauss point coordinates
        //here we assume that the element wall normal points in y direction
        gp[0]=iquad.Point()[0];
        gp[2]=iquad.Point()[1];

        const double* gpconst=gp;
        EvalShapeFuncAndDerivsAtIntPoint(gpconst,iquad.Weight());

        // we assume that every plane parallel to the wall is preserved
        // hence we can compute the jacobian determinant of the 2d cutting
        // element by replacing max-min by one on the diagonal of the
        // jacobi matrix (the two non-diagonal elements are zero)

        //here we still have the bug with the element normal directions
        //but we have to find out the element wall normal first to correct it
        //this part of the code works only if all the normals point in y direction!
        //please change the following lines to xjm_(elenormdirect,normdirect)
        //but you have to get elenormdirect first...
        if(xjm_(normdirect,normdirect)<0)
        {
          xjm_(normdirect,normdirect)=-1.0;
        }
        else
        {
          xjm_(normdirect,normdirect)= 1.0;
        }

        const double det =
          xjm_(0,0)*xjm_(1,1)*xjm_(2,2)
          +
          xjm_(0,1)*xjm_(1,2)*xjm_(2,0)
          +
          xjm_(0,2)*xjm_(1,0)*xjm_(2,1)
          -
          xjm_(0,2)*xjm_(1,1)*xjm_(2,0)
          -
          xjm_(0,0)*xjm_(1,2)*xjm_(2,1)
          -
          xjm_(0,1)*xjm_(1,0)*xjm_(2,2);

        // check for degenerated elements
        if (det <= 0.0)
        {
          dserror("GLOBAL ELEMENT NO.%i\nNEGATIVE JACOBIAN DETERMINANT: %f", ele->Id(), det);
        }

        //interpolated values at gausspoints
        double ugp=0;
        double vgp=0;
        double wgp=0;
        double pgp=0;

        // the computation of this jacobian determinant from the 3d
        // mapping is based on the assumption that we do not deform
        // our elements in wall normal direction!
        const double fac=det*iquad.Weight();

        // increase area of cutting plane in element
        area += fac;

        for(int inode=0;inode<nen_;inode++)
        {
          ugp += funct_(inode)*sol(inode*4  );
          vgp += funct_(inode)*sol(inode*4+1);
          wgp += funct_(inode)*sol(inode*4+2);
          pgp += funct_(inode)*sol(inode*4+3);
        }

        // add contribution to integral
        ubar   += ugp*fac;
        vbar   += vgp*fac;
        wbar   += wgp*fac;
        pbar   += pgp*fac;

        usqbar += ugp*ugp*fac;
        vsqbar += vgp*vgp*fac;
        wsqbar += wgp*wgp*fac;
        uvbar  += ugp*vgp*fac;
        uwbar  += ugp*wgp*fac;
        vwbar  += vgp*wgp*fac;
        psqbar += pgp*pgp*fac;
      } // end loop integration points


      // add increments from this layer to processor local vectors
      (*sumarea)[ele_cart_id[1]*numsublayers+rr] += area;

      (*sumu   )[ele_cart_id[1]*numsublayers+rr] += ubar;
      (*sumv   )[ele_cart_id[1]*numsublayers+rr] += vbar;
      (*sumw   )[ele_cart_id[1]*numsublayers+rr] += wbar;
      (*sump   )[ele_cart_id[1]*numsublayers+rr] += pbar;

      (*sumsqu )[ele_cart_id[1]*numsublayers+rr] += usqbar;
      (*sumsqv )[ele_cart_id[1]*numsublayers+rr] += vsqbar;
      (*sumsqw )[ele_cart_id[1]*numsublayers+rr] += wsqbar;
      (*sumuv  )[ele_cart_id[1]*numsublayers+rr] += uvbar;
      (*sumuw  )[ele_cart_id[1]*numsublayers+rr] += uwbar;
      (*sumvw  )[ele_cart_id[1]*numsublayers+rr] += vwbar;
      (*sumsqp )[ele_cart_id[1]*numsublayers+rr] += psqbar;
    }

  }
  else
  {
    dserror("Unknown element type for mean value evaluation\n");
  }

  return 0;
}

/*-----------------------------------------------------------------------------*
 | Calculate properties for adaptive time step based on CFL number  bk 08/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidEleCalc<distype>::CalcTimeStep(
    DRT::ELEMENTS::Fluid*                ele,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm,
    Epetra_SerialDenseVector&       elevec1
    )
{
  // ---------------------------------------------------------------------------
  // Geometry
  // ---------------------------------------------------------------------------
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // Do ALE specific updates if necessary
  if (ele->IsAle())
    dserror("The adaptive time step is not implemented for Ale flow up to now. How is the CFL number defined with Ale? Update and test yourself, it is easy.");

  // evaluate shape functions and derivatives element center
  EvalShapeFuncAndDerivsAtEleCenter();

  // np_genalpha: additional vector for velocity at time n+1
  LINALG::Matrix<nsd_,nen_> evelnp(true);

  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelnp, NULL,"velnp");

  LINALG::Matrix<nsd_,1> velint(true);
  velint.Multiply(evelnp,funct_);

  //calculate element length via the stream length definition, see corresponding implementation
  //in fluid_ele_calc for calculation of the stabilization parameter
  double h=0.0;
  double vel_norm=velint.Norm2();
  if(vel_norm>1.0e-6)
  {
    LINALG::Matrix<nsd_,1> velino(true);
    velino.Update(1.0/vel_norm,velint);

    // get streamlength using the normed velocity at element centre
    LINALG::Matrix<nen_,1> tmp;
    tmp.MultiplyTN(derxy_,velino);
    const double val = tmp.Norm1();
    h = 2.0/val; // h=streamlength

    elevec1[0] = h/vel_norm;
  }
  else
    elevec1[0] = 1.0e12;

  return 0;
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
