/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_calc_general_service.cpp

\brief general service routines for calculation of fluid element

\maintainer Volker Gravemeier

\level 1
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
#include "../drt_mat/fluid_linear_density_viscosity.H"
#include "../drt_mat/fluid_murnaghantait.H"
#include "../drt_mat/fluidporo.H"

#include "../drt_nurbs_discret/drt_nurbs_utils.H"
#include "Sacado.hpp"

// immersed fsi related
#include "../drt_immersed_problem/immersed_base.H"
#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 * Evaluate supporting methods of the element
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::EvaluateService(
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
    case FLD::calc_press_grad_and_div_eps:
    {
      // calculate pressure gradient and stress term at specified element coordinates
      return CalcPressGradAndDivEps(
          ele,
          params,
          discretization,
          lm,
          elevec1,
          elevec2);
    }
    break;
    case FLD::calc_volume_gaussint:
    {
      // calculate volume integral over the element for fluid fraction
      return ComputeVolumeIntegral(ele, params, discretization, lm, elevec1);
    }
    break;
    case FLD::calc_fluidfrac_projection:
    {
      // project element fluid fraction to nodal level
      return FluidFractionProjection(ele, params, discretization, lm, elemat1, elevec1);
    }
    break;
    case FLD::calc_mass_matrix:
    {
      // compute element mass matrix
      return CalcMassMatrix(ele, discretization, lm, mat, elemat1);
    }
    break;
    case FLD::interpolate_velgrad_to_given_point:
    {
      // interpolate velocity gradient grad(u) to given point
      return InterpolateVelocityGradientAndPressure(ele, discretization, lm, elevec1, elevec2);
    }
    break;
    case FLD::interpolate_velocity_to_given_point_immersed:
    {
      // interpolate structural velocity to given point
      return InterpolateVelocityToNode(params, ele, discretization, lm, elevec1, elevec2);
    }
    break;
    case FLD::search_immersed_boundary_elements:
    {
      // search for immersed boundary fluid elements
      return SearchImmersedBoundaryElements(params, ele, elemat1, elevec1);
    }
    case FLD::least_squares_matrix_rhs_immersed_boundary:
    {
      // get least squares matrix and rhs for immersed boundary elements
      return GetLeastSquaresMatrixImmersedBoundary(params, ele, discretization, lm, elemat1, elevec1);
    }
    break;
    case FLD::update_immersed_information:
    {
      // correct immersed velocities for fluid boundary elements
      return UpdateImmersedInformation(params, ele, discretization,lm);
    }
    break;
    case FLD::correct_immersed_fluid_bound_vel:
    {
      // correct immersed velocities for fluid boundary elements
      return CorrectImmersedBoundVelocities(params, ele, discretization, lm, mat, elevec1, elevec2);
    }
    break;
    case FLD::calc_artificial_velocity_divergence:
    {
      // interpolate structural velocity to given point
      return CalcArtificialVelocityDivergence(params, ele, discretization, lm, elevec1);
    }
    break;
    case FLD::interpolate_velocity_to_given_point:
    {
      // interpolate velocity to given point
      return InterpolateVelocityToPoint(ele, params, discretization, lm, elevec1, elevec2);
    }
    break;
    case FLD::interpolate_pressure_to_given_point:
    {
      // interpolate pressure to given point
      return InterpolatePressureToPoint(ele, params, discretization, lm, elevec1);
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
    case FLD::velgradient_projection:
    {
      // project velocity gradient to nodal level
      return VelGradientProjection(ele, params, discretization, lm, elemat1, elemat2);
    }
    break;
    case FLD::presgradient_projection:
    {
      // project velocity gradient to nodal level
      return PresGradientProjection(ele, params, discretization, lm, elemat1, elemat2);
    }
    break;
    case FLD::calc_dt_via_cfl:
    {
      return CalcTimeStep(ele, discretization, lm, elevec1);
    }
    break;
    case FLD::calc_velgrad_ele_center:
    {
      return CalcVelGradientEleCenter(ele, discretization, lm, elevec1, elevec2);
    }
    break;
    case FLD::calc_mass_flow_periodic_hill:
    {
      // compute element mass matrix
      return CalcMassFlowPeriodicHill(ele, params, discretization, lm, elevec1,mat);
    }
    break;
    case FLD::reset_immersed_ele:
    {
      return ResetImmersedEle(ele,params);
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
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::IntegrateShapeFunction(
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
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::IntegrateShapeFunction(
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
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcDivOp(
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
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcMatDerivAndRotU(
  DRT::ELEMENTS::Fluid*     ele,
  Teuchos::ParameterList&   params,
  DRT::Discretization&      discretization,
  std::vector<int>&         lm,
  Epetra_SerialDenseVector& elevec1,
  Epetra_SerialDenseVector& elevec2,
  Epetra_SerialDenseVector& elevec3)
{
  // fill the local element vector/matrix with the global values
  LINALG::Matrix<nsd_,nen_> evel(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evel, NULL,"vel");
  LINALG::Matrix<nsd_,nen_> eacc(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &eacc, NULL,"acc");

  // coordinates of the current integration point
  xsi_.Update(1.0, params.get<LINALG::Matrix<nsd_,1> >("elecoords"));

  // shape functions and their first derivatives
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

  // get velocities u_n and u_nm at integration point
  velint_.Multiply(evel,funct_);
  LINALG::Matrix<nsd_,1> accint(true);
  accint.Multiply(eacc,funct_);

  for (int isd=0;isd<nsd_;isd++)
  {
    elevec1[isd] = velint_(isd);
  }

  // get gradient of velocity at integration point
  vderxy_.MultiplyNT(evel,deriv_);

  // calculate (u_n * nabla) u_n
  conv_old_.Multiply(vderxy_,velint_);

  // calculate (u_n - u_nm)/dt + (u_n * nabla) u_n
  conv_old_.Update(1.0, accint, 1.0);

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
 | Action type: calc_press_grad_and_div_eps                            |
 | calculate pressure gradient and stress term             ghamm 06/14 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcPressGradAndDivEps(
    DRT::ELEMENTS::Fluid*     ele,
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2)
{
  //--------------------------------------------------------------------------------
  // extract element based or nodal values
  //--------------------------------------------------------------------------------

  // get velocity gradient at the nodes
  const Teuchos::RCP<Epetra_MultiVector> velocity_gradient = params.get< Teuchos::RCP<Epetra_MultiVector> >("velgradient");

  const int nsdsquare = nsd_*nsd_;

  Epetra_SerialDenseVector evelgrad(nen_*nsdsquare);
  DRT::UTILS::ExtractMyNodeBasedValues(ele,evelgrad,velocity_gradient,nsdsquare);

  // insert into element arrays
  for (int i=0;i<nen_;++i)
  {
    // insert velocity gradient field into element array
    for (int idim=0 ; idim < nsdsquare; ++idim)
    {
      viscs2_(idim,i) = evelgrad[idim + i*nsdsquare];
    }
  }

  //----------------------------------------------------------------------------
  //   Extract velocity/pressure from global vectors
  //----------------------------------------------------------------------------

  // fill the local element vector with the global values
  LINALG::Matrix<nen_,1>    epre(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, NULL, &epre,"vel");

  // coordinates of the current integration point
  LINALG::Matrix<nsd_,1> elecoords = params.get<LINALG::Matrix<nsd_,1> >("elecoords");

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
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"disp");

    // get new node positions for isale
     xyze_ += edispnp;
  }

  // the int point considered is the point given from outside
  EvalShapeFuncAndDerivsAtIntPoint(elecoords.A(), -1.0);

  // elevec1 contains the pressure gradient
  gradp_.Multiply(derxy_,epre);
  for (int isd=0; isd<nsd_; ++isd)
  {
    elevec1[isd] = gradp_(isd);
  }

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

  // get second derivatives w.r.t. xyz of velocity at given point
  LINALG::Matrix<nsd_*nsd_,nsd_> evelgrad2(true);
  evelgrad2.MultiplyNT(viscs2_,derxy_);

  /*--- evelgrad2 --------------------------------*/
  /*
     /                        \
     |   u,xx   u,xy   u,xz   |
     |   u,yx   u,yy   u,yz   |
     |   u,zx   u,zy   u,zz   |
     |                        |
     |   v,xx   v,xy   v,xz   |
     |   v,yx   v,yy   v,yz   |
     |   v,zx   v,zy   v,zz   |
     |                        |
     |   w,xx   w,xy   w,xz   |
     |   w,yx   w,yy   w,yz   |
     |   w,zx   w,zy   w,zz   |
     \                        /
                                                  */

  // elevec2 contains div(eps(u))
  elevec2[0] = 2.0 * evelgrad2(0,0) + evelgrad2(1,1) + evelgrad2(3,1) + evelgrad2(2,2) + evelgrad2(6,2);
  elevec2[1] = evelgrad2(3,0) + evelgrad2(1,0) + 2.0 * evelgrad2(4,1) + evelgrad2(5,2) + evelgrad2(7,2);
  elevec2[2] = evelgrad2(6,0) + evelgrad2(2,0) + evelgrad2(7,1) + evelgrad2(5,1) + 2.0 * evelgrad2(8,2);
  elevec2.Scale(0.5);

  // subtraction of div((1/3)*(div u)*I)
  const double onethird = 1.0/3.0;
  elevec2[0] -= onethird * (evelgrad2(0,0) + evelgrad2(4,0) + evelgrad2(8,0));
  elevec2[1] -= onethird * (evelgrad2(0,1) + evelgrad2(4,1) + evelgrad2(8,1));
  elevec2[2] -= onethird * (evelgrad2(0,2) + evelgrad2(4,2) + evelgrad2(8,2));

  return 0;
}


/*---------------------------------------------------------------------*
 | Action type: calc_volume_gaussint                                   |
 | calculate volume integral over the element for void                 |
 | fraction                                                ghamm 01/13 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::ComputeVolumeIntegral(
  DRT::ELEMENTS::Fluid*     ele,
  Teuchos::ParameterList&   params,
  DRT::Discretization&      discretization,
  std::vector<int>&         lm,
  Epetra_SerialDenseVector& elevec1)
{
  const double influence = params.get<double>("influence");
  LINALG::Matrix<3,1> pos = params.get<LINALG::Matrix<3,1> >("particlepos");

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
  if (nsd_ == 3 && distype == DRT::Element::hex8)
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
    }
  }
  else if (nsd_ == 3 && distype == DRT::Element::tet4)
  {
    const int gp_per_dir = params.get<int>("gp_per_dir");
    switch(gp_per_dir)
    {
    case 1:
      gaussrule = DRT::UTILS::intrule_tet_1point;
    break;
    case 2:
      gaussrule = DRT::UTILS::intrule_tet_4point;
    break;
    case 3:
      gaussrule = DRT::UTILS::intrule_tet_5point;
    break;
    case 4:
      gaussrule = DRT::UTILS::intrule_tet_11point;
    break;
    case 5:
      gaussrule = DRT::UTILS::intrule_tet_24point;
    break;
    case 6:
      gaussrule = DRT::UTILS::intrule_tet_45point;
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

    const double integrandtimesfac = integrand*fac_;
    // sum void fraction over all gauss points
    for (int i=0; i<nen_; ++i)
    {
      elevec1[0] += integrandtimesfac*funct_(i);
    }

  }

  return 0;
}


/*---------------------------------------------------------------------*
 | Action type: calc_fluidfrac_projection                               |
 | project element fluid fraction to nodal level            ghamm 04/14 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::FluidFractionProjection(
  DRT::ELEMENTS::Fluid*     ele,
  Teuchos::ParameterList&   params,
  DRT::Discretization&      discretization,
  std::vector<int>&         lm,
  Epetra_SerialDenseMatrix& elemat1,
  Epetra_SerialDenseVector& elevec1)
{
  // get fluid fraction of element
  double elefluidfrac = params.get<double>("elefluidfrac");

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

  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    // fill element matrix (mass matrix) and elevec in each node
    for (int vi=0; vi<nen_; ++vi) // loop rows
    {
      for (int ui=0; ui<nen_; ++ui) // loop columns
      {
        elemat1(vi,ui) += fac_*funct_(ui)*funct_(vi);
      }

     elevec1(vi) += fac_*funct_(vi)*elefluidfrac;
    }
  } // end of integration loop

  return 0;
}


/*---------------------------------------------------------------------*
 | Action type: velgradient_projection                                 |
 | project velocity gradient to nodal level                ghamm 06/14 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::VelGradientProjection(
  DRT::ELEMENTS::Fluid*     ele,
  Teuchos::ParameterList&   params,
  DRT::Discretization&      discretization,
  std::vector<int>&         lm,
  Epetra_SerialDenseMatrix& elemat1,
  Epetra_SerialDenseMatrix& elemat2)
{
  //----------------------------------------------------------------------------
  //   Extract velocity/pressure from global vectors
  //----------------------------------------------------------------------------

  LINALG::Matrix<nsd_,nen_> evel(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evel, NULL,"vel");

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
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"disp");

    // get new node positions for isale
     xyze_ += edispnp;
  }

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------

  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    vderxy_.MultiplyNT(evel,derxy_);

    // fill element matrix (mass matrix) and elevec in each node
    for (int vi=0; vi<nen_; ++vi) // loop rows
    {
      for (int ui=0; ui<nen_; ++ui) // loop columns
      {
        elemat1(vi,ui) += fac_*funct_(ui)*funct_(vi);
      }
    }

    // fill elemat which is of size node x (nsd_*nsd_) with rhs
    for (int i=0; i<nsd_; ++i) // loop rows of vderxy
    {
      for (int j=0; j<nsd_; ++j) // loop columns of vderxy
      {
        for (int vi=0; vi<nen_; ++vi) // loop nodes
        {
          elemat2(vi,i*nsd_+j) += fac_*funct_(vi)*vderxy_(i,j);
        }
      }
    }

  } // end of integration loop

  return 0;
}

/*---------------------------------------------------------------------*
 | Action type: presgradient_projection                                |
 | project pressure gradient to nodal level               winter 09/15 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::PresGradientProjection(
  DRT::ELEMENTS::Fluid*     ele,
  Teuchos::ParameterList&   params,
  DRT::Discretization&      discretization,
  std::vector<int>&         lm,
  Epetra_SerialDenseMatrix& elemat1,
  Epetra_SerialDenseMatrix& elemat2)
{
  //----------------------------------------------------------------------------
  //   Extract velocity/pressure from global vectors
  //----------------------------------------------------------------------------

  LINALG::Matrix<nen_,1> epres(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, NULL, &epres,"pres");

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
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"disp");

    // get new node positions for isale
     xyze_ += edispnp;
  }

  //------------------------------------------------------------------
  //                       INTEGRATION LOOP
  //------------------------------------------------------------------

  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    gradp_.Multiply(derxy_,epres);

    // fill element matrix (mass matrix) and elevec in each node
    for (int vi=0; vi<nen_; ++vi) // loop rows
    {
      for (int ui=0; ui<nen_; ++ui) // loop columns
      {
        elemat1(vi,ui) += fac_*funct_(ui)*funct_(vi);
      }
    }

    // fill elemat which is of size node x (1*nsd_) with rhs
    for (int i=0; i<nsd_; ++i) // loop rows of vderxy
    {
//      for (int j=0; j<nsd_; ++j) // loop columns of vderxy
//      {
        for (int vi=0; vi<nen_; ++vi) // loop nodes
        {
          elemat2(vi,i) += fac_*funct_(vi)*gradp_(i,0);
        }
//      }
    }

  } // end of integration loop

  return 0;
}


/*----------------------------------------------------------------------*
 * Action type: Compute Div u                                 ehrl 12/12|
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::ComputeDivU(
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
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::ComputeError(
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
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::ComputeError(
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
  const int calcerrfunctno = DRT::INPUT::get<int>(params,"error function number");

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

    EvaluateAnalyticSolutionPoint(xyzint, fldparatimint_->Time(), calcerr, calcerrfunctno, mat, u, p, dervel, fldparatimint_->IsFullImplPressureAndCont(), fldparatimint_->Dt());

    if (calcerr == INPAR::FLUID::topoptchannel &&
        !(xyzint(1)>-0.2-1.0e-014 && xyzint(1)<0.2+1.0e-014))
    {
      preint = 0.0;
      u(0) = 0.0;
      u(1) = 0.0;
    }

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


template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::EvaluateAnalyticSolutionPoint (
      const LINALG::Matrix<nsd_,1>      &xyzint,
      const double                      t,
      const INPAR::FLUID::CalcError     calcerr,
      const int                         calcerrfunctno,
      const Teuchos::RCP<MAT::Material> &mat,
      LINALG::Matrix<nsd_,1>            &u,
      double                            &p,
      LINALG::Matrix<nsd_,nsd_>         &dervel,
      bool                              isFullImplPressure,
      double                            deltat
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
      if(not isFullImplPressure)
      {
        p = -a*a/2.0 *
            ( std::exp(2.0*a*xyzint(0))
        + std::exp(2.0*a*xyzint(1))
        + std::exp(2.0*a*xyzint(2))
        + 2.0 * std::sin(a*xyzint(0) + d*xyzint(1)) * std::cos(a*xyzint(2) + d*xyzint(0)) * std::exp(a*(xyzint(1)+xyzint(2)))
        + 2.0 * std::sin(a*xyzint(1) + d*xyzint(2)) * std::cos(a*xyzint(0) + d*xyzint(1)) * std::exp(a*(xyzint(2)+xyzint(0)))
        + 2.0 * std::sin(a*xyzint(2) + d*xyzint(0)) * std::cos(a*xyzint(1) + d*xyzint(2)) * std::exp(a*(xyzint(0)+xyzint(1)))
            )* std::exp(-2.0*visc*d*d*t);
      }
      else //pressure for full implicit OST scheme:
      {
        p = -a*a/2.0 *
            ( std::exp(2.0*a*xyzint(0))
        + std::exp(2.0*a*xyzint(1))
        + std::exp(2.0*a*xyzint(2))
        + 2.0 * std::sin(a*xyzint(0) + d*xyzint(1)) * std::cos(a*xyzint(2) + d*xyzint(0)) * std::exp(a*(xyzint(1)+xyzint(2)))
        + 2.0 * std::sin(a*xyzint(1) + d*xyzint(2)) * std::cos(a*xyzint(0) + d*xyzint(1)) * std::exp(a*(xyzint(2)+xyzint(0)))
        + 2.0 * std::sin(a*xyzint(2) + d*xyzint(0)) * std::cos(a*xyzint(1) + d*xyzint(2)) * std::exp(a*(xyzint(0)+xyzint(1)))
            )* std::exp(-2.0*visc*d*d*(t-0.5*deltat));
      }

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
      p = (xyzint(0)-0.5)*(-50*visc);

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
      const double u_exact_x = DRT::Problem::Instance()->Funct(0).Evaluate(0,position,t);
      const double u_exact_y = DRT::Problem::Instance()->Funct(0).Evaluate(1,position,t);
      u(0) = u_exact_x;
      u(1) = u_exact_y;
    }

  }
  break;
  case INPAR::FLUID::byfunct:
  {
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
      const double u_exact_x = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(0,position,t);
      const double u_exact_y = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(1,position,t);
      const double p_exact   = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(2,position,t);

      u(0) = u_exact_x;
      u(1) = u_exact_y;
      p    = p_exact;

      std::vector<double> uder_exact_x = DRT::Problem::Instance()->Funct(calcerrfunctno-1).EvaluateSpatialDerivative(0,position,t);
      std::vector<double> uder_exact_y = DRT::Problem::Instance()->Funct(calcerrfunctno-1).EvaluateSpatialDerivative(1,position,t);
      //std::vector<double> pder_exact   = DRT::Problem::Instance()->Funct(func_no-1).EvaluateSpatialDerivative(2,position,t,1);

      if(uder_exact_x.size())
      {
        dervel(0,0)=uder_exact_x[0];
        dervel(0,1)=uder_exact_x[1];
      }

      if(uder_exact_y.size())
      {
        dervel(1,0)=uder_exact_y[0];
        dervel(1,1)=uder_exact_y[1];
      }
    }
    else if(nsd_==3)
    {
      const double u_exact_x = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(0,position,t);
      const double u_exact_y = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(1,position,t);
      const double u_exact_z = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(2,position,t);
      const double p_exact   = DRT::Problem::Instance()->Funct(calcerrfunctno-1).Evaluate(3,position,t);

      u(0) = u_exact_x;
      u(1) = u_exact_y;
      u(2) = u_exact_z;
      p    = p_exact;

      std::vector<double> uder_exact_x = DRT::Problem::Instance()->Funct(calcerrfunctno-1).EvaluateSpatialDerivative(0,position,t);
      std::vector<double> uder_exact_y = DRT::Problem::Instance()->Funct(calcerrfunctno-1).EvaluateSpatialDerivative(1,position,t);
      std::vector<double> uder_exact_z = DRT::Problem::Instance()->Funct(calcerrfunctno-1).EvaluateSpatialDerivative(2,position,t);

      if(uder_exact_x.size())
      {
        dervel(0,0)=uder_exact_x[0];
        dervel(0,1)=uder_exact_x[1];
        dervel(0,2)=uder_exact_x[2];
      }

      if(uder_exact_y.size())
      {
        dervel(1,0)=uder_exact_y[0];
        dervel(1,1)=uder_exact_y[1];
        dervel(1,2)=uder_exact_y[2];
      }

      if(uder_exact_z.size())
      {
        dervel(2,0)=uder_exact_z[0];
        dervel(2,1)=uder_exact_z[1];
        dervel(2,2)=uder_exact_z[2];
      }

    }
    else dserror("invalid dimension");

  }
  break;
  case INPAR::FLUID::fsi_fluid_pusher:
  {
    /* Since the fluid pusher solution depends only on time, but not on spatial
     * cooordinates x,y,z, we only compute the L2-error and no H1-error.
     */

    // get pointer to material in order to access density
    if (mat->MaterialType() == INPAR::MAT::m_fluid)
    {
      const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(mat.get());

      // available solution types
      enum SOLUTIONTYPE
      {
        quadratic, // d(t) = -t^2
        cubic, // d(t) = -t^3
        quartic, // d(t) = -t^4
        quintic // d(t) = -t^5
      };

      // choose solution type
      SOLUTIONTYPE solutiontype = quintic;

      // compute analytical velocities and pressure for different prescribed time curves
      // Note: consider offset of coordinate system
      switch (solutiontype)
      {
      case quadratic:
      {
        u(0) = -2.0 * t;
        u(1) = 0.0;
        u(2) = 0.0;
        p = actmat->Density() * 2.0 * (xyzint(0) + 1.5);

        break;
      }
      case cubic:
      {
        u(0) = -3.0 * t * t;
        u(1) = 0.0;
        u(2) = 0.0;
        p = actmat->Density() * 6.0 * t * (xyzint(0) + 1.5);

        break;
      }
      case quartic:
      {
        u(0) = -4.0 * t * t * t;
        u(1) = 0.0;
        u(2) = 0.0;
        p = actmat->Density() * 12.0 * t * t * (xyzint(0) + 1.5);

        break;
      }
      case quintic:
      {
        u(0) = -5.0 * t * t * t * t;
        u(1) = 0.0;
        u(2) = 0.0;
        p = actmat->Density() * 20.0 * t * t * t * (xyzint(0) + 1.5);

        break;
      }
      default:
      {
        dserror("Unknown solution type.");
        break;
      }
      }
    }
    else dserror("Material is not a Newtonian Fluid");
  }
  break;
  case INPAR::FLUID::channel_weakly_compressible:
  {
    // Steady, weakly compressible isothermal flow of a Newtonian fluid
    // with pressure-dependent density according to Murnaghan-Tait law
    // Comparison to analytical solution obtained in:
    // "New analytical solutions for weakly compressible Newtonian
    // Poiseuille flows with pressure-dependent viscosity"
    // Kostas D. Housiadas, Georgios C. Georgiou

    if (mat->MaterialType() == INPAR::MAT::m_fluid_murnaghantait)
    {
      const MAT::MurnaghanTaitFluid* actmat = static_cast<const MAT::MurnaghanTaitFluid*>(mat.get());

      if (actmat->MatParameter() != 1.0)
      {
        dserror("The analytical solution is only valid for material parameter = 1");
      }

      double x = xyzint(0);
      double y = xyzint(1);

      double length = 10.0;
      double radius = 1.0;
      double aspect_ratio = radius/length;
      double mean_velocity_channel_exit = 1.0;
      double viscosity = actmat->Viscosity();
      double reference_pressure = actmat->RefPressure();
      double reference_bulk_modulus = actmat->RefBulkModulus();
      double linear_coefficient_density = (3.0*(1.0/reference_bulk_modulus)*viscosity*length*mean_velocity_channel_exit)/std::pow(radius,2.0);

      if (nsd_ == 2)
      {
        u(0) = 3.0/2.0*(1-std::pow(y/radius,2.0))*(1.0+linear_coefficient_density*(x/length-1.0));
        u(1) = 0.0;
        p    = 1.0-x/length-linear_coefficient_density*(1.0/6.0*std::pow(aspect_ratio,2.0)*(std::pow(y/radius,2.0)-1.0)+1.0/2.0*std::pow(1.0-x/length,2.0));
        dervel(0,0) = 3.0/2.0*(1.0-std::pow(y/radius,2.0))*linear_coefficient_density/length;
        dervel(0,1) = -3.0/std::pow(radius,2.0)*y*(1.0+linear_coefficient_density*(x/length-1.0));
        dervel(1,0) = 0.0;
        dervel(1,1) = 0.0;

        // scaling correctly the variables
        u(0) = u(0)*mean_velocity_channel_exit;
        u(1) = u(1)*mean_velocity_channel_exit*radius/length;
        p    = p*(3.0*viscosity*length*mean_velocity_channel_exit/std::pow(radius,2.0))+reference_pressure;
        dervel(0,0) = dervel(0,0)*mean_velocity_channel_exit;
        dervel(0,1) = dervel(0,1)*mean_velocity_channel_exit;
        dervel(1,0) = dervel(1,0)*mean_velocity_channel_exit*radius/length;
        dervel(1,1) = dervel(1,1)*mean_velocity_channel_exit*radius/length;
      }
      else
        dserror("3D analytical solution is not implemented");
    }
    else if (mat->MaterialType() == INPAR::MAT::m_fluid_linear_density_viscosity)
    {
      const MAT::LinearDensityViscosity* actmat = static_cast<const MAT::LinearDensityViscosity*>(mat.get());

      double x = xyzint(0);
      double y = xyzint(1);

      double length = 10.0;
      double radius = 1.0;
      double aspect_ratio = radius/length;
      double mean_velocity_channel_exit = 1.0;
      double reference_viscosity = actmat->RefViscosity();
      double reference_pressure = actmat->RefPressure();
      double coefficient_density = actmat->CoeffDensity();
      double coefficient_viscosity = actmat->CoeffViscosity();
      double coefficient_density_adim = (3.0*coefficient_density*reference_viscosity*length*mean_velocity_channel_exit)/std::pow(radius,2.0);
      double coefficient_viscosity_adim = (3.0*coefficient_viscosity*reference_viscosity*length*mean_velocity_channel_exit)/std::pow(radius,2.0);

      // parameters according with the paper
      double z = x / length;
      double r = y / radius;
      double alfa = aspect_ratio;
      double beta = coefficient_viscosity_adim;
      double epsilon = coefficient_density_adim;
      double a = aspect_ratio;
      double B = alfa*beta;
      double lambda =     1.0                               +
                      (   1.0 /      5.0) * std::pow(B,2.0) +
                      (  11.0 /    175.0) * std::pow(B,4.0) +
                      ( 533.0 /  23625.0) * std::pow(B,6.0) +
                      (5231.0 / 606375.0) * std::pow(B,8.0);
      double p_0_hat = std::cosh(alfa * beta * lambda * r) /
                       std::cosh(alfa * beta * lambda);
      double u_r1_hat = -(11.0 * r * std::pow(1.0 - std::pow(r,2.0),2.0)) / 40.0 * std::pow(B,2.0) *
                         (
                                  1.0                                                                                                                                 +
                          ((    173.0 -       85.0 * std::pow(r,2.0))                                                              /       (770.0)) * std::pow(B,2.0) +
                          ((   5793.0 -     7190.0 * std::pow(r,2.0) +     3965.0 * std::pow(r,4.0))                               /     (83160.0)) * std::pow(B,4.0) +
                          ((7435723.0 - 16839665.0 * std::pow(r,2.0) + 16836225.0 * std::pow(r,4.0) - 5021275.0 * std::pow(r,6.0)) / (320166000.0)) * std::pow(B,6.0)
                         );
      double u_r1_hat_first = (11.0*std::pow(B,2.0)*std::pow(std::pow(r,2.0)-1.0,2.0)*(((4099.0*std::pow(r,6.0))/261360.0-(32069.0*std::pow(r,4.0))/
          609840.0+(3367933.0*std::pow(r,2.0))/64033200.0-7435723.0/320166000.0)*std::pow(B,6.0)+(-(793.0*std::pow(r,4.0))/
          16632.0+(719.0*std::pow(r,2.0))/8316.0-1931.0/27720.0)*std::pow(B,4.0)+((17.0*std::pow(r,2.0))/154.0-173.0/770.0)*std::pow(B,2.0)-1.0))/40.0 +
          (11.0*std::pow(B,2.0)*std::pow(r,2.0)*(std::pow(r,2.0)-1.0)*(((4099.0*std::pow(r,6.0))/261360.0-(32069.0*std::pow(r,4.0))/609840.0+
          (3367933.0*std::pow(r,2.0))/64033200.0-7435723.0/320166000.0)*std::pow(B,6.0)+(-(793.0*std::pow(r,4.0))/16632.0+(719.0*
          std::pow(r,2.0))/8316.0-1931.0/27720.0)*std::pow(B,4.0)+((17.0*std::pow(r,2.0))/154.0-173.0/770.0)*std::pow(B,2.0)-1.0))/10.0 +
          (11.0*std::pow(B,2.0)*
          r*std::pow(std::pow(r,2.0)-1.0,2.0)*(((4099.0*std::pow(r,5.0))/43560.0-(32069.0*std::pow(r,3.0))/152460.0+(3367933.0*r)/
          32016600.0)*std::pow(B,6.0)+((719.0*r)/4158.0-(793.0*std::pow(r,3.0))/4158.0)*std::pow(B,4.0)+(17.0*r*std::pow(B,2.0))/77.0))/40.0;
      double h = 1.0 / std::pow(beta,2.0) *
                 (
                        -1.0                                                                                                                                                         +
                  ((    11.0 -     10.0 * std::pow(r,2.0))                                                                                          /      (15.0)) * std::pow(B,2.0) +
                  ((   359.0 -    126.0 * std::pow(r,2.0) +      35.0 * std::pow(r,4.0))                                                            /    (1260.0)) * std::pow(B,4.0) +
                  (( 13761.0 -  17790.0 * std::pow(r,2.0) +   34125.0 * std::pow(r,4.0) -   17500.0 * std::pow(r,6.0))                              /   (94500.0)) * std::pow(B,6.0) +
                  ((225311.0 - 614515.0 * std::pow(r,2.0) + 1492755.0 * std::pow(r,4.0) - 1324785.0 * std::pow(r,6.0) + 394350.0 * std::pow(r,8.0)) / (3118500.0)) * std::pow(B,8.0)
                 );
      double h_1 = 1.0 / std::pow(beta,2.0) *
                 (
                        -1.0                                                                                 +
                  ((    11.0 -     10.0)                                    /      (15.0)) * std::pow(B,2.0) +
                  ((   359.0 -    126.0 +      35.0)                        /    (1260.0)) * std::pow(B,4.0) +
                  (( 13761.0 -  17790.0 +   34125.0 -   17500.0)            /   (94500.0)) * std::pow(B,6.0) +
                  ((225311.0 - 614515.0 + 1492755.0 - 1324785.0 + 394350.0) / (3118500.0)) * std::pow(B,8.0)
                 );

      if (nsd_ == 2)
      {
        u(0) = - (3.0 * std::log(p_0_hat)) / (std::pow(B,2.0) * lambda) +
               epsilon *
               (
                (3 * (std::tanh(B * lambda) - r * std::tanh(B * lambda * r)) + std::log(std::pow(p_0_hat,3.0)) / (B * lambda)) /
                (beta * (3 * std::tanh(B * lambda) - 2 * B)) +
                (std::exp(lambda * beta * (1.0 - z))) / (lambda * beta) *
                ((p_0_hat * std::log(std::pow(p_0_hat,3.0))) / (std::pow(B,2.0)) + u_r1_hat_first)
               );
        u(1) = epsilon * u_r1_hat * std::exp(lambda * beta * (1.0 - z));
        p    = (p_0_hat * std::exp(lambda * beta * (1.0 - z)) - 1.0) / beta +
               epsilon * p_0_hat * std::exp(lambda * beta * (1.0 - z)) *
               (
                (lambda * a * (1.0 - z + a * (r * std::tanh(B * lambda * r) - std::tanh(B * lambda)))) /
                (3.0 * std::tanh(B * lambda) - 2.0 * B) +
                p_0_hat * h * std::exp(lambda * beta * (1.0- z )) -
                h_1
               );

        // scaling correctly the variables
        u(0) = u(0)*mean_velocity_channel_exit;
        u(1) = u(1)*mean_velocity_channel_exit*radius/length;
        p    = p*(3.0*reference_viscosity*length*mean_velocity_channel_exit/std::pow(radius,2.0))+reference_pressure;
      }
      else
        dserror("3D analytical solution is not implemented");
    }
    else
    {
      dserror("The analytical solution is not implemented for this material");
    }

  }
  break;
  case INPAR::FLUID::channel_weakly_compressible_fourier_3:
    {
      // Steady, weakly compressible isothermal flow of a Newtonian fluid
      // with pressure-dependent density according to Murnaghan-Tait law
      // Comparison to analytical solution obtained in:
      // "New analytical solutions for weakly compressible Newtonian
      // Poiseuille flows with pressure-dependent viscosity"
      // Kostas D. Housiadas, Georgios C. Georgiou
      // Solution obtained with a Fourier expansion up to 3rd order

      if (mat->MaterialType() == INPAR::MAT::m_fluid_murnaghantait)
      {
        const MAT::MurnaghanTaitFluid* actmat = static_cast<const MAT::MurnaghanTaitFluid*>(mat.get());

        if (actmat->MatParameter() != 1.0)
        {
          dserror("The analytical solution is only valid for material parameter = 1");
        }

        double x = xyzint(0);
        double y = xyzint(1);

        double L = 10.0;
        double R = 1.0;
        double U = 1.0;
        double mu = actmat->Viscosity();
        double p0 = actmat->RefPressure();
        double K0 = actmat->RefBulkModulus();

        if (nsd_ == 2)
        {
          u(0) = (2.0*L*R*U-(3.0*(std::pow(L,2.0))*(std::pow(U,2.0))*mu)/(K0*R))/(2.0*L*R)+(3.0*U*std::cos((y*M_PI)/R)*(2.0*K0*(std::pow(R,2.0))-3.0*L*U*mu))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))-(3.0*U*std::cos((2.0*y*M_PI)/R)*(2.0*K0*(std::pow(R,2.0))-3.0*L*U*mu))/(4.0*K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))+(U*std::cos((3.0*y*M_PI)/R)*(2.0*K0*(std::pow(R,2.0))-3.0*L*U*mu))/(3.0*K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))-(3.0*L*(std::pow(U,2.0))*mu*std::sin((2.0*x*M_PI)/L))/(K0*(std::pow(R,2.0))*M_PI)-(3.0*L*(std::pow(U,2.0))*mu*std::sin((4.0*x*M_PI)/L))/(2.0*K0*(std::pow(R,2.0))*M_PI)-(L*(std::pow(U,2.0))*mu*std::sin((6.0*x*M_PI)/L))/(K0*(std::pow(R,2.0))*M_PI)-(18.0*L*(std::pow(U,2.0))*mu*std::cos((y*M_PI)/R)*std::sin((2.0*x*M_PI)/L))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,3.0)))-(9.0*L*(std::pow(U,2.0))*mu*std::cos((y*M_PI)/R)*std::sin((4.0*x*M_PI)/L))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,3.0)))+(9.0*L*(std::pow(U,2.0))*mu*std::cos((2.0*y*M_PI)/R)*std::sin((2.0*x*M_PI)/L))/(2.0*K0*(std::pow(R,2.0))*(std::pow(M_PI,3.0)))-(2.0*L*(std::pow(U,2.0))*mu*std::cos((3.0*y*M_PI)/R)*std::sin((2.0*x*M_PI)/L))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,3.0)))-(6.0*L*(std::pow(U,2.0))*mu*std::cos((y*M_PI)/R)*std::sin((6.0*x*M_PI)/L))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,3.0)))+(9.0*L*(std::pow(U,2.0))*mu*std::cos((2.0*y*M_PI)/R)*std::sin((4.0*x*M_PI)/L))/(4.0*K0*(std::pow(R,2.0))*(std::pow(M_PI,3.0)))-(L*(std::pow(U,2.0))*mu*std::cos((3.0*y*M_PI)/R)*std::sin((4.0*x*M_PI)/L))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,3.0)))+(3.0*L*(std::pow(U,2.0))*mu*std::cos((2.0*y*M_PI)/R)*std::sin((6.0*x*M_PI)/L))/(2.0*K0*(std::pow(R,2.0))*(std::pow(M_PI,3.0)))-(2.0*L*(std::pow(U,2.0))*mu*std::cos((3.0*y*M_PI)/R)*std::sin((6.0*x*M_PI)/L))/(3.0*K0*(std::pow(R,2.0))*(std::pow(M_PI,3.0)));
          u(1) = 0.0;
          p    = (2.0*L*R*p0+(L*(std::pow(R,2.0))*(2.0*(std::pow(U,2.0))*(std::pow(mu,2.0))+3.0*K0*L*U*mu)-3.0*(std::pow(L,3.0))*(std::pow(U,2.0))*(std::pow(mu,2.0)))/(K0*(std::pow(R,3.0))))/(2.0*L*R)+(6.0*(std::pow(U,2.0))*(std::pow(mu,2.0))*std::cos((y*M_PI)/R))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))-(3.0*(std::pow(U,2.0))*(std::pow(mu,2.0))*std::cos((2.0*y*M_PI)/R))/(2.0*K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))+(2.0*(std::pow(U,2.0))*(std::pow(mu,2.0))*std::cos((3.0*y*M_PI)/R))/(3.0*K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))-(9.0*(std::pow(L,2.0))*(std::pow(U,2.0))*(std::pow(mu,2.0))*std::cos((2.0*x*M_PI)/L))/(2.0*K0*(std::pow(R,4.0))*(std::pow(M_PI,2.0)))-(9.0*(std::pow(L,2.0))*(std::pow(U,2.0))*(std::pow(mu,2.0))*std::cos((4.0*x*M_PI)/L))/(8.0*K0*(std::pow(R,4.0))*(std::pow(M_PI,2.0)))-((std::pow(L,2.0))*(std::pow(U,2.0))*(std::pow(mu,2.0))*std::cos((6.0*x*M_PI)/L))/(2.0*K0*(std::pow(R,4.0))*(std::pow(M_PI,2.0)))+(3.0*L*U*mu*std::sin((2.0*x*M_PI)/L)*(2.0*K0*(std::pow(R,2.0))-3.0*L*U*mu))/(2.0*K0*(std::pow(R,4.0))*M_PI)+(3.0*L*U*mu*std::sin((4.0*x*M_PI)/L)*(2.0*K0*(std::pow(R,2.0))-3.0*L*U*mu))/(4.0*K0*(std::pow(R,4.0))*M_PI)+(L*U*mu*std::sin((6.0*x*M_PI)/L)*(2.0*K0*(std::pow(R,2.0))-3.0*L*U*mu))/(2.0*K0*(std::pow(R,4.0))*M_PI);
          dervel(0,0) = (9.0*(std::pow(U,2.0))*mu*std::cos((2.0*x*M_PI)/L)*std::cos((2.0*y*M_PI)/R))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))-(6.0*(std::pow(U,2.0))*mu*std::cos((4.0*x*M_PI)/L))/(K0*(std::pow(R,2.0)))-(6.0*(std::pow(U,2.0))*mu*std::cos((6.0*x*M_PI)/L))/(K0*(std::pow(R,2.0)))-(36.0*(std::pow(U,2.0))*mu*std::cos((2.0*x*M_PI)/L)*std::cos((y*M_PI)/R))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))-(6.0*(std::pow(U,2.0))*mu*std::cos((2.0*x*M_PI)/L))/(K0*(std::pow(R,2.0)))-(36.0*(std::pow(U,2.0))*mu*std::cos((4.0*x*M_PI)/L)*std::cos((y*M_PI)/R))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))-(4.0*(std::pow(U,2.0))*mu*std::cos((2.0*x*M_PI)/L)*std::cos((3.0*y*M_PI)/R))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))+(9.0*(std::pow(U,2.0))*mu*std::cos((4.0*x*M_PI)/L)*std::cos((2.0*y*M_PI)/R))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))-(36.0*(std::pow(U,2.0))*mu*std::cos((6.0*x*M_PI)/L)*std::cos((y*M_PI)/R))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))-(4.0*(std::pow(U,2.0))*mu*std::cos((4.0*x*M_PI)/L)*std::cos((3.0*y*M_PI)/R))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))+(9.0*(std::pow(U,2.0))*mu*std::cos((6.0*x*M_PI)/L)*std::cos((2.0*y*M_PI)/R))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)))-(4.0*(std::pow(U,2.0))*mu*std::cos((6.0*x*M_PI)/L)*std::cos((3.0*y*M_PI)/R))/(K0*(std::pow(R,2.0))*(std::pow(M_PI,2.0)));
          dervel(0,1) = (3.0*U*std::sin((2.0*y*M_PI)/R)*(2.0*K0*(std::pow(R,2.0))-3.0*L*U*mu))/(2.0*K0*(std::pow(R,3.0))*M_PI)-(3.0*U*std::sin((y*M_PI)/R)*(2.0*K0*(std::pow(R,2.0))-3.0*L*U*mu))/(K0*(std::pow(R,3.0))*M_PI)-(U*std::sin((3.0*y*M_PI)/R)*(2.0*K0*(std::pow(R,2.0))-3.0*L*U*mu))/(K0*(std::pow(R,3.0))*M_PI)+(18.0*L*(std::pow(U,2.0))*mu*std::sin((2.0*x*M_PI)/L)*std::sin((y*M_PI)/R))/(K0*(std::pow(R,3.0))*(std::pow(M_PI,2.0)))-(9.0*L*(std::pow(U,2.0))*mu*std::sin((2.0*x*M_PI)/L)*std::sin((2.0*y*M_PI)/R))/(K0*(std::pow(R,3.0))*(std::pow(M_PI,2.0)))+(9.0*L*(std::pow(U,2.0))*mu*std::sin((4.0*x*M_PI)/L)*std::sin((y*M_PI)/R))/(K0*(std::pow(R,3.0))*(std::pow(M_PI,2.0)))+(6.0*L*(std::pow(U,2.0))*mu*std::sin((2.0*x*M_PI)/L)*std::sin((3.0*y*M_PI)/R))/(K0*(std::pow(R,3.0))*(std::pow(M_PI,2.0)))-(9.0*L*(std::pow(U,2.0))*mu*std::sin((4.0*x*M_PI)/L)*std::sin((2.0*y*M_PI)/R))/(2.0*K0*(std::pow(R,3.0))*(std::pow(M_PI,2.0)))+(6.0*L*(std::pow(U,2.0))*mu*std::sin((6.0*x*M_PI)/L)*std::sin((y*M_PI)/R))/(K0*(std::pow(R,3.0))*(std::pow(M_PI,2.0)))+(3.0*L*(std::pow(U,2.0))*mu*std::sin((4.0*x*M_PI)/L)*std::sin((3.0*y*M_PI)/R))/(K0*(std::pow(R,3.0))*(std::pow(M_PI,2.0)))-(3.0*L*(std::pow(U,2.0))*mu*std::sin((6.0*x*M_PI)/L)*std::sin((2.0*y*M_PI)/R))/(K0*(std::pow(R,3.0))*(std::pow(M_PI,2.0)))+(2.0*L*(std::pow(U,2.0))*mu*std::sin((6.0*x*M_PI)/L)*std::sin((3.0*y*M_PI)/R))/(K0*(std::pow(R,3.0))*(std::pow(M_PI,2.0)));
          dervel(1,0) = 0.0;
          dervel(1,1) = 0.0;
        }
        else
          dserror("3D analytical solution is not implemented");
      }
      else
      {
        dserror("The analytical solution is not implemented for this material");
      }
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
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::ExtractValuesFromGlobalVector( const DRT::Discretization&   discretization, ///< discretization
                                    const std::vector<int>&      lm,             ///<
                                    FLD::RotationallySymmetricPeriodicBC<distype,nsd_+1,enrtype> & rotsymmpbc, ///<
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
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcDissipation(
  Fluid*                     ele,
  Teuchos::ParameterList&    params,
  DRT::Discretization&       discretization,
  std::vector<int>&          lm,
  Teuchos::RCP<MAT::Material> mat)
{
  //----------------------------------------------------------------------
  // get all nodal values
  // ---------------------------------------------------------------------
  if(not fldparatimint_->IsGenalpha())
    dserror("this routine supports only GenAlpha currently");
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
  // af_genalpha: velocity/pressure at time n+alpha_F and n+alpha_M
  // np_genalpha: velocity at time n+alpha_F, pressure at time n+1
  // ost:         velocity/pressure at time n+1
  LINALG::Matrix<nsd_,nen_> evelaf(true);
  LINALG::Matrix<nen_,1>    epreaf(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelaf, &epreaf,"velaf");

  LINALG::Matrix<nsd_,nen_> evelam(true);
  LINALG::Matrix<nen_,1>    epream(true);
  if (fldpara_->PhysicalType() == INPAR::FLUID::weakly_compressible && fldparatimint_->IsGenalpha())
  {
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelam, &epream,"velam");
  }
  if (fldpara_->PhysicalType() == INPAR::FLUID::weakly_compressible_stokes && fldparatimint_->IsGenalpha())
  {
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelam, &epream,"velam");
  }

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

  if (fldpara_->IsReconstructDer())
  {
    const Teuchos::RCP<Epetra_MultiVector> velafgrad = params.get< Teuchos::RCP<Epetra_MultiVector> >("velafgrad");
    DRT::UTILS::ExtractMyNodeBasedValues(ele,evelafgrad_,velafgrad,nsd_*nsd_);
  }

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
  for(int inode=0;inode<ele->NumNode();inode++)
    center+=xyze_(1,inode);

  center/=(double)ele->NumNode();

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

    GetMaterialParams(mat,evelaf,epreaf,epream,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam,vol);

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
        //GetMaterialParams(material,evelaf,epreaf,epream,escaaf,escaam,thermpressaf,thermpressam,thermpressdtam,vol);
        GetMaterialParams(mat,evelaf,epreaf,epream,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam,vol);

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
      GetMaterialParams(mat,evelaf,epreaf,epream,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam,vol);

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
          GetMaterialParams(mat,evelaf,epreaf,epream,escaaf,escaam,escabofoaf,thermpressaf,thermpressam,thermpressdtaf,thermpressdtam,vol);

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
          UpdateMaterialParams(mat,evelaf,epreaf,epream,escaaf,escaam,thermpressaf,thermpressam,mfssgscaint_);
        else
          UpdateMaterialParams(mat,evelaf,epreaf,epream,escaaf,escaam,thermpressaf,thermpressam,sgscaint_);
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
      // add consistency terms for MFS if applicable
      MultfracSubGridScalesConsistentResidual();
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
      // add consistency terms for MFS if applicable
      MultfracSubGridScalesConsistentResidual();
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
    // The stream length is not correct for xwall as calculated here, because the virtual (enriched) dofs don't have any geometric meaning
    if(enrtype==DRT::ELEMENTS::Fluid::xwall)
      strle=10000000;
    else
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

  Teuchos::RCP<std::vector<double> > incrmk            = params.get<Teuchos::RCP<std::vector<double> > >("incrmk"           );

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

  // element mk in stabilisation parameter
  (*incrmk       )[nlayer] += GetMK();

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
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::FDcheck(
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
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
void DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::InflowElement(DRT::Element* ele)
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
 | Calculate element mass matrix                              la spina 06/2017 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcMassMatrix(
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

  GetMaterialParams(mat, mat1, mat2, mat2, mat2, mat2, mat2, 0.0, 0.0, 0.0, 0.0, 0.0);

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

  // add terms associated to pressure dofs for weakly_compressible flows
  if (fldpara_->PhysicalType() == INPAR::FLUID::weakly_compressible)
  {
    dserror("Evaluation of the mass matrix for pressure dofs");
    // check fluid material
    if (mat->MaterialType() != INPAR::MAT::m_fluid_murnaghantait)
    {
      dserror("The evaluation of the mass matrix for pressure dofs is implemented only for Murnaghan-Tait equation of state");
    }

    // extract fluid material parameters
    const MAT::MurnaghanTaitFluid* actmat = static_cast<const MAT::MurnaghanTaitFluid*>(mat.get());
    double RefPressure           = actmat->RefPressure();        // reference pressure
    double RefBulkModulus        = actmat->RefBulkModulus();     // reference bulk modulus
    double MatParameter          = actmat->MatParameter();       // material parameter according to Murnaghan-Tait

    // evaluation of the "compressibility factor"
    double compr_fac = 1.0 / (RefBulkModulus + MatParameter * (preaf_ - RefPressure));

    // definition of matrices
    LINALG::Matrix<nen_,nen_> ppmat(true);

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
          ppmat(vi,ui) += funct_(vi)*funct_(ui)*fac_*compr_fac;
        } // end for (vi)
      } // end for (ui)
    } // end of integration loop

    // ---------------------------------------------------------------------------
    // Add pressure-pressure part to matrix
    // ---------------------------------------------------------------------------
    for (int ui=0; ui<nen_; ++ui)
    {
      const int numdof_ui = numdofpernode_*ui;

      for (int vi=0; vi<nen_; ++vi)
      {
        const int numdof_vi = numdofpernode_*vi;

        elemat1_epetra(numdof_vi+nsd_, numdof_ui+nsd_) += ppmat(vi,ui);
      } // end for (vi)
    } // end for (ui)
  }

  return 0;
}

/*-----------------------------------------------------------------------------*
 | Interpolate velocity gradient                                 rauch 05/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::InterpolateVelocityGradientAndPressure(
    DRT::ELEMENTS::Fluid*                ele,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm,
    Epetra_SerialDenseVector&            elevec1_epetra, // vectofill
    Epetra_SerialDenseVector&            elevec2_epetra  // given point in parameter space coordinates
    )
{
  // declare and initialize matrix for shapefunction evaluation
  LINALG::Matrix<nen_,1> shapefunct;
  // derivatives of shapefunctions at int point
  LINALG::Matrix<nsd_,nen_> pderiv_loc;
  // velocity gradient
  LINALG::Matrix<nsd_,nsd_> dudxi;
  // du/dxi * dxi/dx
  LINALG::Matrix<nsd_,nsd_> dudxioJinv;
  // material coord. of element
  LINALG::Matrix<nsd_,nen_> xrefe;
  // current coord. of element
  LINALG::Matrix<nsd_,nen_> xcurr;
  // element velocity at time n+1
  LINALG::Matrix<nsd_,nen_> evelnp;
  // element pressure at time n+1
  LINALG::Matrix<nen_,1> eprenp;
  // velocity at int point
  LINALG::Matrix<nsd_,1> velint;
  // pressure at int point
  LINALG::Matrix<1,1> pressint;
  // dx/dxi
  LINALG::Matrix<nsd_,nsd_>xjm;
  // dxi/dx
  LINALG::Matrix<nsd_,nsd_> xji;
  //cauchystress
  LINALG::Matrix<nsd_,nsd_> cauchystress(true);

  // get dynamic viscosity
  Teuchos::RCP<MAT::Material> currentmaterial;
  currentmaterial = ele->Material(0);
  double fluiddynamicviscosity=-1234;
  if(discretization.Name()=="fluid")
    fluiddynamicviscosity = Teuchos::rcp_dynamic_cast<MAT::NewtonianFluid>(currentmaterial)->Viscosity();
  else if(discretization.Name()=="porofluid")
    fluiddynamicviscosity = Teuchos::rcp_dynamic_cast<MAT::FluidPoro>(currentmaterial)->Viscosity();
  else
    dserror("no support for discretization guaranteed. check for valid material.");

  // determine whether fluid mesh is deformable or not
  static int isALE =
      (DRT::Problem::Instance()->ImmersedMethodParams().get<std::string>("DEFORM_BACKGROUND_MESH")=="yes");

  // resize vector to the size of the nsd_ times nsd_ independent entries of the velocity gradient du/dx and the pressure
  // causes seg fault -> therefore commented; needs to be investigated. elevec1 is directly built with a size of 10 for now
  //elevec1_epetra.Resize(nsd_*nsd_+1);

  // save point anew for safety -> check later if elevec2_epetra can be use directly
  LINALG::Matrix<nsd_,1> xi;
  for(int i=0;i<nsd_;++i)
    xi(i)=elevec2_epetra(i);

  // evaluate shapefunctions at given point in reference coordinates
  DRT::UTILS::shape_function<distype>(xi,shapefunct);
  // evaluate derivatives of element shape functions at given point in reference configuration
  DRT::UTILS::shape_function_deriv1<distype>(xi,pderiv_loc);
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> state = discretization.GetState("velnp");

#ifdef DEBUG
  if(state == Teuchos::null)
    dserror("Cannot get state vector %s", "velnp");
#endif

  if(isALE)
  {
    // update fluid displacements
    LINALG::Matrix<nsd_,nen_>       edispnp(true);
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"dispnp");

    DRT::Node** nodes = ele->Nodes();
    for (int inode=0;inode<nen_;++inode)
    {
      for (int idof=0;idof<nsd_;++idof)
      {
        const double* x = nodes[inode]->X();
        xrefe(idof,inode) = x[idof];
        xcurr(idof,inode) = x[idof]+edispnp(idof,inode);
      }
    }
  }
  else
  {
    // do not update element geometry (here xrefe=X and no xcurr present)
    {
      DRT::Node** nodes = ele->Nodes();
      for (int inode=0;inode<nen_;++inode)
      {
        for (int idof=0;idof<nsd_;++idof)
        {
          const double* x = nodes[inode]->X();
          xrefe(idof,inode) = x[idof];
        }
      }
    }
    xcurr = xrefe;
  }

  // get Jacobian matrix and determinant w.r.t. spatial configuration
  //
  // |J| = det(xjm) * det(Jmat^-1) = det(xjm) * 1/det(Jmat)
  //
  //    _                     _
  //   |  x_1,1  x_2,1  x_3,1  |           d x_i
  //   |  x_1,2  x_2,2  x_3,2  | = xjm  = --------
  //   |_ x_1,3  x_2,3  x_3,3 _|           d s_j
  //    _
  xjm.MultiplyNT(pderiv_loc,xcurr); // xcurr=xrefe -> dX/ds

  // inverse of transposed jacobian "ds/dx" (xjm) -> here: ds/dX

  //    _                     _
  //   |  s_1,1  s_2,1  s_3,1  |           d s_i
  //   |  s_1,2  s_2,2  s_3,2  | = xji  = -------- ;  [xji] o [xjm] = I
  //   |_ s_1,3  s_2,3  s_3,3 _|           d x_j
  //    _
  xji.Invert(xjm);

  // fill locationarray
  DRT::Element::LocationArray la(1);
  ele->LocationVector(discretization,la,false);
  // extract local values of the global vectors
  std::vector<double> myvalues(la[0].lm_.size());
  DRT::UTILS::ExtractMyValues(*state,myvalues,la[0].lm_);

  // split velocity and pressure
  for (int inode=0; inode<nen_; ++inode)  // number of nodes
  {
    // fill a vector field via a pointer
    for(int idim=0; idim<nsd_; ++idim) // number of dimensions
    {
      evelnp(idim,inode) = myvalues[idim+(inode*numdofpernode_)];
    }  // end for(idim)

    // fill a scalar field via a pointer
    eprenp(inode,0) = myvalues[nsd_+(inode*numdofpernode_)];
  }

  // velocity at int point
  //      _   _
  //     | u_0 |
  //     | u_1 |
  //     | u_2 |
  //     |_   _|
  //
  velint.Multiply(evelnp,shapefunct);

  // pressure at int point
  // scalar value
  // pseudo matrix "pressint" for ease of implementation
  //
  pressint.MultiplyTN(eprenp,shapefunct);

  //                                         _              _
  //                                        | u1,1 u1,2 u1,3 |
  // dudxi = u_i,alhpa = N_A,alpha u^A_i =  | u2,1 u2,2 u2,3 |
  //                                        |_u3,1 u3,2 u3,3_|
  //
  dudxi.MultiplyNT(evelnp,pderiv_loc);
  //                                            l=_  1     2     3  _
  //         -1                               i=1| u1,x1 u1,x2 u1,x3 |
  // dudxi o J  = N_A,alpha u^A_i xi_alpha,l =  2| u2,x1 u2,x2 u2,x3 | = gradu
  //                                            3|_u3,x1 u3,x2 u3,x3_|
  //
  dudxioJinv.MultiplyNT(dudxi,xji);

  // cauchystress = (gradu)^T
  cauchystress.UpdateT(1.0,dudxioJinv);
  // chauchystress = gradu + (gradu)^T
  cauchystress.Update(1.0,dudxioJinv,1.0);
  // cauchystress = tau
  cauchystress.Scale(fluiddynamicviscosity);
  // cauchystress finished
  cauchystress(0,0) += -pressint(0,0);
  cauchystress(1,1) += -pressint(0,0);
  cauchystress(2,2) += -pressint(0,0);

  // save cauchystress row for row as vector [11 22 33 12 23 13]
  elevec1_epetra(0)=cauchystress(0,0);
  elevec1_epetra(1)=cauchystress(1,1);
  elevec1_epetra(2)=cauchystress(2,2);
  elevec1_epetra(3)=cauchystress(0,1);
  elevec1_epetra(4)=cauchystress(1,2);
  elevec1_epetra(5)=cauchystress(0,2);

  if(elevec1_epetra.M()==7)
  {
    elevec1_epetra(6)=-pressint(0,0);
  }

    return 0;
}

/*-----------------------------------------------------------------------------*
 | Update 'IsImmersed' info in nodes and elements                rauch 05/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::UpdateImmersedInformation(
    Teuchos::ParameterList&              params,
    DRT::ELEMENTS::Fluid*                ele,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm
    )
{
  DRT::Problem* globalproblem = DRT::Problem::Instance();
  if(globalproblem->ProblemType() != prb_immersed_cell)
    dserror("UpdateImmersedInformation() not intended to be used with ProblemType other than prb_immersed_cell.\n"
        "    First check implementation before using this method in another context!");

  //-------------------------------------------------------------------------------
  //  This method provides the fluid discretization with the information about the
  //  overlapping immersed and background domains.  This is needed, when the mesh
  //  positions of at least one participating discretization has changed.
  //-------------------------------------------------------------------------------

  DRT::ELEMENTS::FluidImmersedBase* immersedele = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(ele);

  std::string backgrddisname(discretization.Name());
  std::string immerseddisname(params.get<std::string>("immerseddisname"));

  static double searchradiusfac = globalproblem->ImmersedMethodParams().get<double>("FLD_SRCHRADIUS_FAC");

  const Teuchos::RCP<DRT::Discretization> backgrddis  = globalproblem->GetDis(backgrddisname);
  const Teuchos::RCP<DRT::Discretization> immerseddis = globalproblem->GetDis(immerseddisname);


  // determine whether fluid mesh is deformable or not.
  // true for immersed cell migration. Not used in any other case so far.
  static int isALE = (globalproblem->ImmersedMethodParams().get<std::string>("DEFORM_BACKGROUND_MESH")=="yes");

  // initialize vectors for interpolation
  std::vector<double> dummy(4);


  std::vector<double> targeteledisp(nsd_*nen_);

  // update fluid displacements
  if(isALE)
  {
    LINALG::Matrix<nsd_,nen_>       edispnp(true);
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL, "dispnp");

    for(int node=0;node<nen_;++node)
      for(int dof=0; dof<nsd_;++dof)
        targeteledisp[node*nsd_+dof] = edispnp(dof,node);
  }
  else
  {
    // fluid target elements have no displacements in standard case. fill dummy vector.
    for(int i=0;i<nsd_*nen_;++i)
      targeteledisp[i]=0.0;
  }

  // parameter space coordinates of nodes 0 to 7 according to global report
  std::vector<std::vector<double> >nodalrefcoords(8);
  nodalrefcoords[0].push_back(-1.0); nodalrefcoords[0].push_back(-1.0); nodalrefcoords[0].push_back(-1.0);
  nodalrefcoords[1].push_back(1.0);  nodalrefcoords[1].push_back(-1.0); nodalrefcoords[1].push_back(-1.0);
  nodalrefcoords[2].push_back(1.0);  nodalrefcoords[2].push_back(1.0);  nodalrefcoords[2].push_back(-1.0);
  nodalrefcoords[3].push_back(-1.0); nodalrefcoords[3].push_back(1.0);  nodalrefcoords[3].push_back(-1.0);
  nodalrefcoords[4].push_back(-1.0); nodalrefcoords[4].push_back(-1.0); nodalrefcoords[4].push_back(1.0);
  nodalrefcoords[5].push_back(1.0);  nodalrefcoords[5].push_back(-1.0); nodalrefcoords[5].push_back(1.0);
  nodalrefcoords[6].push_back(1.0);  nodalrefcoords[6].push_back(1.0);  nodalrefcoords[6].push_back(1.0);
  nodalrefcoords[7].push_back(-1.0); nodalrefcoords[7].push_back(1.0);  nodalrefcoords[7].push_back(1.0);

  // get immersed structure search tree
  Teuchos::RCP<GEO::SearchTree> struct_searchtree = params.get<Teuchos::RCP<GEO::SearchTree> >("structsearchtree_rcp");

  // search tree related stuff
  std::map<int,LINALG::Matrix<3,1> >* currpositions_struct = params.get<std::map<int,LINALG::Matrix<3,1> >* >("currpositions_struct");

  // subset of strucutral elements immersed near the current fluid element
  std::map<int,std::set<int> > curr_subset_of_structdis;

  {
    // search radius (diagonal length of targetele (current fluid ele) )
    double radius = 0.0;
    LINALG::Matrix<3,1> searchcenter; // center of fluid ele
    if(isALE)
    {
      LINALG::Matrix<nsd_,nen_>       edispnp(true);
      ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"dispnp");

      radius = sqrt(pow((ele->Nodes()[1]->X()[0]+edispnp(0,1))-(ele->Nodes()[7]->X()[0]+edispnp(0,7)),2)+pow((ele->Nodes()[1]->X()[1]+edispnp(1,1))-(ele->Nodes()[7]->X()[1]+edispnp(1,7)),2)+pow((ele->Nodes()[1]->X()[2]+edispnp(2,1))-(ele->Nodes()[7]->X()[2]+edispnp(2,7)),2));
      searchcenter(0)=(ele->Nodes()[1]->X()[0]+edispnp(0,1))+((ele->Nodes()[7]->X()[0]+edispnp(0,7))- (ele->Nodes()[1]->X()[0]+edispnp(0,1)))*0.5;
      searchcenter(1)=(ele->Nodes()[1]->X()[1]+edispnp(1,1))+((ele->Nodes()[7]->X()[1]+edispnp(1,7))- (ele->Nodes()[1]->X()[1]+edispnp(1,1)))*0.5;
      searchcenter(2)=(ele->Nodes()[1]->X()[2]+edispnp(2,1))+((ele->Nodes()[7]->X()[2]+edispnp(2,7))- (ele->Nodes()[1]->X()[2]+edispnp(2,1)))*0.5;
    }
    else
    {
      radius = sqrt(pow(ele->Nodes()[1]->X()[0]-ele->Nodes()[7]->X()[0],2)+pow(ele->Nodes()[1]->X()[1]-ele->Nodes()[7]->X()[1],2)+pow(ele->Nodes()[1]->X()[2]-ele->Nodes()[7]->X()[2],2));
      searchcenter(0)=ele->Nodes()[1]->X()[0]+(ele->Nodes()[7]->X()[0] - ele->Nodes()[1]->X()[0])*0.5;
      searchcenter(1)=ele->Nodes()[1]->X()[1]+(ele->Nodes()[7]->X()[1] - ele->Nodes()[1]->X()[1])*0.5;
      searchcenter(2)=ele->Nodes()[1]->X()[2]+(ele->Nodes()[7]->X()[2] - ele->Nodes()[1]->X()[2])*0.5;
    }
    // search for immersed elements within a certain radius around the searchcenter node
    curr_subset_of_structdis = struct_searchtree->searchElementsInRadius(*immerseddis,*currpositions_struct,searchcenter,radius*searchradiusfac,0);
  }

  bool match = false;
  int  matchnum = 0;

  if(curr_subset_of_structdis.size()>0)
  {
    for (int node=0; node<nen_; node++)
    {
      std::vector<double> backgrdxi(nsd_);
      backgrdxi[0] = nodalrefcoords[node][0];
      backgrdxi[1] = nodalrefcoords[node][1];
      backgrdxi[2] = nodalrefcoords[node][2];

      if(static_cast<IMMERSED::ImmersedNode* >(ele->Nodes()[node])->IsMatched())
      {
        match = true;
      }

      IMMERSED::InterpolateToBackgrdPoint
      <DRT::Element::hex8,          // source/structure
      DRT::Element::hex8>           // target/fluid
      (curr_subset_of_structdis,
          immerseddis,              // source/structure
          backgrddis,               // target/fluid
          *ele,
          backgrdxi,
          targeteledisp,
          "none",
          dummy,                   // result
          match,
          false,
          false                    // do no communication. immerseddis is ghosted. every proc finds an immersed element
      );                           // to interpolate to its backgrd nodes.

      if(match)
      {
        matchnum++;
        static_cast<IMMERSED::ImmersedNode* >(ele->Nodes()[node])->SetIsMatched(1);

      } // if match

      // reset match to false and check next node in the following loop execution
      match = false;

    } // loop over all nodes of this element
  } // if immersed elements are in vicinity of ele

  // set ele "IsImmersed" if all nodes lie underneath the immersed dis (i.e. matched = true)
  if(matchnum==nen_)
    immersedele -> SetIsImmersed(1);
  // set ele "IsBoundaryImmersed" if 1<=x<8 nodes lie underneath the immersed dis ("HasProjectedDirichlet" is conjunction of "IsImmersed" and "IsBoundaryImmersed")
  else if (matchnum < nen_ and matchnum > 0)
  {
    immersedele -> SetBoundaryIsImmersed(1);

    // loop over nodes of this ele and set IsBoundaryImmersed
    for (int node=0; node<nen_; node++)
      static_cast<IMMERSED::ImmersedNode* >(ele->Nodes()[node])->SetIsBoundaryImmersed(1);

  }

  return 0;
}

/*-----------------------------------------------------------------------------*
 | Interpolate velocity                                          rauch 05/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::InterpolateVelocityToNode(
    Teuchos::ParameterList&              params,
    DRT::ELEMENTS::Fluid*                ele,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm,
    Epetra_SerialDenseVector&            elevec1_epetra, // vectofill
    Epetra_SerialDenseVector&            elevec2_epetra  // given point in parameter space coordinates
    )
{
  //----------------------------------------------------------------------
  //  two major things are done here
  //  1) interpolation of structural velocity to fluid nodes covered by immersed structure
  //  2) interpolation of structural divergence to fluid integration points covered by immersed structure
  //----------------------------------------------------------------------

  DRT::ELEMENTS::FluidImmersedBase* immersedele = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(ele);

  DRT::Problem* globalproblem = DRT::Problem::Instance();

  // check if fluid interacton is switched ON
  // if NOT : just mark isimmersed and isboundaryimmersed elements
  int isfluidinteraction=1;
  if(globalproblem->ProblemType()==prb_immersed_cell)
    isfluidinteraction = (globalproblem->CellMigrationParams().get<std::string>("FLUID_INTERACTION") == "yes");

  std::string backgrddisname(discretization.Name());
  std::string immerseddisname(params.get<std::string>("immerseddisname"));

  static double searchradiusfac = globalproblem->ImmersedMethodParams().get<double>("FLD_SRCHRADIUS_FAC");

  const Teuchos::RCP<DRT::Discretization> backgrddis  = globalproblem->GetDis(backgrddisname);
  const Teuchos::RCP<DRT::Discretization> immerseddis = globalproblem->GetDis(immerseddisname);

  // numgp in cut boundary elements
  static int num_gp_fluid_bound = globalproblem->ImmersedMethodParams().get<int>("NUM_GP_FLUID_BOUND");
  // degree of gp in cut boundary elements
  int degree_gp_fluid_bound = params.get("intpoints_fluid_bound",0);

  // determine whether fluid mesh is deformable or not
  static int isALE =
      (globalproblem->ImmersedMethodParams().get<std::string>("DEFORM_BACKGROUND_MESH")=="yes");

  // initialize vectors for interpolation
  std::vector<double> vel(numdofpernode_); // dofs 0,1,2
  std::vector<double> div(numdofpernode_); // dof  3


  std::vector<double> targeteledisp(nsd_*nen_,0.0);

  // update fluid displacements
  if(isALE)
  {
    LINALG::Matrix<nsd_,nen_>       edispnp(true);
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"dispnp");

    for(int node=0;node<nen_;++node)
      for(int dof=0; dof<nsd_;++dof)
        targeteledisp[node*nsd_+dof] = edispnp(dof,node);
  }

  const std::string action="interpolate_velocity_to_given_point";

  // parameter space coordinates of nodes 0 to 7 according to global report
  std::vector<std::vector<double> >nodalrefcoords(8);
  nodalrefcoords[0].push_back(-1.0); nodalrefcoords[0].push_back(-1.0); nodalrefcoords[0].push_back(-1.0);
  nodalrefcoords[1].push_back(1.0);  nodalrefcoords[1].push_back(-1.0); nodalrefcoords[1].push_back(-1.0);
  nodalrefcoords[2].push_back(1.0);  nodalrefcoords[2].push_back(1.0);  nodalrefcoords[2].push_back(-1.0);
  nodalrefcoords[3].push_back(-1.0); nodalrefcoords[3].push_back(1.0);  nodalrefcoords[3].push_back(-1.0);
  nodalrefcoords[4].push_back(-1.0); nodalrefcoords[4].push_back(-1.0); nodalrefcoords[4].push_back(1.0);
  nodalrefcoords[5].push_back(1.0);  nodalrefcoords[5].push_back(-1.0); nodalrefcoords[5].push_back(1.0);
  nodalrefcoords[6].push_back(1.0);  nodalrefcoords[6].push_back(1.0);  nodalrefcoords[6].push_back(1.0);
  nodalrefcoords[7].push_back(-1.0); nodalrefcoords[7].push_back(1.0);  nodalrefcoords[7].push_back(1.0);

  // get immersed structure search tree
  Teuchos::RCP<GEO::SearchTree> struct_searchtree = params.get<Teuchos::RCP<GEO::SearchTree> >("structsearchtree_rcp");

  // search tree related stuff
  std::map<int,LINALG::Matrix<3,1> >* currpositions_struct = params.get<std::map<int,LINALG::Matrix<3,1> >* >("currpositions_struct");

  // subset of strucutral elements immersed near the current fluid element
  std::map<int,std::set<int> > curr_subset_of_structdis;

  {
    //TEUCHOS_FUNC_TIME_MONITOR("DRT::ELEMENTS::InterpolateVelocityToNode() - search in searchtree");
    // search radius (diagonal length of targetele (current fluid ele) )
    double radius = 0.0;
    LINALG::Matrix<3,1> searchcenter; // center of fluid ele
    if(isALE)
    {
      LINALG::Matrix<nsd_,nen_>       edispnp(true);
      ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"dispnp");

      radius = sqrt(pow((ele->Nodes()[1]->X()[0]+edispnp(0,1))-(ele->Nodes()[7]->X()[0]+edispnp(0,7)),2)+pow((ele->Nodes()[1]->X()[1]+edispnp(1,1))-(ele->Nodes()[7]->X()[1]+edispnp(1,7)),2)+pow((ele->Nodes()[1]->X()[2]+edispnp(2,1))-(ele->Nodes()[7]->X()[2]+edispnp(2,7)),2));
      searchcenter(0)=(ele->Nodes()[1]->X()[0]+edispnp(0,1))+((ele->Nodes()[7]->X()[0]+edispnp(0,7))- (ele->Nodes()[1]->X()[0]+edispnp(0,1)))*0.5;
      searchcenter(1)=(ele->Nodes()[1]->X()[1]+edispnp(1,1))+((ele->Nodes()[7]->X()[1]+edispnp(1,7))- (ele->Nodes()[1]->X()[1]+edispnp(1,1)))*0.5;
      searchcenter(2)=(ele->Nodes()[1]->X()[2]+edispnp(2,1))+((ele->Nodes()[7]->X()[2]+edispnp(2,7))- (ele->Nodes()[1]->X()[2]+edispnp(2,1)))*0.5;
    }
    else
    {
      radius = sqrt(pow(ele->Nodes()[1]->X()[0]-ele->Nodes()[7]->X()[0],2)+pow(ele->Nodes()[1]->X()[1]-ele->Nodes()[7]->X()[1],2)+pow(ele->Nodes()[1]->X()[2]-ele->Nodes()[7]->X()[2],2));
      searchcenter(0)=ele->Nodes()[1]->X()[0]+(ele->Nodes()[7]->X()[0] - ele->Nodes()[1]->X()[0])*0.5;
      searchcenter(1)=ele->Nodes()[1]->X()[1]+(ele->Nodes()[7]->X()[1] - ele->Nodes()[1]->X()[1])*0.5;
      searchcenter(2)=ele->Nodes()[1]->X()[2]+(ele->Nodes()[7]->X()[2] - ele->Nodes()[1]->X()[2])*0.5;
    }
    // search for immersed elements within a certain radius around the searchcenter node
    curr_subset_of_structdis = struct_searchtree->searchElementsInRadius(*immerseddis,*currpositions_struct,searchcenter,radius*searchradiusfac,0);
  }

  bool match = false;
  int  matchnum = 0;

  /********************************************************************************/
  // 1) Interpolation of structural velocity
  //    (loop over all nodes of this element)
  /********************************************************************************/

  // in first step only velocity is interpolated
  bool vel_calculation = true;

  if(curr_subset_of_structdis.size()>0)
  {
    for (int node=0; node<nen_; node++)
    {
      std::vector<double> backgrdxi(nsd_);
      backgrdxi[0] = nodalrefcoords[node][0];
      backgrdxi[1] = nodalrefcoords[node][1];
      backgrdxi[2] = nodalrefcoords[node][2];

      if(static_cast<IMMERSED::ImmersedNode* >(ele->Nodes()[node])->IsMatched())
      {
        match = true;
      }

      IMMERSED::InterpolateToBackgrdPoint
      <DRT::Element::hex8,          // source/structure
      DRT::Element::hex8>           // target/fluid
      (curr_subset_of_structdis,
          immerseddis,              // source/structure
          backgrddis,               // target/fluid
          *ele,
          backgrdxi,
          targeteledisp,
          action,
          vel,                      // result
          match,
          vel_calculation,
          false                     // do no communication. immerseddis is ghosted. every proc finds an immersed element
      );                            // to interpolate to its backgrd nodes.

      // under fsi structure NOT under immersed structure !
      if(vel[0]<-12344.0 and vel[1]<-12344.0)
        match=false;

      if(match)
      {
        matchnum++;
        static_cast<IMMERSED::ImmersedNode* >(ele->Nodes()[node])->SetIsMatched(1);
        immersedele->SetHasProjectedDirichlet(1);

        for(int i=0; i<nsd_;++i)
        {
          elevec1_epetra((node*numdofpernode_)+i) += vel[i];
        }
      } // if match

      // reset match to false and check next node in the following loop execution
      match = false;

    } // loop over all nodes of this element
  } // if immersed elements are in vicinity of ele

  // set ele "IsImmersed" if all nodes lie underneath the immersed dis (i.e. matched = true)
  if(matchnum==nen_)
    immersedele -> SetIsImmersed(1);
  // set ele "IsBoundaryImmersed" if 1<=x<8 nodes lie underneath the immersed dis ("HasProjectedDirichlet" is conjunction of "IsImmersed" and "IsBoundaryImmersed")
  else if (matchnum < nen_ and matchnum > 0)
  {
    // inform background ele about immersed boundary
    immersedele -> SetBoundaryIsImmersed(1);

    // loop over nodes of this ele and set IsBoundaryImmersed
    for (int node=0; node<nen_; node++)
      static_cast<IMMERSED::ImmersedNode* >(ele->Nodes()[node])->SetIsBoundaryImmersed(1);

    if (isfluidinteraction)
    {
      immersedele -> ConstructElementRCP(num_gp_fluid_bound);

    // DEBUG test
#ifdef DEBUG
    if(immersedele->GetRCPProjectedIntPointDivergence()==Teuchos::null)
      dserror("construction of ProjectedIntPointDivergence failed");
    if((int)(immersedele->GetRCPProjectedIntPointDivergence()->size())!=num_gp_fluid_bound)
      dserror("size of ProjectedIntPointDivergence should be equal numgp in cut element = %d",num_gp_fluid_bound);

    // DEBUG test
    if(immersedele->GetRCPIntPointHasProjectedDivergence()==Teuchos::null)
      dserror("construction of IntPointHasProjectedDivergence failed");
    if((int)(immersedele->GetRCPIntPointHasProjectedDivergence()->size())!=num_gp_fluid_bound)
      dserror("size of IntPointHasProjectedDivergence should be equal numgp in cut element = %d",num_gp_fluid_bound);
# endif
    } // these vectors only need to be constructed when fluid interaction is switched ON
  }

  // interpolate divergence of immerseddis velocity to backgrd. int. points only if fluid interaction is switched ON
  if (isfluidinteraction)
  {
    /********************************************************************************/
    // 2) Interpolation of structural divergence
    //    (loop over all int points of elements set as "BoundaryIsImmersed")
    /********************************************************************************/

    // only velocity divergence needs to be calculated and interpolated here
    vel_calculation = false;
    // get integration rule of fluid element
    const DRT::UTILS::GaussIntegration intpoints_fluid_bound(distype,degree_gp_fluid_bound);

    if (degree_gp_fluid_bound)
    {
      if (immersedele->IsBoundaryImmersed())
      {
        for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints_fluid_bound.begin(); iquad!=intpoints_fluid_bound.end(); ++iquad )
        {
          std::vector<double> backgrdxi(nsd_);
          backgrdxi[0] = iquad.Point()[0];
          backgrdxi[1] = iquad.Point()[1];
          backgrdxi[2] = iquad.Point()[2];

          bool gp_has_projected_divergence = false;
          if(immersedele->GetRCPIntPointHasProjectedDivergence() != Teuchos::null)
            if(immersedele->GetRCPIntPointHasProjectedDivergence()->size() > 0)
              gp_has_projected_divergence = (int)immersedele->IntPointHasProjectedDivergence(*iquad);

          if(gp_has_projected_divergence)
          {
            match = true;
          }

          IMMERSED::InterpolateToBackgrdPoint
          <DRT::Element::hex8,          // source/structure
          DRT::Element::hex8>           // target/fluid
          (curr_subset_of_structdis,
              immerseddis,              // source/structure
              backgrddis,               // target/fluid
              *ele,
              backgrdxi,
              targeteledisp,
              action,
              div,                      // result (in dof 3)
              match,
              vel_calculation,
              false                     // do no communication. immerseddis is ghosted. every proc finds an immersed element
          );                            // to interpolate to its backgrd nodes.

          // under fsi structure NOT under immersed structure !
          if(div[nsd_]<-12344.0)
            match=false;

          if(match)
          {
            immersedele->SetIntPointHasProjectedDivergence(*iquad,1);
            immersedele->StoreProjectedIntPointDivergence(*iquad,div[nsd_]);
          }

          // reset match to false and check next int point in the following loop execution
          match = false;

        } // loop over int points
      } // only if IsBoundaryImmersed
    } // degree_gp_fluid_bound > 0
    else
      dserror("In case of fluid interaction a proper value for NUM_GP_FLUID_BOUND must be set in your .dat file.\n"
              "(valid parameters are 8, 64, 125, 343, 729 and 1000). Check also if you forgot to set the value in \n"
              "the parameter list provided for this action.");
  }
  return 0;
}

/*-----------------------------------------------------------------------------*
 | search for immersed boundary fluid elements                  sfuchs 11/2016 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::SearchImmersedBoundaryElements(
    Teuchos::ParameterList&              params,
    DRT::ELEMENTS::Fluid*                ele,
    Epetra_SerialDenseMatrix&            elemat1,              // dummy matrix
    Epetra_SerialDenseVector&            elevec1               // dummy vector
    )
{
  // throw an error if no hex8 fluid distype
  if(distype != DRT::Element::hex8)
    dserror("immersed partitioned fsi with membrane finite element currently just implemented for hex8 fluid discretisation type!");

  DRT::Problem* globalproblem = DRT::Problem::Instance();

  // throw an error if problem type is not provided
  if(globalproblem->ProblemType()!=prb_immersed_membrane_fsi)
    dserror("Function SearchImmersedBoundaryElements() currently just implemented for problem type 'prb_immersed_membrane_fsi'!");

  // get structure discretization
  std::string structdisname(params.get<std::string>("immerseddisname"));
  const Teuchos::RCP<DRT::Discretization> structdis = globalproblem->GetDis(structdisname);

  // setup parameters for search of (potentially) immersed structural elements
  static double searchradiusfac = globalproblem->ImmersedMethodParams().get<double>("FLD_SRCHRADIUS_FAC");

  // get immersed structure search tree
  Teuchos::RCP<GEO::SearchTree> struct_searchtree = params.get<Teuchos::RCP<GEO::SearchTree> >("structsearchtree_rcp");

  // search tree related stuff
  std::map<int,LINALG::Matrix<3,1> >* currpositions_struct = params.get<std::map<int,LINALG::Matrix<3,1> >* >("currpositions_struct");

  // subset of structural elements immersed near the current fluid element
  std::map<int,std::set<int> > curr_subset_of_structdis;

  // search radius (diagonal length of current fluid ele)
  double radius = 0.0;
  // center of fluid ele
  LINALG::Matrix<3,1> searchcenter;

  radius = sqrt(pow(ele->Nodes()[1]->X()[0]-ele->Nodes()[7]->X()[0],2)+pow(ele->Nodes()[1]->X()[1]-ele->Nodes()[7]->X()[1],2)+pow(ele->Nodes()[1]->X()[2]-ele->Nodes()[7]->X()[2],2));
  searchcenter(0)=ele->Nodes()[1]->X()[0]+(ele->Nodes()[7]->X()[0] - ele->Nodes()[1]->X()[0])*0.5;
  searchcenter(1)=ele->Nodes()[1]->X()[1]+(ele->Nodes()[7]->X()[1] - ele->Nodes()[1]->X()[1])*0.5;
  searchcenter(2)=ele->Nodes()[1]->X()[2]+(ele->Nodes()[7]->X()[2] - ele->Nodes()[1]->X()[2])*0.5;

  // search for immersed elements within a certain radius around the searchcenter node
  curr_subset_of_structdis = struct_searchtree->searchElementsInRadius(*structdis,*currpositions_struct,searchcenter,radius*searchradiusfac,0);

  // check if a structural node or gauss point lies inside the current fluid element
  // for immersed boundary fluid elements set the flags SetIsImmersed(1) and SetHasProjectedDirichlet(1)
  IMMERSED::GetInFluidEleImmersedPoints<DRT::Element::quad4,    // structure
                                        DRT::Element::hex8>     // fluid
      (curr_subset_of_structdis,   // relevant structural elements
       structdis,                  // structure discretization
       *ele,                       // current fluid element
       elemat1,                    // element least squares matrix
       elevec1,                    // element least squares vector
       false                       // true: fill least squares matrix and vector, false: just check if ele is immersed
      );

  // Remark:
  // SetBoundaryIsImmersed(1) is not set, as function call ProjectedIntPointDivergence()
  // in DRT::ELEMENTS::FluidEleCalcImmersed<distype>::ContinuityGalPart() is not intended
  //
  // the flag SetIsMatched(1) for corresponding immersed fluid nodes is set in
  // IMMERSED::ImmersedPartitionedFSIDirichletNeumannMembrane::CalcArtificialVelocity()
  // once the nodal dofs are added to the DofRowMap of immersed fluid nodes

  return 0;
}

/*-----------------------------------------------------------------------------*
 | get least squares matrix and rhs for immersed boundary ele   sfuchs 11/2016 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::GetLeastSquaresMatrixImmersedBoundary(
    Teuchos::ParameterList&              params,
    DRT::ELEMENTS::Fluid*                ele,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm,
    Epetra_SerialDenseMatrix&            elemat1,              // element least squares matrix
    Epetra_SerialDenseVector&            elevec1               // element least squares vector
    )
{
  // throw an error if no hex8 fluid distype
  if(distype != DRT::Element::hex8)
    dserror("immersed partitioned fsi with membrane finite element currently just implemented for hex8 fluid discretisation type!");

  DRT::Problem* globalproblem = DRT::Problem::Instance();

  // throw an error if problem type is not provided
  if(globalproblem->ProblemType()!=prb_immersed_membrane_fsi)
    dserror("Function GetLeastSquaresMatrixImmersedBoundary() currently just implemented for problem type 'prb_immersed_membrane_fsi'!");

  // get structure discretization
  std::string structdisname(params.get<std::string>("immerseddisname"));
  const Teuchos::RCP<DRT::Discretization> structdis = globalproblem->GetDis(structdisname);

  // setup parameters for search of (potentially) immersed structural elements
  static double searchradiusfac_struct = globalproblem->ImmersedMethodParams().get<double>("FLD_SRCHRADIUS_FAC");

  // get immersed structure search tree
  Teuchos::RCP<GEO::SearchTree> struct_searchtree = params.get<Teuchos::RCP<GEO::SearchTree> >("structsearchtree_rcp");

  // search tree related stuff
  std::map<int,LINALG::Matrix<3,1> >* currpositions_struct = params.get<std::map<int,LINALG::Matrix<3,1> >* >("currpositions_struct");

  // subset of structural elements immersed near the current fluid element
  std::map<int,std::set<int> > curr_subset_of_structdis;

  // search radius (diagonal length of targetele (current fluid ele) )
  double radius = 0.0;
  LINALG::Matrix<3,1> searchcenter; // center of fluid ele

  radius = sqrt(pow(ele->Nodes()[1]->X()[0]-ele->Nodes()[7]->X()[0],2)+pow(ele->Nodes()[1]->X()[1]-ele->Nodes()[7]->X()[1],2)+pow(ele->Nodes()[1]->X()[2]-ele->Nodes()[7]->X()[2],2));
  searchcenter(0)=ele->Nodes()[1]->X()[0]+(ele->Nodes()[7]->X()[0] - ele->Nodes()[1]->X()[0])*0.5;
  searchcenter(1)=ele->Nodes()[1]->X()[1]+(ele->Nodes()[7]->X()[1] - ele->Nodes()[1]->X()[1])*0.5;
  searchcenter(2)=ele->Nodes()[1]->X()[2]+(ele->Nodes()[7]->X()[2] - ele->Nodes()[1]->X()[2])*0.5;

  // search for immersed elements within a certain radius around the searchcenter node
  curr_subset_of_structdis = struct_searchtree->searchElementsInRadius(*structdis,*currpositions_struct,searchcenter,radius*searchradiusfac_struct,0);

  // fill entries of element least squares matrix and vector containing structural point position and velocity
  IMMERSED::GetInFluidEleImmersedPoints<DRT::Element::quad4,    // structure
                                        DRT::Element::hex8>     // fluid
      (curr_subset_of_structdis,   // relevant structural elements
       structdis,                  // structure discretization
       *ele,                       // current fluid element
       elemat1,                    // element least squares matrix
       elevec1,                    // element least squares vector
       true                        // true: fill least squares matrix and vector, false: just check if ele is immersed
      );

  // get fluid nodal parameter space coordinates
  Epetra_SerialDenseMatrix locations = DRT::UTILS::getEleNodeNumbering_nodes_paramspace(distype);

  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // double loop over fluid nodes
  for (int k=0; k<nen_; ++k)
  {
    // evaluate shape functions and derivatives at fluid node k
    EvalShapeFuncAndDerivsAtIntPoint(locations[k],1.0);

    for (int l=0; l<nen_; ++l)
    {
      // loop over space dimension
      for (int dim=0; dim<nsd_; ++dim)
      {
        // fill entries of element least squares matrix containing shape function derivatives
        elemat1( (numdofpernode_*k)+3 , (numdofpernode_*l)+dim ) += derxy_(dim,l);
        elemat1( (numdofpernode_*l)+dim , (numdofpernode_*k)+3 ) += derxy_(dim,l);
      }
    }
  }

  return 0;
}

/*-----------------------------------------------------------------------------*
 | Correct Immersed Boundary Velocities                          rauch 07/2015 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CorrectImmersedBoundVelocities(
    Teuchos::ParameterList&              params,
    DRT::ELEMENTS::Fluid*                ele,
    DRT::Discretization &                discretization, // fluid
    const std::vector<int> &             lm,
    Teuchos::RCP<MAT::Material> &        mat,
    Epetra_SerialDenseVector&            elevec1_epetra, // vectofill
    Epetra_SerialDenseVector&            elevec2_epetra
    )
{
  // cast fluid element to immersed element to get/store immersed information
  DRT::ELEMENTS::FluidImmersedBase* immersedele = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(ele);

  // do only if immersed ele is cut by boundary
  if(immersedele->IsBoundaryImmersed())
  {
  // get global problem
  DRT::Problem* globalproblem = DRT::Problem::Instance();

  // get factor for search radius
  static double searchradiusfac = globalproblem->ImmersedMethodParams().get<double>("FLD_SRCHRADIUS_FAC");

  // get discretizations
  Teuchos::RCP<DRT::Discretization> fluid_dis  = Teuchos::null;
  Teuchos::RCP<DRT::Discretization> struct_dis = Teuchos::null;

  if(globalproblem->ProblemType()==prb_immersed_cell)
  {
    // get discretizations
    fluid_dis  = globalproblem->GetDis("porofluid");
    struct_dis = globalproblem->GetDis("cell");
  }
  else
  {
    // get discretizations
    fluid_dis  = globalproblem->GetDis("fluid");
    struct_dis = globalproblem->GetDis("structure");
  }

  // determine whether fluid mesh is deformable or not
  static int isALE =
      (globalproblem->ImmersedMethodParams().get<std::string>("DEFORM_BACKGROUND_MESH")=="yes");

  // get element velocity at time n+1
  LINALG::Matrix<nsd_,nen_> evelnp;
  Teuchos::RCP<const Epetra_Vector> state;
  std::vector<double> myvalues(1);
  state = discretization.GetState("velnp");

  // fill locationarray
  DRT::Element::LocationArray la(1);
  ele->LocationVector(discretization,la,false);

  // extract local values of the global vectors
  myvalues.resize(la[0].lm_.size());
  DRT::UTILS::ExtractMyValues(*state,myvalues,la[0].lm_);

  // split velocity and pressure
  for (int inode=0; inode<nen_; ++inode)
  {
    for(int idim=0; idim<nsd_; ++idim)
    {
      evelnp(idim,inode) = myvalues[idim+(inode*numdofpernode_)];
    }
  }

  // initialize local position and velocity of closest point on structural surface that needs to be found
  std::vector<double> closest_point_xi(nsd_);
  std::vector<double> vel(nsd_);


  std::vector<double> targeteledisp(nsd_*nen_);

  // update fluid displacements
  if(isALE)
  {
    LINALG::Matrix<nsd_,nen_>       edispnp(true);
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"dispnp");

    for(int node=0;node<nen_;++node)
      for(int dof=0; dof<nsd_;++dof)
        targeteledisp[node*nsd_+dof] = edispnp(dof,node);
  }
  else
  {
    for(int i=0;i<nsd_*nen_;++i)
      targeteledisp[i]=0.0;
  }

  // set action for structure evaluation
  const std::string action="interpolate_velocity_to_given_point";

  // parameter space coordinates of nodes 0 to 7 according to global report
  std::vector<std::vector<double> >nodalrefcoords(8);
  nodalrefcoords[0].push_back(-1.0); nodalrefcoords[0].push_back(-1.0); nodalrefcoords[0].push_back(-1.0);
  nodalrefcoords[1].push_back(1.0);  nodalrefcoords[1].push_back(-1.0); nodalrefcoords[1].push_back(-1.0);
  nodalrefcoords[2].push_back(1.0);  nodalrefcoords[2].push_back(1.0);  nodalrefcoords[2].push_back(-1.0);
  nodalrefcoords[3].push_back(-1.0); nodalrefcoords[3].push_back(1.0);  nodalrefcoords[3].push_back(-1.0);
  nodalrefcoords[4].push_back(-1.0); nodalrefcoords[4].push_back(-1.0); nodalrefcoords[4].push_back(1.0);
  nodalrefcoords[5].push_back(1.0);  nodalrefcoords[5].push_back(-1.0); nodalrefcoords[5].push_back(1.0);
  nodalrefcoords[6].push_back(1.0);  nodalrefcoords[6].push_back(1.0);  nodalrefcoords[6].push_back(1.0);
  nodalrefcoords[7].push_back(-1.0); nodalrefcoords[7].push_back(1.0);  nodalrefcoords[7].push_back(1.0);

  // get structure search tree
  Teuchos::RCP<GEO::SearchTree> struct_searchtree = params.get<Teuchos::RCP<GEO::SearchTree> >("structsearchtree_rcp");

  // search tree related stuff
  std::map<int,LINALG::Matrix<3,1> >* currpositions_struct = params.get<std::map<int,LINALG::Matrix<3,1> >* >("currpositions_struct");

  // get relevant structure elements
  std::map<int,std::set<int> > curr_subset_of_structdis;
  {
    // search radius (diagonal length of fluid element)
    double radius = sqrt(pow(ele->Nodes()[1]->X()[0]-ele->Nodes()[7]->X()[0],2)+pow(ele->Nodes()[1]->X()[1]-ele->Nodes()[7]->X()[1],2)+pow(ele->Nodes()[1]->X()[2]-ele->Nodes()[7]->X()[2],2));
    LINALG::Matrix<3,1> searchcenter; // center of fluid element
    searchcenter(0)=ele->Nodes()[1]->X()[0]+(ele->Nodes()[7]->X()[0] - ele->Nodes()[1]->X()[0])*0.5;
    searchcenter(1)=ele->Nodes()[1]->X()[1]+(ele->Nodes()[7]->X()[1] - ele->Nodes()[1]->X()[1])*0.5;
    searchcenter(2)=ele->Nodes()[1]->X()[2]+(ele->Nodes()[7]->X()[2] - ele->Nodes()[1]->X()[2])*0.5;

    // search for immersed elements within a certain radius around the search center node
    curr_subset_of_structdis = struct_searchtree->searchElementsInRadius(*struct_dis,*currpositions_struct,searchcenter,radius*searchradiusfac,0);
  }

  //*********************************
  // loop over all nodes
  //*********************************

  bool match = false;
  for (int node=0; node<this->nen_; node++)
    {
      std::vector<double> backgrdfluidxi(nsd_);
      backgrdfluidxi[0] = nodalrefcoords[node][0];
      backgrdfluidxi[1] = nodalrefcoords[node][1];
      backgrdfluidxi[2] = nodalrefcoords[node][2];

    if(static_cast<IMMERSED::ImmersedNode* >(ele->Nodes()[node])->IsMatched())
    {
      match = true;
    }

    IMMERSED::FindClosestStructureSurfacePoint  <DRT::Element::quad4,                 // structure
                                                 DRT::Element::hex8>                  // fluid
                                                         (curr_subset_of_structdis,   // relevant struct elements
                                                          struct_dis,                 // structure discretization
                                                          fluid_dis,                  // fluid discretization
                                                          *ele,                       // fluid element
                                                          backgrdfluidxi,             // space coordinate of node
                                                          targeteledisp,              // fluid displacements (zero)
                                                          action,                     // action for structure evaluation
                                                          vel,                        // velocity result
                                                          closest_point_xi,           // xi position of closest point
                                                          match,                      // found a closest point
                                                          false                       // do no communication. struct_dis is ghosted. every proc finds an immersed element to interpolate to its backgrd nodes
                                                          );

    // only if closest point lies in this element match=true, otherwise this node is matched by another element
    if(match)
    {
      // if closest point to node lying in this element is found, node is set matched to indicate that now has an dirichlet value
      static_cast<IMMERSED::ImmersedNode* >(ele->Nodes()[node])->SetIsMatched(1);
      // this is done before anyway, but doesn't hurt here
      immersedele->SetHasProjectedDirichlet(1);

      // evaluate shape function of respective node in the closest structure point (has to  lie in the current fluid element)
      LINALG::Matrix<nen_,1> shapefunct;
      double weight = 0.0;

      // get position of closest point in local coordinates of fluid element
      LINALG::Matrix<nsd_,1> xi;
      for(int i=0;i<nsd_;++i)
        xi(i) = closest_point_xi[i];

      // evaluate shape functions at closest point
      DRT::UTILS::shape_function<distype>(xi,shapefunct);
      weight = shapefunct(node,0);

      // calculate new node velocities by weighting the influence of Navier Stokes solution and interpolation via distance of closest point to this node
      for(int idim=0; idim<nsd_;++idim)
      {
        elevec1_epetra((node*numdofpernode_)+idim) = weight*vel[idim] + (1.0 - weight)*evelnp(idim,node);
      }
    } // end if match

    // reset match to false and check next node in the following loop execution
    match = false;

    } // end loop over all nodes
  } // if fluid element has immersed boundary
  return 0;

} // CorrectImmersedBoundVelocities()


/*---------------------------------------------------------------------*
 | Action type: interpolate_velocity_to_given_point                    |
 | calculate velocity at given point                       ghamm 12/15 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::InterpolateVelocityToPoint(
    DRT::ELEMENTS::Fluid*     ele,
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2)
{
  // coordinates of the current integration point
  LINALG::Matrix<nsd_,1> elecoords = params.get<LINALG::Matrix<nsd_,1> >("elecoords");

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
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"disp");

    // get new node positions for isale
     xyze_ += edispnp;
  }

  // the int point considered is the point given from outside
  EvalShapeFuncAndDerivsAtIntPoint(elecoords.A(), -1.0);

  //----------------------------------------------------------------------------
  //   Extract velocity from global vectors and compute velocity at point
  //----------------------------------------------------------------------------

  static LINALG::Matrix<nsd_,nen_> evel;
  // fill the local element vector with the global values
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evel, NULL,"vel");
  velint_.Multiply(evel,funct_);

  for (int isd=0;isd<nsd_;isd++)
  {
    elevec1[isd] = velint_(isd);
  }


  return 0;
}


/*---------------------------------------------------------------------*
 | Action type: interpolate_pressure_to_given_point                    |
 | calculate pressure at given point                       ghamm 06/15 |
 *---------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::InterpolatePressureToPoint(
    DRT::ELEMENTS::Fluid*     ele,
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1)
{
  // coordinates of the current integration point
  LINALG::Matrix<nsd_,1> elecoords = params.get<LINALG::Matrix<nsd_,1> >("elecoords");

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
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"disp");

    // get new node positions for isale
     xyze_ += edispnp;
  }

  // the int point considered is the point given from outside
  EvalShapeFuncAndDerivsAtIntPoint(elecoords.A(), -1.0);

  //----------------------------------------------------------------------------
  //   Extract pressure from global vectors and compute pressure at point
  //----------------------------------------------------------------------------

  static LINALG::Matrix<nen_,1> epre;

  if(discretization.HasState("vel"))
  {
    // fill the local element vector with the global values
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, NULL, &epre,"vel");
    elevec1[0] = funct_.Dot(epre);
  }

  if(discretization.HasState("velnp"))
  {
    // fill the local element vector with the global values
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, NULL, &epre,"velnp");

    if (elevec1.Length() != 2)
      dserror("velnp is set, there must be a vel as well");

    elevec1[1] = funct_.Dot(epre);
  }

  return 0;
}

/*-----------------------------------------------------------------------------*
 | Integrate velocity divergence                                 rauch 05/2014 |
 | Works only for Cell Migration Simulation (poro-immersed)                    |
 |    CALC:   /                                                                |
 |           /                                                                 |
 |           |                                                                 |
 |           | div(u)*phi dv                                                   |
 |           |                                                                 |
 |          /                                                                  |
 |         / d V_t                                                             |
 |                                                                             |
 | phi: porosity                                                               |
 | u  : fluid velocity                                                         |
 |                                                                             |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcArtificialVelocityDivergence(
    Teuchos::ParameterList&              params,
    DRT::ELEMENTS::Fluid*                ele,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm,
    Epetra_SerialDenseVector&            elevec1_epetra) // vectofill
    {

  // cast to immersed ele
  DRT::ELEMENTS::FluidImmersedBase* immersedele = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(ele);

  // only if all nodes are underneath the immersed dis (source/sink is fixed)
  if(immersedele->IsImmersed())
  {
    // set element id
    eid_ = ele->Id();

    // extract element porosities from structure discretization
    const Teuchos::RCP<DRT::Discretization> structuredis  = DRT::Problem::Instance()->GetDis("structure");
    const Teuchos::RCP<const Epetra_Vector> dispnp = structuredis->GetState("dispnp");
    if(dispnp == Teuchos::null)
      dserror("Could not get state 'dispnp' from structural discretization");
    std::vector<double> mydispnp(lm.size());
    DRT::Element* structele = structuredis->gElement(eid_);
    DRT::Element::LocationArray la(1);
    structele->LocationVector(*structuredis,la,false);

    DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,la[0].lm_);

    // extract fluid state
    LINALG::Matrix<nsd_,nen_> evelaf(true);
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelaf, NULL,"velnp");


    // get geometry and update if necessary

    // get node coordinates
    GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

//    if (ele->IsAle())
//    {
//      LINALG::Matrix<nsd_,nen_>       edispnp(true);
//      ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &edispnp, NULL,"dispnp");
//
//      // get new node positions for isale
//      xyze_ += edispnp;
//    }

    double porosity = 1234.0;
    static double rhsfac = DRT::ELEMENTS::FluidEleParameterTimInt::Instance()->TimeFacRhs();

    ///////////////////////////////////////////////////////////////////////////////////////
    //
    // Integrate source/sink for continuity equation over the element
    //
    ///////////////////////////////////////////////////////////////////////////////////////

    // loop over all integration points
    // ---------------------------------------------------------------------------
    // Integration loop
    // ---------------------------------------------------------------------------
    for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
    {
      // evaluate shape functions and derivatives at integration point
      EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

      porosity = 0.0;
      // evaluate porosity at current integration point
      for(int node=0; node<nen_;++node)
      {
        porosity += funct_(node) * mydispnp[node*4+3];
      }

      // get velocity at integration point
      // (values at n+alpha_F for generalized-alpha scheme, n+1 otherwise)
      vderxy_.MultiplyNT(evelaf,derxy_);

      vdiv_= 0.0;
      for (int idim = 0; idim <nsd_; ++idim)
      {
        vdiv_ += vderxy_(idim, idim);
      }

      for(int node=0;node<nen_;++node)
      {
        elevec1_epetra((node*numdofpernode_)+nsd_) += vdiv_*fac_*funct_(node)*porosity*rhsfac;
      }

    }

  }// if isimmersed

  return 0;
    }

/*-----------------------------------------------------------------------------*
 | Calculate channel statistics                                     bk 05/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcChannelStatistics(
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
    // get new node positions of ALE mesh
    GetGridDispALE(discretization, lm, edispnp);

     for (int inode=0; inode<nen_; inode++)
     {
       if(abs(edispnp(normdirect,inode))>1e-6)
       {
         dserror("no sampling possible if homogeneous planes are not conserved\n");
       }
     }
  }

  //we have to fill the virtual nodes of the xyze_ matrix, since we only get
  //coordinates of the real nodes
  if(enrtype == DRT::ELEMENTS::Fluid::xwall)
  {
    for (int inode=0; inode<nen_/2; inode++)
    {
      for(int sdm=0;sdm<nsd_;++sdm)
      {
        xyze_(sdm,inode+nen_/2)=xyze_(sdm,inode);//watch out, this is not in the correct order
        //but otherwise we would have to introduce a new copy which is costly
      }
    }
  }

  if(distype == DRT::Element::hex8
     ||
     distype == DRT::Element::hex27
     ||
     distype == DRT::Element::hex20)
  {
    //decide first, if this element is taken into account!
    double inflowmax = params.get<double>("INFLOW_CHA_SIDE");
    // get the minimum x coordinate of this element and compare with the maximum of the inflow channel
    double minx = 9998.0;
    for(int inode=0;inode<nen_;inode++)
    {
      if(minx > xyze_(0,inode))
      {
        minx = xyze_(0,inode);
      }
    }
    if(inflowmax < minx)
      return 0;


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
//    if(elenormdirect==0)
//    {
//      std::cout << "numplanes in ele:   " << numplanesinele<< std::endl;
//      for(std::set<int>::const_iterator id = planesinele.begin();id!=planesinele.end() ;++id)
//        std::cout <<"planesinele id:  "<< *id << std::endl;
//
//    }
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
        {
          const double* econst=e;

          // evaluate shape functions and derivatives at integration point
          EvalShapeFuncAndDerivsAtIntPoint(econst,iquad.Weight());
        }

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
        if(enrtype == DRT::ELEMENTS::Fluid::xwall)
        {
          if(inode<nen_/2)
          {//funct and xyze not in the same order
            x[0]+=funct_(inode*2)*xyze_(0,inode);
            x[1]+=funct_(inode*2)*xyze_(1,inode);
            x[2]+=funct_(inode*2)*xyze_(2,inode);
          }
        }
        else
        {
          x[0]+=funct_(inode)*xyze_(0,inode);
          x[1]+=funct_(inode)*xyze_(1,inode);
          x[2]+=funct_(inode)*xyze_(2,inode);
        }
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
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcTimeStep(
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
  // ---------------------------------------------------------------------
  // get additional state vectors for ALE case: grid displacement and vel.
  // ---------------------------------------------------------------------
  LINALG::Matrix<nsd_, nen_> edispnp(true);
  LINALG::Matrix<nsd_, nen_> egridv(true);
  if (ele->IsAle())
    GetGridDispVelALE(discretization, lm, edispnp, egridv);


  // evaluate shape functions and derivatives element center
  EvalShapeFuncAndDerivsAtEleCenter();

  // np_genalpha: additional vector for velocity at time n+1
  LINALG::Matrix<nsd_,nen_> evelnp(true);

  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelnp, NULL,"velnp");

  convvelint_.Multiply(evelnp,funct_);

  //calculate element length via the stream length definition, see corresponding implementation
  //in fluid_ele_calc for calculation of the stabilization parameter
  double h=0.0;
  double vel_norm=0.0;

  if(ele->IsAle())
  {
    gridvelint_.Multiply(egridv,funct_);
    convvelint_.Update(-1.0,gridvelint_,1.0);
  }

  vel_norm = convvelint_.Norm2();

  if(vel_norm>1.0e-6)
  {
    LINALG::Matrix<nsd_,1> velino(true);
    velino.Update(1.0/vel_norm,convvelint_);

    // get streamlength using the normed velocity at element centre
    LINALG::Matrix<nen_,1> tmp;
    //enriched dofs are not interpolatory with respect to geometry
    if(enrtype == DRT::ELEMENTS::Fluid::xwall)
    {
      LINALG::Matrix<nsd_,nen_> derxy_copy(derxy_);
      for(int inode=1;inode<nen_;inode+=2)
      {
        for (int idim=0; idim<nsd_;idim++)
          derxy_copy(idim,inode)=0.0;
      }
      tmp.MultiplyTN(derxy_copy,velino);
    }
    else
      tmp.MultiplyTN(derxy_,velino);

    const double val = tmp.Norm1();
    h = 2.0/val; // h=streamlength

    elevec1[0] = h/vel_norm;
  }
  else
    elevec1[0] = 1.0e12;

  return 0;
}

/*-----------------------------------------------------------------------------*
 | Calculate properties for adaptive forcing of periodic hill       bk 12/2014 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcMassFlowPeriodicHill(
    DRT::ELEMENTS::Fluid*                ele,
    Teuchos::ParameterList&   params,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm,
    Epetra_SerialDenseVector&       elevec1,
    Teuchos::RCP<MAT::Material> &        mat
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

  GetMaterialParams(mat, mat1, mat2, mat2, mat2, mat2, mat2, 0.0, 0.0, 0.0, 0.0, 0.0);

  // ---------------------------------------------------------------------------
  // Geometry
  // ---------------------------------------------------------------------------
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);

  // Do ALE specific updates if necessary
  if (ele->IsAle())
    dserror("no ale for periodic hill");

  LINALG::Matrix<nsd_,nen_> evelnp(true);
  ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evelnp, NULL,"velnp");

  // definition of matrices
  LINALG::Matrix<nen_*nsd_,nen_*nsd_> estif_u(true);

  //length of whole domain
  double length = params.get<double>("length");

  // ---------------------------------------------------------------------------
  // Integration loop
  // ---------------------------------------------------------------------------
  double massf=0.0;
  for ( DRT::UTILS::GaussIntegration::iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
  {
    // evaluate shape functions and derivatives at integration point
    EvalShapeFuncAndDerivsAtIntPoint(iquad.Point(),iquad.Weight());

    // create dummy matrices
    LINALG::Matrix<nsd_, nen_> mat1(true);
    LINALG::Matrix<nen_, 1> mat2(true);

    GetMaterialParams(mat, mat1, mat2, mat2, mat2, mat2, mat2, 0.0, 0.0, 0.0, 0.0, 0.0);

    velint_.Multiply(evelnp,funct_);

    massf += velint_(0)*densaf_*fac_;
  }
  massf /= length;
  elevec1[0] = massf;

  return 0;
}

  /*-----------------------------------------------------------------------------*
   | Reset Immersed Ele                                            rauch 05/2014 |
   *-----------------------------------------------------------------------------*/
  template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
  int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::ResetImmersedEle(
      DRT::ELEMENTS::Fluid*           ele,
      Teuchos::ParameterList&         params
     )
  {
  DRT::ELEMENTS::FluidImmersedBase* immersedele = dynamic_cast<DRT::ELEMENTS::FluidImmersedBase*>(ele);

  // reset element information
  immersedele->SetIsImmersed(0);
  immersedele->SetBoundaryIsImmersed(0);
  immersedele->SetHasProjectedDirichlet(0);

  // reset node information
  DRT::Node** nodes = immersedele->Nodes();
  for(int i=0;i<immersedele->NumNode();++i)
  {
    static_cast<IMMERSED::ImmersedNode* >(nodes[i])->SetIsMatched(0);
    static_cast<IMMERSED::ImmersedNode* >(nodes[i])->SetIsBoundaryImmersed(0);
  }

  // reset element int point information
  if(immersedele->GetRCPProjectedIntPointDivergence() != Teuchos::null)
    immersedele->DestroyElementRCP();

  return 0;
}

/*-----------------------------------------------------------------------------*
 | Calculate coordinates and velocities and element center          bk 01/2015 |
 *-----------------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, DRT::ELEMENTS::Fluid::EnrichmentType enrtype>
int DRT::ELEMENTS::FluidEleCalc<distype,enrtype>::CalcVelGradientEleCenter(
    DRT::ELEMENTS::Fluid*                ele,
    DRT::Discretization &                discretization,
    const std::vector<int> &             lm,
    Epetra_SerialDenseVector&            elevec1,
    Epetra_SerialDenseVector&            elevec2
    )
{
  if(distype!=DRT::Element::hex8 && distype!=DRT::Element::tet4 && distype!=DRT::Element::quad4 && distype!=DRT::Element::tri3)
    dserror("this is currently only implemented for linear elements");
  // get node coordinates
  GEO::fillInitialPositionArray<distype,nsd_, LINALG::Matrix<nsd_,nen_> >(ele,xyze_);
  // Do ALE specific updates if necessary
  if (ele->IsAle())
  {
    LINALG::Matrix<nsd_,nen_> edisp(true);
    ExtractValuesFromGlobalVector(discretization, lm, *rotsymmpbc_, &edisp, NULL, "disp");

    // get new node positions of ALE mesh
     xyze_ += edisp;
  }

  // evaluate shape functions and derivatives element center
  EvalShapeFuncAndDerivsAtEleCenter();

  if(discretization.HasState("vel"))
  {
    // extract element velocities
    LINALG::Matrix<nsd_,nen_> evel(true);
    ExtractValuesFromGlobalVector(discretization,lm, *rotsymmpbc_, &evel, NULL,"vel");

    // get gradient of velocity at element center
    vderxy_.MultiplyNT(evel,derxy_);

    // write values into element vector (same order as in VelGradientProjection())
    for (int i=0; i<nsd_; ++i) // loop rows of vderxy
    {
      for (int j=0; j<nsd_; ++j) // loop columns of vderxy
      {
        elevec1[i*nsd_+j] = vderxy_(i,j);
      }
    }
  }

  // get position of element centroid
  LINALG::Matrix<nsd_,1> x_centroid(true);
  x_centroid.Multiply(xyze_,funct_);
  for (int i=0; i<nsd_; ++i)
  {
    elevec2[i] = x_centroid(i);
  }

  return 0;
}



// template classes
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex8,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex8,DRT::ELEMENTS::Fluid::xwall>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet4,DRT::ELEMENTS::Fluid::xwall>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex20,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::hex27,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet4,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tet10,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::wedge6,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::wedge15,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::pyramid5,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad4,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad8,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::quad9,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tri3,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::tri6,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::nurbs9,DRT::ELEMENTS::Fluid::none>;
template class DRT::ELEMENTS::FluidEleCalc<DRT::Element::nurbs27,DRT::ELEMENTS::Fluid::none>;
