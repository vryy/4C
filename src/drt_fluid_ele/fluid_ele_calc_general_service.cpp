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
