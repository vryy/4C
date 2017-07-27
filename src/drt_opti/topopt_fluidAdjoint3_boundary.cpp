/*!------------------------------------------------------------------------------------------------*
\file topopt_fluidAdjoint3_boundary.cpp

\brief boundary element implementation of fluid adjoint equations for topology optimization

\level 2

<pre>
\maintainer Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_fluidAdjoint3_boundary.H"
#include "topopt_fluidAdjoint3_impl_parameter.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_lib/drt_element_integration_select.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_geometry/position_array.H"
#include "../drt_inpar/inpar_topopt.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/standardtypes_cpp.H"
#include "../drt_mat/newtonianfluid.H"

#include "../drt_cut/cut_position.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidAdjoint3BoundaryImplInterface* DRT::ELEMENTS::FluidAdjoint3BoundaryImplInterface::Impl(const DRT::Element* ele)
{
  switch (ele->Shape())
  {
  case DRT::Element::quad4:
  {
    return FluidAdjoint3BoundaryImpl<DRT::Element::quad4>::Instance();
  }
  case DRT::Element::quad8:
  {
    return FluidAdjoint3BoundaryImpl<DRT::Element::quad8>::Instance();
  }
  case DRT::Element::quad9:
  {
    return FluidAdjoint3BoundaryImpl<DRT::Element::quad9>::Instance();
  }
//  case DRT::Element::tri3:
//  {
//    return FluidAdjoint3BoundaryImpl<DRT::Element::tri3>::Instance();
//  }
//  case DRT::Element::tri6:
//  {
//    return FluidAdjoint3BoundaryImpl<DRT::Element::tri6>::Instance();
//  }
  case DRT::Element::line2:
  {
    return FluidAdjoint3BoundaryImpl<DRT::Element::line2>::Instance();
  }
  case DRT::Element::line3:
  {
    return FluidAdjoint3BoundaryImpl<DRT::Element::line3>::Instance();
  }
  default:
  {
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
    break;
  }
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidAdjoint3BoundaryImpl<distype> * DRT::ELEMENTS::FluidAdjoint3BoundaryImpl<distype>::Instance(
    bool create)
{
  static FluidAdjoint3BoundaryImpl<distype> * instance;
  if (create)
  {
    if (instance==NULL)
      instance = new FluidAdjoint3BoundaryImpl<distype>();
  }
  else
  {
    if (instance!=NULL)
      delete instance;
    instance = NULL;
  }
  return instance;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3BoundaryImpl<distype>::Done()
{
  // delete this pointer! Afterwards we have to go! But since this is a
  // cleanup call, we can do it this way.
  Instance( false );
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidAdjoint3BoundaryImpl<distype>::FluidAdjoint3BoundaryImpl()
  : intpoints_(distype),
    xyze_(true),
    funct_(true),
    deriv_(true),
    unitnormal_(true),
    velint_(true),
    velint_old_(true),
    fluidvelint_(true),
    fluidvelint_old_(true),
    pressint_(0.0),
    pressint_old_(0.0),
    drs_(0.0),
    fac_(0.0),
    visc_(0.0),
    dens_(1.0)
{
  // pointer to class FluidImplParameter (access to the general parameter)
  fldAdPara_ = DRT::ELEMENTS::FluidAdjoint3ImplParameter::Instance();

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition     winklmaier 03/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidAdjoint3BoundaryImpl<distype>::EvaluateNeumann(
                              DRT::ELEMENTS::FluidBoundary*  ele,
                              Teuchos::ParameterList&        params,
                              DRT::Discretization&           discretization,
                              const std::vector<int>&        lm,
                              Epetra_SerialDenseVector&      elevec)
{
  // reshape element vector
  elevec.Shape(numdofpernode_*nen_,1);
  // initialize to zero
  elevec.Scale(0.0);

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

  Teuchos::RCP<DRT::Condition> condition = params.get<Teuchos::RCP<DRT::Condition> >("condition",Teuchos::null);

  // get values, switches and spatial functions from the condition
  // (assumed to be constant on element boundary)
  const std::vector<int>*    onoff = condition->Get<std::vector<int> >   ("onoff");

  // get time factor for Neumann term
  const double timefac = fldAdPara_->Timefac();
  const double timefacrhs = fldAdPara_->TimefacRhs();

  // get Gaussrule
  const DRT::UTILS::IntPointsAndWeights<bdrynsd_> intpoints(DRT::ELEMENTS::DisTypeToOptGaussRule<distype>::rule);

  // get local node coordinates
  // (we have a nsd_ dimensional domain, since nsd_ determines the dimension of the parent element!)
  GEO::fillInitialPositionArray<distype,nsd_,LINALG::Matrix<nsd_,nen_> >(ele,xyze_);


  /*----------------------------------------------------------------------*
  |               start loop over integration points                     |
  *----------------------------------------------------------------------*/
  for ( DRT::UTILS::GaussIntegration::const_iterator iquad=intpoints_.begin(); iquad!=intpoints_.end(); ++iquad )
  {
    if (params.get<INPAR::TOPOPT::AdjointCase>("special test case") == INPAR::TOPOPT::adjointtest_no)
    {
      ; // boundary terms are currently independent of velocity -> no entry here
    }
    else // special cases
    {
      // evaluate shape functions and their derivatives,
      // compute unit normal vector and infinitesimal area element drs
      EvalShapeFuncAtBouIntPoint(iquad,ele->Id());

      // time factors
      const double timefacfac = fac_*timefac;
      const double timefacfacrhs = fac_*timefacrhs;

      INPAR::TOPOPT::AdjointCase testcase = params.get<INPAR::TOPOPT::AdjointCase>("special test case");

      // get global coordinates of gauss point
      double x = 0.0;
      double y = 0.0;
      LINALG::Matrix<nsd_,1> coords(true);
      coords.Multiply(xyze_,funct_);
      x = coords(0);
      y = coords(1); // z-component currently not required in tests

      LINALG::Matrix<nsd_,1> values(true);
      LINALG::Matrix<nsd_,1> values_old(true);

      switch (testcase)
      {
      case INPAR::TOPOPT::adjointtest_stat_const_vel_lin_pres:
      {
        values(0) = values_old(0) = -9;
        values(1) = values_old(1) = 0;
        break;
      }
      case INPAR::TOPOPT::adjointtest_stat_lin_vel_quad_pres:
      {
        values(0) = values_old(0) = 10 + 0.5*(125010*x*x + 175000*y*y + 100004*x*y);
        values(1) = values_old(1) = 3*x*x + 7*x*y + 5;
        break;
      }
      case INPAR::TOPOPT::adjointtest_stat_quad_vel_lin_pres:
      {
        values(0) = values_old(0) = x*x*x + 2*x*x*y - 2*y;
        values(1) = values_old(1) = - 3*x*x*x + 4*y*y*y - 6*x*x*y + 2*x*y*y - 6*x;
        break;
      }
      case INPAR::TOPOPT::adjointtest_stat_all_terms_all_constants:
      {
        values(0) = values_old(0) = 2*x*x*x + 6*x*x*y + 4*x*y*y - 3*y;
        values(1) = values_old(1) = - 6*x*x*x - 12*x*x*y + 4*x*y*y + 8*y*y*y - 15*x;
        break;
      }
      case INPAR::TOPOPT::adjointtest_instat_varying_theta:
      {
        double t = fldAdPara_->Time();
        values(0) = 3*x - 3*t + 4 - 5*x*t + 10*y*t;
        values(1) = -3*y + 6*t;

        t += fldAdPara_->Dt(); // old time = t + dt
        values_old(0) = 3*x - 3*t + 4 - 5*x*t + 10*y*t;
        values_old(1) = -3*y + 6*t;
        break;
      }
      case INPAR::TOPOPT::adjointtest_instat_all_terms_all_constants:
      {
        double t = fldAdPara_->Time();
        values(0) = 2*x*x*x + 4*x*x*y + 2*x*x*y*t + 4*x*y*y*t + 18*y + 3*y*t - 24*y*t*t;
        values(1) = - 6*x*x*x - 12*x*x*y + 4*x*y*y*t*t + 8*y*y*y*t*t - 18*x + 3*x*t;

        t += fldAdPara_->Dt(); // old time = t + dt
        values_old(0) = 2*x*x*x + 4*x*x*y + 2*x*x*y*t + 4*x*y*y*t + 18*y + 3*y*t - 24*y*t*t;
        values_old(1) = - 6*x*x*x - 12*x*x*y + 4*x*y*y*t*t + 8*y*y*y*t*t - 18*x + 3*x*t;
        break;
      }
      case INPAR::TOPOPT::adjointtest_instat_primal_and_dual:
      {
        double t = fldAdPara_->Time();
        values(0) = 3*x*y*t*t + 3*x*x*y*t + 5.5*x*x + 6*x*y*y*t + 6*x*y + 2*y*t + 4 - 1.5*y*y + 3*y*t*t + 5.5*y*y*t;
        values(1) = 3*y*t*t - 3*x*t*t + 6*y*y*t - 3*x*y*t - 2*t + 2*x*t - 3*x*x*t;

        t += fldAdPara_->Dt(); // old time = t + dt
        values_old(0) = 3*x*y*t*t + 3*x*x*y*t + 5.5*x*x + 6*x*y*y*t + 6*x*y + 2*y*t + 4 - 1.5*y*y + 3*y*t*t + 5.5*y*y*t;
        values_old(1) = 3*y*t*t - 3*x*t*t + 6*y*y*t - 3*x*y*t - 2*t + 2*x*t - 3*x*x*t;
        break;
      }
      case INPAR::TOPOPT::adjointtest_primal:
        break;
      default:
      {
        dserror("no dirichlet condition implemented for special test case");
        break;
      }
      }


      for (int jdim=0;jdim<nsd_;++jdim)
      {
        if((*onoff)[jdim])  // Is this dof activated
        {
          const double functval = timefacfac*values(jdim);

          for (int vi=0;vi<nen_;++vi)
          {
            elevec(vi*numdofpernode_+jdim) += funct_(vi)*functval;
          }

          if (not fldAdPara_->IsStationary())
          {
            const double functval_old = timefacfacrhs*values_old(jdim);

            for(int vi=0; vi < nen_; ++vi )
            {
              elevec(numdofpernode_*vi+jdim) += funct_(vi)*functval_old;
            }
          }
        }
      }
    }
  }

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3BoundaryImpl<distype>::EvalShapeFuncAtBouIntPoint(
    DRT::UTILS::GaussIntegration::iterator & iquad,       // actual integration point
    const int                              eleid        // element ID
)
{
  // local gauss point coordinates
  xsi_ = LINALG::Matrix<bdrynsd_,1>(iquad.Point());

  // shape functions and their first derivatives of boundary element
  DRT::UTILS::shape_function<distype>(xsi_,funct_);
  DRT::UTILS::shape_function_deriv1<distype>(xsi_,deriv_);

  // compute measure tensor for surface element, infinitesimal area element drs
  // and (outward-pointing) unit normal vector
  LINALG::Matrix<bdrynsd_,bdrynsd_> metrictensor(true);
  DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>(xyze_,deriv_,metrictensor,drs_,&unitnormal_);

  // compute integration factor
  fac_ = iquad.Weight()*drs_;

  return;
}


/*----------------------------------------------------------------------*
 |  call material parameters                           winklmaier 03/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3BoundaryImpl<distype>::GetMaterialParams(
  Teuchos::RCP<const MAT::Material>    material
)
{
// initially set density and density factor for Neumann boundary conditions to 1.0
// (the latter only changed for low-Mach-number flow problems)
if (material->MaterialType() == INPAR::MAT::m_fluid)
{
  const MAT::NewtonianFluid* actmat = static_cast<const MAT::NewtonianFluid*>(material.get());

  // get constant viscosity
  visc_ = actmat->Viscosity();

  dens_ = actmat->Density();
}
else
  dserror("Material type is not supported for boundary element!");

// check whether there is zero or negative (physical) viscosity
if (visc_ < EPS15) dserror("zero or negative (physical) diffusivity");

return;
} // FluidBoundaryImpl::GetMaterialParams



/*!
 * \brief fill elment matrix and vectors with the global values
 */
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3BoundaryImpl<distype>::ExtractValuesFromGlobalVector(
    const DRT::Discretization&      discretization, ///< discretization
    const std::vector<int>&              lm,             ///<
    LINALG::Matrix<nsd_,nen_>*      matrixtofill,   ///< vector field
    LINALG::Matrix<nen_,1>*         vectortofill,   ///< scalar field
    const std::string               state          ///< state of the global vector
) const
{
  // get state of the global vector
  Teuchos::RCP<const Epetra_Vector> matrix_state = discretization.GetState(state);
  if(matrix_state == Teuchos::null)
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
