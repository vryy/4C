/*!------------------------------------------------------------------------------------------------*
\file topopt_fluidAdjoint3_boundary.cpp

\brief boundary element implementation of fluid adjoint equations for topology optimization

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include "topopt_fluidAdjoint3_boundary.H"
#include "topopt_fluidAdjoint3_impl_parameter.H"
#include "../drt_fluid_ele/fluid_ele.H"
#include "../drt_fluid_ele/fluid_ele_utils.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_geometry/position_array.H"
#include "../drt_lib/drt_utils.H"
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
    dserror("Element shape %d (%d nodes) not activated. Just do it.", ele->Shape(), ele->NumNode());
  }
  return NULL;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
template<DRT::Element::DiscretizationType distype>
DRT::ELEMENTS::FluidAdjoint3BoundaryImpl<distype> * DRT::ELEMENTS::FluidAdjoint3BoundaryImpl<distype>::Instance()
{
  static FluidAdjoint3BoundaryImpl<distype> * instance;
  if ( instance==NULL )
    instance = new FluidAdjoint3BoundaryImpl<distype>();
  return instance;
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
  // pointer to class Fluid3ImplParameter (access to the general parameter)
  fluidAdjoint3Parameter_ = DRT::ELEMENTS::FluidAdjoint3ImplParameter::Instance();

  return;
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface Neumann boundary condition     winklmaier 03/12 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype>
int DRT::ELEMENTS::FluidAdjoint3BoundaryImpl<distype>::EvaluateNeumann(
                              DRT::ELEMENTS::Fluid3Boundary* ele,
                              ParameterList&                 params,
                              DRT::Discretization&           discretization,
                              vector<int>&                   lm,
                              Epetra_SerialDenseVector&      elevec)
{
  // the vectors have been allocated outside in EvaluateConditionUsingParentData()
  RCP<vector<int> > plm = params.get<RCP<vector<int> > >("plm");
  RCP<vector<int> > plmowner = params.get<RCP<vector<int> > >("plmowner");
  RCP<vector<int> > plmstride = params.get<RCP<vector<int> > >("plmstride");
  ele->LocationVector(discretization,*plm,*plmowner,*plmstride);

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
  const vector<int>*    onoff = condition->Get<vector<int> >   ("onoff");

  // get time factor for Neumann term
  const double timefac = fluidAdjoint3Parameter_->timefac_;
  const double timefacrhs = fluidAdjoint3Parameter_->timefacrhs_;

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
    // evaluate shape functions and their derivatives,
    // compute unit normal vector and infinitesimal area element drs
    // (evaluation of nurbs-specific stuff not activated here)
    EvalShapeFuncAtBouIntPoint(iquad,ele->Id());

    // get the required material information
    Teuchos::RCP<MAT::Material> material = ele->ParentElement()->Material();

    // get material parameters
    // (evaluation always at integration point, in contrast to parent element)
    GetMaterialParams(material);

    // time factors
    const double timefacfac = fac_*timefac;
    const double timefacfacrhs = fac_*timefacrhs;

    // evaluate state vectors on gauss point
    velint_.Multiply(evelnp,funct_);
    velint_old_.Multiply(eveln,funct_);

    fluidvelint_.Multiply(efluidvelnp,funct_);
    fluidvelint_old_.Multiply(efluidveln,funct_);

    pressint_ = epren.Dot(funct_);
    pressint_old_ = eprenp.Dot(funct_);

    // manual setting for example
    LINALG::Matrix<nsd_,1> coords(true);
    coords.Multiply(xyze_,funct_);

    LINALG::Matrix<nsd_,1> values(true);
    // dens 2, visc 3
    values(0) = 2*coords(0)*coords(0)*coords(0) + 6*coords(0)*coords(0)*coords(1)
                + 4*coords(0)*coords(1)*coords(1) - 3*coords(1);
    values(1) = - 6*coords(0)*coords(0)*coords(0) - 12*coords(0)*coords(0)*coords(1)
                + 4*coords(0)*coords(1)*coords(1) + 8*coords(1)*coords(1)*coords(1)
                - 15*coords(0);
    // quadratic velocity, linear pressure, exact example
//    values(0) = coords(0)*coords(0)*coords(0) + 2*coords(0)*coords(0)*coords(1) - 2*coords(1);
//    values(1) = - 3*coords(0)*coords(0)*coords(0) + 4*coords(1)*coords(1)*coords(1)
//                - 6*coords(0)*coords(0)*coords(1) + 2*coords(0)*coords(1)*coords(1)
//                - 6*coords(0);
    // quadratic exact example
//    values(0) = 10.0 + 0.5*(125010.0*coords(0)*coords(0)
//                          +175000.0*coords(1)*coords(1)
//                          +100004.0*coords(0)*coords(1));
//    values(1) = 3.0*coords(0)*coords(0) + 7.0*coords(0)*coords(01) + 5.0;
    // basic example : alpha = 0, no convection, no diffusion
//    values(0) = -9.0;
//    values(1) = 0.0;
    // most fully example
//    values(0) = -12500*(coords(1)*coords(1)-coords(0)*coords(0)) + 2.0;
//    values(1) = -coords(1);

    LINALG::Matrix<nsd_,1> values_old(true); // TODO fill for testing

    for (int jdim=0;jdim<nsd_;++jdim)
    {
      if((*onoff)[jdim])  // Is this dof activated
      {
        const double functval = timefacfac*values(jdim);

        for (int vi=0;vi<nen_;++vi)
        {
          elevec(vi*numdofpernode_+jdim) += funct_(vi)*functval;
        }

        if (not fluidAdjoint3Parameter_->is_stationary_)
        {
          const double functval_old = timefacfacrhs*values_old(jdim);

          for(int vi=0; vi < nen_; ++vi )
          {
            elevec(numdofpernode_*vi+jdim) += funct_(vi)*functval_old;
          }
        }
      }
    }
    /* convective boundary term */
    /*
                     /                       \
                    |   /  n     \            |
            - rho * |  |  u  o n  | Dv ,   w  |
                    |   \        /    (i)     |
                     \                       /
    */
//    for (int idim=0;idim<nsd_;++idim)
//    {
//      if((*onoff)[idim])  // Is this dof activated
//      {
//        // system matrix entry
//        double value = 0.0;
//
//        for (int dim=0;dim<nsd_;++dim)
//          value += fluidvelint_(dim)*unitnormal_(dim);
//
//        value*=timefacfac*dens_;
//
//        for(int ui=0; ui < nen_; ++ui )
//        {
//          double functval = value*funct_(ui);
//
//          for (int vi=0;vi<nen_;vi++)
//          {
//            elemat(vi*numdofpernode_+idim,ui*numdofpernode_+idim) -= funct_(vi)*functval;
//
//          }
//        }
//
//        // right hand side at new time step
//        for(int vi=0; vi < nen_; ++vi )
//        {
//          elevec(numdofpernode_*vi+idim,0) += funct_(vi)*value*velint_(idim);
//        }
//
//        // right hand side at old time step
//        if (not fluidAdjoint3Parameter_->is_stationary_)
//        {
//          value = 0.0;
//
//          for (int dim=0;dim<nsd_;++dim)
//            value += fluidvelint_old_(dim)*unitnormal_(dim);
//
//          value*=timefacfacrhs;
//
//          for(int vi=0; vi < nen_; ++vi )
//          {
//            elevec(numdofpernode_*vi+idim,0) += funct_(vi)*value*velint_old_(idim);
//          }
//        }
//      }  // if (*onoff)
//    }
//
//
    /* pressure boundary term */
    /*
                     /          \
                    |            |
                  + |   Dq ,  w  |
                    |            |
                     \          /
    */
//    for(int vi=0; vi < nen_; ++vi )
//    {
//      double value = timefacfac*funct_(vi);
//
//      for (int idim = 0; idim<nsd_;idim++)
//      {
//        if((*onoff)[idim])  // Is this dof activated
//        {
//          double functval = value*unitnormal_(idim);
//
//          for (int ui=0;ui<nen_;ui++)
//          {
//            elemat(vi*numdofpernode_+idim,ui*numdofpernode_+idim) += functval*funct_(ui);
//          }
//
//          elevec(numdofpernode_*vi+idim,0) -= functval*pressint_;
//
//          if (not fluidAdjoint3Parameter_->is_stationary_)
//          {
//            elevec(numdofpernode_*vi+idim,0) -= timefacfacrhs*funct_(vi)*unitnormal_(idim)*pressint_old_;
//          }
//        }
//      }
//    }
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
// (the latter only changed for low-Mach-number flow/combustion problems)
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
} // Fluid3BoundaryImpl::GetMaterialParams



/*!
 * \brief fill elment matrix and vectors with the global values
 */
template <DRT::Element::DiscretizationType distype>
void DRT::ELEMENTS::FluidAdjoint3BoundaryImpl<distype>::ExtractValuesFromGlobalVector(
    const DRT::Discretization&      discretization, ///< discretization
    const vector<int>&              lm,             ///<
    LINALG::Matrix<nsd_,nen_>*  matrixtofill,   ///< vector field
    LINALG::Matrix<nen_,1>*     vectortofill,   ///< scalar field
    const std::string               state          ///< state of the global vector
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
