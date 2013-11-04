/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_boundary_evaluate.cpp

\brief Evaluate boundary conditions for fluid problems

<pre>
Maintainer: Ursula Rasthofer & Volker Gravemeier
            {rasthofer,vgravem}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236/-245
</pre>
 */
/*----------------------------------------------------------------------*/


#include "fluid_ele.H"
#include "fluid_ele_action.H"
#include "fluid_ele_boundary_calc.H"
#include "fluid_ele_boundary_parent_calc.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_opti/topopt_fluidAdjoint3_boundary.H"


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidBoundary::Evaluate(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::BoundaryAction act = DRT::INPUT::get<FLD::BoundaryAction>(params,"action");

  // get status of Ale
  const bool isale = this->ParentElement()->IsAle();

  switch(act)
  {
  case FLD::integrate_Shapefunction:
  {
    RCP<const Epetra_Vector> dispnp;
    std::vector<double> mydispnp;

    if (isale)
    {
      dispnp = discretization.GetState("dispnp");
      if (dispnp!=Teuchos::null)
      {
        mydispnp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
      }
    }

    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->IntegrateShapeFunction(
        this,
        params,
        discretization,
        lm,
        elevec1,
        mydispnp);
    break;
  }
  case FLD::calc_area:
  {
    if (this->Owner() == discretization.Comm().MyPID())
      DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->AreaCalculation(
          this,
          params,
          discretization,
          lm);
    break;
  }
  case FLD::calc_pressure_bou_int:
  {
    if(this->Owner() == discretization.Comm().MyPID())
      DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->PressureBoundaryIntegral(
          this,
          params,
          discretization,
          lm);
    break;
  }
  // general action to calculate the flow rate
  case FLD::calc_flowrate:
  {
    // pointer to class FluidEleParameter (access to the general parameter)
    Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance();
    switch(fldpara->PhysicalType())
    {
    case INPAR::FLUID::incompressible:
    {
      DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->ComputeFlowRate(
          this,
          params,
          discretization,
          lm,
          elevec1);
      break;
    }
    case INPAR::FLUID::poro:
    {
      DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->PoroFlowRate(
          this,
          params,
          discretization,
          lm,
          elevec1);
      break;
    }
    default:
      dserror("action 'calc_flowrate' not implemented for this physical type");
    break;
    }
    break;
  }
  case FLD::flowratederiv:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->FlowRateDeriv(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elemat2,
        elevec1,
        elevec2,
        elevec3);
    break;
  }
  case FLD::Outletimpedance:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->ImpedanceIntegration(
        this,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case FLD::ba_calc_node_normal:
  {
    RCP<const Epetra_Vector> dispnp;
    std::vector<double> mydispnp;

    if (isale)
    {
      dispnp = discretization.GetState("dispnp");
      if (dispnp!=Teuchos::null)
      {
        mydispnp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
      }
    }
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->ElementNodeNormal(
        this,
        params,
        discretization,
        lm,
        elevec1,
        mydispnp);
    break;
  }
  case FLD::calc_node_curvature:
  {
    RCP<const Epetra_Vector> dispnp;
    std::vector<double> mydispnp;

    if (isale)
    {
      dispnp = discretization.GetState("dispnp");
      if (dispnp!=Teuchos::null)
      {
        mydispnp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
      }
    }

    RCP<const Epetra_Vector> normals;
    std::vector<double> mynormals;

    normals = discretization.GetState("normals");
    if (normals!=Teuchos::null)
    {
      mynormals.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*normals,mynormals,lm);
    }

    // what happens, if the mynormals vector is empty? (ehrl)
    dserror("the action calc_node_curvature has not been called by now. What happens, if the mynormal vector is empty");

    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->ElementMeanCurvature(
        this,
        params,
        discretization,
        lm,
        elevec1,
        mydispnp,
        mynormals);
    break;
  }
  case FLD::flow_dep_pressure_bc:
  {
    DRT::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->FlowDepPressureBC(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case FLD::enforce_weak_dbc:
  {
    DRT::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->EvaluateWeakDBC(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case FLD::evaluate_nitsche_par:
  {
    DRT::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->EvaluateNitschePar(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elemat2);
    break;
  }
  case FLD::mixed_hybrid_dbc:
  {
    DRT::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->MixHybDirichlet(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case FLD::conservative_outflow_bc:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->ConservativeOutflowConsistency(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case FLD::calc_Neumann_inflow:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->NeumannInflow(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case FLD::calc_surface_tension:
  {
    // employs the divergence theorem acc. to Saksono eq. (24) and does not
    // require second derivatives.

    RCP<const Epetra_Vector> dispnp;
    std::vector<double> mydispnp;

    dispnp = discretization.GetState("dispnp");
    if (dispnp!=Teuchos::null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }

    // mynormals and mycurvature are not used in the function
    std::vector<double> mynormals;
    std::vector<double> mycurvature;

    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->ElementSurfaceTension(
        this,
        params,
        discretization,
        lm,
        elevec1,
        mydispnp,
        mynormals,
        mycurvature);
    break;
  }
  case FLD::center_of_mass_calc:
  {
    // evaluate center of mass
    if(this->Owner() == discretization.Comm().MyPID())
      DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->CenterOfMassCalculation(
          this,
          params,
          discretization,
          lm);
    break;
  }
  case FLD::traction_velocity_component:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->CalcTractionVelocityComponent(
        this,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case FLD::traction_Uv_integral_component:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->ComputeNeumannUvIntegral(
        this,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case FLD::no_penetration:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->NoPenetration(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elemat2,
        elevec1,
        elevec2);
    break;
  }
  case FLD::no_penetrationIDs:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->NoPenetrationIDs(
        this,
        params,
        discretization,
        elevec1,
        lm);
    break;
  }
  case FLD::poro_boundary:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->PoroBoundary(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case FLD::poro_prescoupl:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->PressureCoupling(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case FLD::ba_calc_adjoint_neumann:
  {
    DRT::ELEMENTS::FluidAdjoint3BoundaryImplInterface::Impl(this)->EvaluateNeumann(
        this,
        params,
        discretization,
        lm,
        elevec1
    );
    break;
  }
  case FLD::fpsi_coupling:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->FPSICoupling(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  default:
    dserror("Unknown type of action for Fluid_Boundary!");
    break;
  } // end of switch(act)

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition       gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidBoundary::EvaluateNeumann(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  return DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->EvaluateNeumann(
      this,
      params,
      discretization,
      condition,
      lm,
      elevec1,
      elemat1);
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidBoundary::LocationVector(
    const Discretization&    dis,
    LocationArray&           la,
    bool                     doDirichlet,
    const std::string&       condstring,
    Teuchos::ParameterList&  params
    ) const
{
  // get the action required
  const FLD::BoundaryAction act = DRT::INPUT::get<FLD::BoundaryAction>(params,"action");

  switch(act)
  {
  case FLD::enforce_weak_dbc:
  case FLD::poro_boundary:
  case FLD::mixed_hybrid_dbc:
  case FLD::flow_dep_pressure_bc:
  case FLD::fpsi_coupling:
    // special cases: the boundary element assembles also into
    // the inner dofs of its parent element
    // note: using these actions, the element will get the parent location vector
    //       as input in the respective evaluate routines
    parent_->LocationVector(dis,la,doDirichlet);
    break;
  case FLD::calc_flowrate:
  {
    // pointer to class FluidEleParameter (access to the general parameter)
    Teuchos::RCP<DRT::ELEMENTS::FluidEleParameter> fldpara = DRT::ELEMENTS::FluidEleParameter::Instance();
    switch(fldpara->PhysicalType())
    {
    case INPAR::FLUID::incompressible:
    {
      DRT::Element::LocationVector(dis,la,doDirichlet);
      break;
    }
    case INPAR::FLUID::poro:
    {
      // special cases: the boundary element assembles also into
      // the inner dofs of its parent element
      // note: using these actions, the element will get the parent location vector
      //       as input in the respective evaluate routines
      parent_->LocationVector(dis,la,doDirichlet);
      break;
    }
    default:
      dserror("action 'calc_flowrate' not implemented for this physical type");
    break;
    }
    break;
  }
  case FLD::ba_none:
    dserror("No action supplied");
    break;
  default:
    // standard case: element assembles into its own dofs only
    DRT::Element::LocationVector(dis,la,doDirichlet);
    break;
  }
  return;
}

