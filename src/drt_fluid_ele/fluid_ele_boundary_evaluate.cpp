/*----------------------------------------------------------------------*/
/*!
\file fluid_ele_boundary_evaluate.cpp

\brief Evaluate boundary conditions for fluid problems

<pre>
Maintainer: Volker Gravemeier & Andreas Ehrl
            {vgravem,ehrl}@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15245/15252
</pre>
 */
/*----------------------------------------------------------------------*/


#include "fluid_ele.H"
#include "fluid_ele_action.H"
#include "fluid_ele_boundary_calc.H"
#include "fluid_ele_calc_weak_dbc.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_opti/topopt_fluidAdjoint3_boundary.H"


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidBoundary::Evaluate(
    ParameterList&            params,
    DRT::Discretization&      discretization,
    vector<int>&              lm,
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
    RefCountPtr<const Epetra_Vector> dispnp;
    vector<double> mydispnp;

    if (isale)
    {
      dispnp = discretization.GetState("dispnp");
      if (dispnp!=null)
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
  case FLD::areacalc:
  {
    if (this->Owner() == discretization.Comm().MyPID())
      DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->AreaCaculation(
          this,
          params,
          discretization,
          lm);
    break;
  }
  case FLD::integ_pressure_calc:
  {
    if(this->Owner() == discretization.Comm().MyPID())
      DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->IntegratedPressureParameterCalculation(
          this,
          params,
          discretization,
          lm);
    break;
  }
  // general action to calculate the flow rate
  case FLD::calc_flowrate:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->ComputeFlowRate(
        this,
        params,
        discretization,
        lm,
        elevec1);
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
    RefCountPtr<const Epetra_Vector> dispnp;
    vector<double> mydispnp;

    if (isale)
    {
      dispnp = discretization.GetState("dispnp");
      if (dispnp!=null)
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
    RefCountPtr<const Epetra_Vector> dispnp;
    vector<double> mydispnp;

    if (isale)
    {
      dispnp = discretization.GetState("dispnp");
      if (dispnp!=null)
      {
        mydispnp.resize(lm.size());
        DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
      }
    }

    RefCountPtr<const Epetra_Vector> normals;
    vector<double> mynormals;

    normals = discretization.GetState("normals");
    if (normals!=null)
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
  case FLD::enforce_weak_dbc:
  {
    return DRT::ELEMENTS::FluidBoundaryWeakDBCInterface::Impl(this)->EvaluateWeakDBC(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case FLD::mixed_hybrid_dbc:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->MixHybDirichlet(
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
    vector<double> mydispnp;

    dispnp = discretization.GetState("dispnp");
    if (dispnp!=null)
    {
      mydispnp.resize(lm.size());
      DRT::UTILS::ExtractMyValues(*dispnp,mydispnp,lm);
    }

    // mynormals and mycurvature are not used in the function
    vector<double> mynormals;
    vector<double> mycurvature;

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
        elevec1);
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
  default:
    dserror("Unknown type of action for Fluid_Boundary!");
  } // end of switch(act)

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition       gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidBoundary::EvaluateNeumann(
    ParameterList&            params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
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


