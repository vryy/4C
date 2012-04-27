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
#include "fluid_ele_boundary_calc.H"
#include "fluid_ele_calc_weak_dbc.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_opti/topopt_fluidAdjoint3_boundary.H"


/*---------------------------------------------------------------------*
|  converts a string into an action for this element                   |
 *----------------------------------------------------------------------*/
DRT::ELEMENTS::FluidBoundary::ActionType DRT::ELEMENTS::FluidBoundary::convertStringToActionType(
    const string& action) const
{
  dsassert(action != "none", "No action supplied");

  DRT::ELEMENTS::FluidBoundary::ActionType act = FluidBoundary::none;
  if (action == "none") dserror("No action supplied");
  else if (action == "integrate_Shapefunction")
    act = FluidBoundary::integrate_Shapefunction;
  else if (action == "area calculation")
    act = FluidBoundary::areacalc;
  else if (action == "calc_flowrate")
    act = FluidBoundary::calc_flowrate;
  else if (action == "flowrate_deriv")
    act = FluidBoundary::flowratederiv;
  else if (action == "Outlet impedance")
    act = FluidBoundary::Outletimpedance;
  else if (action == "calc_node_normal")
    act = FluidBoundary::calc_node_normal;
  else if (action == "calc_node_curvature")
    act = FluidBoundary::calc_node_curvature;
  else if (action == "calc_surface_tension")
    act = FluidBoundary::calc_surface_tension;
  else if (action == "enforce_weak_dbc")
    act = FluidBoundary::enforce_weak_dbc;
  else if (action == "MixedHybridDirichlet")
    act = FluidBoundary::mixed_hybrid_dbc;
  else if (action == "conservative_outflow_bc")
    act = FluidBoundary::conservative_outflow_bc;
  else if (action == "calc_Neumann_inflow")
    act = FluidBoundary::calc_Neumann_inflow;
  else if (action == "calculate integrated pressure")
    act = FluidBoundary::integ_pressure_calc;
  else if (action == "center of mass calculation")
    act = FluidBoundary::center_of_mass_calc;
  else if (action == "calculate traction velocity component")
    act = FluidBoundary::traction_velocity_component;
  else if (action == "calculate Uv integral component")
    act = FluidBoundary::traction_Uv_integral_component;
  else if (action == "no penetration")
    act = FluidBoundary::no_penetration;
  else if (action == "AdjointNeumannBoundaryCondition")
    act = FluidBoundary::adjoint_neumann;
  else
    dserror("Unknown type of action for Fluid_Boundary: %s",action.c_str());

  return act;
}


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
  const string action = params.get<string>("action","none");
  const DRT::ELEMENTS::FluidBoundary::ActionType act = convertStringToActionType(action);



  // get status of Ale
  const bool isale = this->ParentElement()->IsAle();

  switch(act)
  {
  case FluidBoundary::integrate_Shapefunction:
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
  case FluidBoundary::areacalc:
  {
    if (this->Owner() == discretization.Comm().MyPID())
      DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->AreaCaculation(
          this,
          params,
          discretization,
          lm);
    break;
  }
  case FluidBoundary::integ_pressure_calc:
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
  case FluidBoundary::calc_flowrate:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->ComputeFlowRate(
        this,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case FluidBoundary::flowratederiv:
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
  case FluidBoundary::Outletimpedance:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->ImpedanceIntegration(
        this,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case FluidBoundary::calc_node_normal:
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
  case FluidBoundary::calc_node_curvature:
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
  case FluidBoundary::enforce_weak_dbc:
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
  case FluidBoundary::mixed_hybrid_dbc:
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
  case FluidBoundary::conservative_outflow_bc:
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
  case FluidBoundary::calc_Neumann_inflow:
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
  case FluidBoundary::calc_surface_tension:
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
  case FluidBoundary::center_of_mass_calc:
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
  case FluidBoundary::traction_velocity_component:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->CalcTractionVelocityComponent(
        this,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case FluidBoundary::traction_Uv_integral_component:
  {
    DRT::ELEMENTS::FluidBoundaryImplInterface::Impl(this)->ComputeNeumannUvIntegral(
        this,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case FluidBoundary::no_penetration:
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
  case FluidBoundary::adjoint_neumann:
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
    dserror("Unknown type of action for Fluid_Boundary: %s",action.c_str());
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


