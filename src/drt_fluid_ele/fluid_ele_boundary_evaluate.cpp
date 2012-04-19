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
DRT::ELEMENTS::Fluid3Boundary::ActionType DRT::ELEMENTS::Fluid3Boundary::convertStringToActionType(
    const string& action) const
{
  dsassert(action != "none", "No action supplied");

  DRT::ELEMENTS::Fluid3Boundary::ActionType act = Fluid3Boundary::none;
  if (action == "none") dserror("No action supplied");
  else if (action == "integrate_Shapefunction")
    act = Fluid3Boundary::integrate_Shapefunction;
  else if (action == "area calculation")
    act = Fluid3Boundary::areacalc;
  else if (action == "calc_flowrate")
    act = Fluid3Boundary::calc_flowrate;
  else if (action == "flowrate_deriv")
    act = Fluid3Boundary::flowratederiv;
  else if (action == "Outlet impedance")
    act = Fluid3Boundary::Outletimpedance;
  else if (action == "calc_node_normal")
    act = Fluid3Boundary::calc_node_normal;
  else if (action == "calc_node_curvature")
    act = Fluid3Boundary::calc_node_curvature;
  else if (action == "calc_surface_tension")
    act = Fluid3Boundary::calc_surface_tension;
  else if (action == "enforce_weak_dbc")
    act = Fluid3Boundary::enforce_weak_dbc;
  else if (action == "MixedHybridDirichlet")
    act = Fluid3Boundary::mixed_hybrid_dbc;
  else if (action == "conservative_outflow_bc")
    act = Fluid3Boundary::conservative_outflow_bc;
  else if (action == "calc_Neumann_inflow")
    act = Fluid3Boundary::calc_Neumann_inflow;
  else if (action == "calculate integrated pressure")
    act = Fluid3Boundary::integ_pressure_calc;
  else if (action == "center of mass calculation")
    act = Fluid3Boundary::center_of_mass_calc;
  else if (action == "calculate traction velocity component")
    act = Fluid3Boundary::traction_velocity_component;
  else if (action == "calculate Uv integral component")
    act = Fluid3Boundary::traction_Uv_integral_component;
  else if (action == "no penetration")
    act = Fluid3Boundary::no_penetration;
  else if (action == "AdjointNeumannBoundaryCondition")
    act = Fluid3Boundary::adjoint_neumann;
  else
    dserror("Unknown type of action for Fluid3_Boundary: %s",action.c_str());

  return act;
}


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid3Boundary::Evaluate(
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
  const DRT::ELEMENTS::Fluid3Boundary::ActionType act = convertStringToActionType(action);



  // get status of Ale
  const bool isale = this->ParentElement()->IsAle();

  switch(act)
  {
  case Fluid3Boundary::integrate_Shapefunction:
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

    DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->IntegrateShapeFunction(
        this,
        params,
        discretization,
        lm,
        elevec1,
        mydispnp);
    break;
  }
  case Fluid3Boundary::areacalc:
  {
    if (this->Owner() == discretization.Comm().MyPID())
      DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->AreaCaculation(
          this,
          params,
          discretization,
          lm);
    break;
  }
  case Fluid3Boundary::integ_pressure_calc:
  {
    if(this->Owner() == discretization.Comm().MyPID())
      DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->IntegratedPressureParameterCalculation(
          this,
          params,
          discretization,
          lm);
    break;
  }
  // general action to calculate the flow rate
  case Fluid3Boundary::calc_flowrate:
  {
    DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->ComputeFlowRate(
        this,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case Fluid3Boundary::flowratederiv:
  {
    DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->FlowRateDeriv(
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
  case Fluid3Boundary::Outletimpedance:
  {
    DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->ImpedanceIntegration(
        this,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case Fluid3Boundary::calc_node_normal:
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
    DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->ElementNodeNormal(
        this,
        params,
        discretization,
        lm,
        elevec1,
        mydispnp);
    break;
  }
  case Fluid3Boundary::calc_node_curvature:
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

    DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->ElementMeanCurvature(
        this,
        params,
        discretization,
        lm,
        elevec1,
        mydispnp,
        mynormals);
    break;
  }
  case Fluid3Boundary::enforce_weak_dbc:
  {
    return DRT::ELEMENTS::Fluid3BoundaryWeakDBCInterface::Impl(this)->EvaluateWeakDBC(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case Fluid3Boundary::mixed_hybrid_dbc:
  {
    DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->MixHybDirichlet(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case Fluid3Boundary::conservative_outflow_bc:
  {
    DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->ConservativeOutflowConsistency(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case Fluid3Boundary::calc_Neumann_inflow:
  {
    DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->NeumannInflow(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elevec1);
    break;
  }
  case Fluid3Boundary::calc_surface_tension:
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

    DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->ElementSurfaceTension(
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
  case Fluid3Boundary::center_of_mass_calc:
  {
    // evaluate center of mass
    if(this->Owner() == discretization.Comm().MyPID())
      DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->CenterOfMassCalculation(
          this,
          params,
          discretization,
          lm);
    break;
  }
  case Fluid3Boundary::traction_velocity_component:
  {
    DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->CalcTractionVelocityComponent(
        this,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case Fluid3Boundary::traction_Uv_integral_component:
  {
    DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->ComputeNeumannUvIntegral(
        this,
        params,
        discretization,
        lm,
        elevec1);
    break;
  }
  case Fluid3Boundary::no_penetration:
  {
    DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->NoPenetration(
        this,
        params,
        discretization,
        lm,
        elemat1,
        elemat2,
        elevec1);
    break;
  }
  case Fluid3Boundary::adjoint_neumann:
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
    dserror("Unknown type of action for Fluid3_Boundary: %s",action.c_str());
  } // end of switch(act)

  return 0;
}


/*----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition       gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Fluid3Boundary::EvaluateNeumann(
    ParameterList&            params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  return DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->EvaluateNeumann(
      this,
      params,
      discretization,
      condition,
      lm,
      elevec1,
      elemat1);
}


