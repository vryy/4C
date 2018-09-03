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
#include "fluid_ele_boundary_factory.H"
#include "fluid_ele_boundary_parent_calc.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_opti/topopt_fluidAdjoint3_boundary.H"


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidBoundary::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::BoundaryAction act = DRT::INPUT::get<FLD::BoundaryAction>(params, "action");

  switch (act)
  {
    case FLD::integrate_Shapefunction:
    case FLD::calc_area:
    case FLD::calc_flowrate:
    case FLD::flowratederiv:
    case FLD::Outletimpedance:
    case FLD::dQdu:
    case FLD::ba_calc_node_normal:
    case FLD::calc_node_curvature:
    case FLD::calc_surface_tension:
    case FLD::conservative_outflow_bc:
    case FLD::calc_Neumann_inflow:
    case FLD::calc_pressure_bou_int:
    case FLD::center_of_mass_calc:
    case FLD::traction_velocity_component:
    case FLD::traction_Uv_integral_component:
    {
      DRT::ELEMENTS::FluidBoundaryFactory::ProvideImpl(Shape(), "std")
          ->EvaluateAction(
              this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    case FLD::fpsi_coupling:
    {
      // Note: the FluidEleBoundaryCalcPoro class is called here, as the method FPSICoupling() is
      // implemented there.
      // Todo: One could think about splitting this method in pure fluid and poro fluid part...
      // vuong 11/13
      DRT::ELEMENTS::FluidBoundaryFactory::ProvideImpl(Shape(), "poro")
          ->EvaluateAction(
              this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    case FLD::enforce_weak_dbc:
    {
      DRT::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->EvaluateWeakDBC(
          this, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::estimate_Nitsche_trace_maxeigenvalue_:
    {
      DRT::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->EstimateNitscheTraceMaxEigenvalue(
          this, params, discretization, lm, elemat1, elemat2);
      break;
    }
    case FLD::mixed_hybrid_dbc:
    {
      DRT::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->MixHybDirichlet(
          this, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::flow_dep_pressure_bc:
    {
      DRT::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->FlowDepPressureBC(
          this, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::slip_supp_bc:
    {
      DRT::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->SlipSuppBC(
          this, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::navier_slip_bc:
    {
      DRT::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->NavierSlipBC(
          this, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::ba_calc_adjoint_neumann:
    {
      DRT::ELEMENTS::FluidAdjoint3BoundaryImplInterface::Impl(this)->EvaluateNeumann(
          this, params, discretization, lm, elevec1);
      break;
    }
    default:
    {
      dserror("Unknown type of action '%i' for Fluid Boundary", act);
      break;
    }
  }  // end of switch(act)

  return 0;
}

/*----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition       gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidBoundary::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  return DRT::ELEMENTS::FluidBoundaryFactory::ProvideImpl(Shape(), "std")
      ->EvaluateNeumann(this, params, discretization, condition, lm, elevec1, elemat1);
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::FluidBoundary::LocationVector(const Discretization& dis, LocationArray& la,
    bool doDirichlet, const std::string& condstring, Teuchos::ParameterList& params) const
{
  // get the action required
  const FLD::BoundaryAction act = DRT::INPUT::get<FLD::BoundaryAction>(params, "action");

  switch (act)
  {
    case FLD::enforce_weak_dbc:
    case FLD::mixed_hybrid_dbc:
    case FLD::flow_dep_pressure_bc:
    case FLD::slip_supp_bc:
    case FLD::navier_slip_bc:
    case FLD::fpsi_coupling:
      // special cases: the boundary element assembles also into
      // the inner dofs of its parent element
      // note: using these actions, the element will get the parent location vector
      //       as input in the respective evaluate routines
      ParentElement()->LocationVector(dis, la, doDirichlet);
      break;
    case FLD::ba_none:
      dserror("No action supplied");
      break;
    default:
      // standard case: element assembles into its own dofs only
      DRT::Element::LocationVector(dis, la, doDirichlet);
      break;
  }
  return;
}
