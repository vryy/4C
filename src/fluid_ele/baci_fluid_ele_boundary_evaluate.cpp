/*----------------------------------------------------------------------*/
/*! \file

\brief outer evaluation class of fluid terms at integration points of boundaries

\level 1


*/
/*----------------------------------------------------------------------*/


#include "baci_fluid_ele.hpp"
#include "baci_fluid_ele_action.hpp"
#include "baci_fluid_ele_boundary_calc.hpp"
#include "baci_fluid_ele_boundary_factory.hpp"
#include "baci_fluid_ele_boundary_parent_calc.hpp"
#include "baci_lib_discret.hpp"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::FluidBoundary::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::BoundaryAction act = CORE::UTILS::GetAsEnum<FLD::BoundaryAction>(params, "action");

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
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
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
  const FLD::BoundaryAction act = CORE::UTILS::GetAsEnum<FLD::BoundaryAction>(params, "action");

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

BACI_NAMESPACE_CLOSE
