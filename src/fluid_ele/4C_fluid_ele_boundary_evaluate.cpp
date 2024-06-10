/*----------------------------------------------------------------------*/
/*! \file

\brief outer evaluation class of fluid terms at integration points of boundaries

\level 1


*/
/*----------------------------------------------------------------------*/


#include "4C_fem_discretization.hpp"
#include "4C_fluid_ele.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_ele_boundary_calc.hpp"
#include "4C_fluid_ele_boundary_factory.hpp"
#include "4C_fluid_ele_boundary_parent_calc.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::FluidBoundary::Evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, std::vector<int>& lm,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // get the action required
  const FLD::BoundaryAction act = Core::UTILS::GetAsEnum<FLD::BoundaryAction>(params, "action");

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
      Discret::ELEMENTS::FluidBoundaryFactory::ProvideImpl(Shape(), "std")
          ->evaluate_action(
              this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    case FLD::fpsi_coupling:
    {
      // Note: the FluidEleBoundaryCalcPoro class is called here, as the method fpsi_coupling() is
      // implemented there.
      // Todo: One could think about splitting this method in pure fluid and poro fluid part...
      // vuong 11/13
      Discret::ELEMENTS::FluidBoundaryFactory::ProvideImpl(Shape(), "poro")
          ->evaluate_action(
              this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);
      break;
    }
    case FLD::enforce_weak_dbc:
    {
      Discret::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->EvaluateWeakDBC(
          this, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::estimate_Nitsche_trace_maxeigenvalue_:
    {
      Discret::ELEMENTS::FluidBoundaryParentInterface::Impl(this)
          ->estimate_nitsche_trace_max_eigenvalue(
              this, params, discretization, lm, elemat1, elemat2);
      break;
    }
    case FLD::mixed_hybrid_dbc:
    {
      Discret::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->MixHybDirichlet(
          this, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::flow_dep_pressure_bc:
    {
      Discret::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->FlowDepPressureBC(
          this, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::slip_supp_bc:
    {
      Discret::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->SlipSuppBC(
          this, params, discretization, lm, elemat1, elevec1);
      break;
    }
    case FLD::navier_slip_bc:
    {
      Discret::ELEMENTS::FluidBoundaryParentInterface::Impl(this)->NavierSlipBC(
          this, params, discretization, lm, elemat1, elevec1);
      break;
    }
    default:
    {
      FOUR_C_THROW("Unknown type of action '%i' for Fluid Boundary", act);
      break;
    }
  }  // end of switch(act)

  return 0;
}

/*----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition       gjb 01/09 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::FluidBoundary::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  return Discret::ELEMENTS::FluidBoundaryFactory::ProvideImpl(Shape(), "std")
      ->evaluate_neumann(this, params, discretization, condition, lm, elevec1, elemat1);
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 *----------------------------------------------------------------------*/
void Discret::ELEMENTS::FluidBoundary::LocationVector(const Core::FE::Discretization& dis,
    LocationArray& la, bool doDirichlet, const std::string& condstring,
    Teuchos::ParameterList& params) const
{
  // get the action required
  const FLD::BoundaryAction act = Core::UTILS::GetAsEnum<FLD::BoundaryAction>(params, "action");

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
      parent_element()->LocationVector(dis, la, doDirichlet);
      break;
    case FLD::ba_none:
      FOUR_C_THROW("No action supplied");
      break;
    default:
      // standard case: element assembles into its own dofs only
      Core::Elements::Element::LocationVector(dis, la, doDirichlet);
      break;
  }
  return;
}

FOUR_C_NAMESPACE_CLOSE
