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

  // switch between different physical types as used below
  std::string impltype = "std";
  switch(params.get<int>("physical type",INPAR::FLUID::incompressible))
  {
  case INPAR::FLUID::loma:    impltype = "std";  break;
  case INPAR::FLUID::poro:    impltype = "poro"; break;
  case INPAR::FLUID::poro_p1: impltype = "poro"; break;
  case INPAR::FLUID::poro_p2: impltype = "poro"; break;
  }

  switch(act)
  {
  case FLD::integrate_Shapefunction:
  case FLD::calc_area:
  case FLD::calc_flowrate:
  case FLD::flowratederiv:
  case FLD::Outletimpedance:
  case FLD::ba_calc_node_normal:
  case FLD::calc_node_curvature:
  case FLD::calc_surface_tension:
  case FLD::conservative_outflow_bc:
  case FLD::calc_Neumann_inflow:
  case FLD::calc_pressure_bou_int:
  case FLD::center_of_mass_calc:
  case FLD::traction_velocity_component:
  case FLD::traction_Uv_integral_component:
  case FLD::no_penetration:
  case FLD::no_penetrationIDs:
  case FLD::poro_boundary:
  case FLD::poro_prescoupl:
  case FLD::fpsi_coupling:
  {
  DRT::ELEMENTS::FluidBoundaryFactory::ProvideImpl(Shape(),impltype)->EvaluateAction(
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
  {
    dserror("Unknown type of action for Fluid_Boundary!");
    break;
  }
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
  return DRT::ELEMENTS::FluidBoundaryFactory::ProvideImpl(Shape(),"std")->EvaluateNeumann(
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
    //
    // todo:
    // the following hacky request for the discretization name is temporary!
    // it will be replaced as soon as a specialized poro boundary element is
    // implemented ! // rauch 11/13
    //
  case FLD::calc_flowrate:
  {
    if(dis.Name() == "fluid")
    {
      DRT::Element::LocationVector(dis,la,doDirichlet);
    }
    else if (dis.Name() == "porofluid")
    {
      // special cases: the boundary element assembles also into
      // the inner dofs of its parent element
      // note: using these actions, the element will get the parent location vector
      //       as input in the respective evaluate routines
      parent_->LocationVector(dis,la,doDirichlet);
    }
    else
      dserror("action 'calc_flowrate' not implemented for this physical type");

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

