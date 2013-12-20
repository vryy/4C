/*----------------------------------------------------------------------*/
/*!
\file meshfree_fluid_cell_boundary_evaluate.cpp

\brief Evaluate boundary conditions for meshfree fluid problems

<pre>
Maintainer: Keijo Nissen
            nissen@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15253
</pre>

*---------------------------------------------------------------------------*/

#include "meshfree_fluid_cell.H"
#include "meshfree_fluid_cell_boundary_interface.H"

#include "../drt_fluid_ele/fluid_ele_boundary_factory.H"
#include "../drt_fluid_ele/fluid_ele_action.H"

#include "../drt_inpar/inpar_fluid.H"

/*----------------------------------------------------------------------*
 |  evaluate the meshfree boundary element (public)           nis Nov13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::MeshfreeFluidBoundary::Evaluate(
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

  switch(act)
  {
  default:
  {
    dserror("Unknown type of action '%i' for Meshfree Fluid Boundary", act);
    break;
  }
  } // end of switch(act)

  return 0;
}

/*----------------------------------------------------------------------*
 |  Integrate a surface/line Neumann boundary condition       nis Nov13 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::MeshfreeFluidBoundary::EvaluateNeumann(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  return DRT::ELEMENTS::FluidBoundaryFactory::ProvideImplMeshfree(Shape(),"std_meshfree")->EvaluateNeumann(
      this,
      params,
      discretization,
      condition,
      lm,
      elevec1,
      elemat1);
}
