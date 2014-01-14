/*!----------------------------------------------------------------------
\file scatra_element_boundary_evaluate.cpp
\brief

Evaluate boundary conditions for scalar transport problems

<pre>
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 *----------------------------------------------------------------------*/

#include "scatra_ele.H"
#include "scatra_ele_boundary_impl.H"

#include "scatra_ele_action.H"

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::Evaluate(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    std::vector<int>&         lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // the type of scalar transport problem has to be provided for all actions!
  const INPAR::SCATRA::ScaTraType scatratype = DRT::INPUT::get<INPAR::SCATRA::ScaTraType>(params, "scatratype");
  if (scatratype == INPAR::SCATRA::scatratype_undefined)
    dserror("Element parameter SCATRATYPE has not been set!");

  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Transport
  // boundary element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return DRT::ELEMENTS::ScaTraBoundaryImplInterface::Impl(this,scatratype)->Evaluate(
      this,
      params,
      discretization,
      lm,
      elemat1,
      elemat2,
      elevec1,
      elevec2,
      elevec3
      );
}


/*----------------------------------------------------------------------*
 |  Integrate a Surface/Line Neumann boundary condition       gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::EvaluateNeumann(
    Teuchos::ParameterList&   params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    std::vector<int>&         lm,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseMatrix* elemat1)
{
  // the type of scalar transport problem has to be provided for all actions!
  const INPAR::SCATRA::ScaTraType scatratype = DRT::INPUT::get<INPAR::SCATRA::ScaTraType>(params, "scatratype");
  if (scatratype == INPAR::SCATRA::scatratype_undefined)
    dserror("Element parameter SCATRATYPE has not been set!");

  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Transport
  // boundary element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return DRT::ELEMENTS::ScaTraBoundaryImplInterface::Impl(this,scatratype)->EvaluateNeumann(
      this,
      params,
      discretization,
      condition,
      lm,
      elevec1
      );
}

/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::LocationVector(
    const Discretization&   dis,
    LocationArray&          la,
    bool                    doDirichlet,
    const std::string&      condstring,
    Teuchos::ParameterList& params
    ) const
{
  // check for the action parameter
  const SCATRA::BoundaryAction action = DRT::INPUT::get<SCATRA::BoundaryAction>(params,"action");
  switch (action)
  {
  case SCATRA::bd_calc_weak_Dirichlet:
  case SCATRA::bd_reinitialize_levelset:
  case SCATRA::bd_calc_TG_outflow:
    // special cases: the boundary element assembles also into
    // the inner dofs of its parent element
    // note: using these actions, the element will get the parent location vector
    //       as input in the respective evaluate routines
    ParentElement()->LocationVector(dis,la,doDirichlet);
    break;
  default:
    DRT::Element::LocationVector(dis,la,doDirichlet);
    break;
  }
  return;
}
