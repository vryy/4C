/*!----------------------------------------------------------------------
\file scatra_element_boundary_evaluate.cpp
\brief

Evaluate boundary conditions for scalar transport problems

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>
 *----------------------------------------------------------------------*/
#if defined(D_FLUID2) || defined(D_FLUID3)
#ifdef CCADISCRET

#include "scatra_ele_boundary_impl.H"

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::Evaluate(
    ParameterList&            params,
    DRT::Discretization&      discretization,
    vector<int>&              lm,
    Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Transport
  // boundary element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return DRT::ELEMENTS::ScaTraBoundaryImplInterface::Impl(this)->Evaluate(
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
    ParameterList&            params,
    DRT::Discretization&      discretization,
    DRT::Condition&           condition,
    vector<int>&              lm,
    Epetra_SerialDenseVector& elevec1)
{
  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Transport
  // boundary element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return DRT::ELEMENTS::ScaTraBoundaryImplInterface::Impl(this)->EvaluateNeumann(
      this,
      params,
      discretization,
      condition,
      lm,
      elevec1
      );
}


#endif  // #ifdef CCADISCRET
#endif  // #if defined(D_FLUID2) || defined(D_FLUID3)
