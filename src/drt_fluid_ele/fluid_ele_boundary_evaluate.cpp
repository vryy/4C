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
  return DRT::ELEMENTS::Fluid3BoundaryImplInterface::Impl(this)->Evaluate(
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
      elemat1
      );
}


