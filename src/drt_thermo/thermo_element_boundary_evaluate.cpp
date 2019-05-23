/*----------------------------------------------------------------------*/
/*!
\brief

Evaluate boundary conditions for thermo problems

\level 1

\maintainer Christoph Meier
*/

/*----------------------------------------------------------------------*
 | definitions                                               dano 09/09 |
 *----------------------------------------------------------------------*/
#ifdef D_THERMO


/*----------------------------------------------------------------------*
 | headers                                                   dano 09/09 |
 *----------------------------------------------------------------------*/
#include "thermo_ele_boundary_impl.H"


/*----------------------------------------------------------------------*
 | evaluate the element for volume coupling (public)         dano 02/10 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ThermoBoundary::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1, Epetra_SerialDenseMatrix& elemat2,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseVector& elevec2,
    Epetra_SerialDenseVector& elevec3)
{
  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Thermo
  // boundary element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return DRT::ELEMENTS::TemperBoundaryImplInterface::Impl(this)->Evaluate(
      this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
}  // Evaluate in case of multiple dofsets


/*----------------------------------------------------------------------*
 | integrate a Surface/Line Neumann boundary condition       dano 09/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ThermoBoundary::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Thermo
  // boundary element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return DRT::ELEMENTS::TemperBoundaryImplInterface::Impl(this)->EvaluateNeumann(
      this, params, discretization, condition, lm, elevec1);
}


/*----------------------------------------------------------------------*/
#endif  // #ifdef D_THERMO
