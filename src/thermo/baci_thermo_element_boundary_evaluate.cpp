/*----------------------------------------------------------------------*/
/*! \file
\brief

Evaluate boundary conditions for thermo problems

\level 1

*/

/*----------------------------------------------------------------------*
 | headers                                                   dano 09/09 |
 *----------------------------------------------------------------------*/
#include "baci_thermo_ele_boundary_impl.hpp"

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | evaluate the element for volume coupling (public)         dano 02/10 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::ThermoBoundary::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
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
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
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

BACI_NAMESPACE_CLOSE
