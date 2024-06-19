/*----------------------------------------------------------------------*/
/*! \file
\brief

Evaluate boundary conditions for thermo problems

\level 1

*/

/*----------------------------------------------------------------------*
 | headers                                                   dano 09/09 |
 *----------------------------------------------------------------------*/
#include "4C_thermo_ele_boundary_impl.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | evaluate the element for volume coupling (public)         dano 02/10 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::ThermoBoundary::evaluate(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Elements::Element::LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Thermo
  // boundary element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return Discret::ELEMENTS::TemperBoundaryImplInterface::Impl(this)->evaluate(
      this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
}  // Evaluate in case of multiple dofsets


/*----------------------------------------------------------------------*
 | integrate a Surface/Line Neumann boundary condition       dano 09/09 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::ThermoBoundary::evaluate_neumann(Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Thermo
  // boundary element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return Discret::ELEMENTS::TemperBoundaryImplInterface::Impl(this)->evaluate_neumann(
      this, params, discretization, condition, lm, elevec1);
}

FOUR_C_NAMESPACE_CLOSE
