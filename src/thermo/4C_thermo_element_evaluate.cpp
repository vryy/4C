/*----------------------------------------------------------------------*/
/*! \file
\brief element evaluation routines

\level 1

*/

/*----------------------------------------------------------------------*
 | headers                                                    gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "4C_thermo_ele_impl.hpp"
#include "4C_thermo_element.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | evaluate the element for volume coupling (public)         dano 02/10 |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::Thermo::Evaluate(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, LocationArray& la,
    Core::LinAlg::SerialDenseMatrix& elemat1, Core::LinAlg::SerialDenseMatrix& elemat2,
    Core::LinAlg::SerialDenseVector& elevec1, Core::LinAlg::SerialDenseVector& elevec2,
    Core::LinAlg::SerialDenseVector& elevec3)
{
  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Thermo element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return Discret::ELEMENTS::TemperImplInterface::Impl(this)->Evaluate(
      this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
}  // Evaluate


/*----------------------------------------------------------------------*
 | do nothing (public)                                       dano 09/09 |
 |                                                                      |
 | The function is just a dummy. For the thermo elements, the           |
 | integration of the volume neumann (body forces) loads takes place    |
 | in the element. We need it there for the stabilisation terms!        |
 *----------------------------------------------------------------------*/
int Discret::ELEMENTS::Thermo::evaluate_neumann(Teuchos::ParameterList& params,
    Discret::Discretization& discretization, Core::Conditions::Condition& condition,
    std::vector<int>& lm, Core::LinAlg::SerialDenseVector& elevec1,
    Core::LinAlg::SerialDenseMatrix* elemat1)
{
  return Discret::ELEMENTS::TemperImplInterface::Impl(this)->evaluate_neumann(
      this, params, discretization, lm, elevec1, elemat1);
}

FOUR_C_NAMESPACE_CLOSE
