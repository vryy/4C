/*----------------------------------------------------------------------*/
/*! \file
\brief element evaluation routines

\level 1

*/

/*----------------------------------------------------------------------*
 | headers                                                    gjb 01/08 |
 *----------------------------------------------------------------------*/
#include "baci_thermo_ele_impl.hpp"
#include "baci_thermo_element.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | evaluate the element for volume coupling (public)         dano 02/10 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Thermo::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Thermo element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return DRT::ELEMENTS::TemperImplInterface::Impl(this)->Evaluate(
      this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
}  // Evaluate


/*----------------------------------------------------------------------*
 | do nothing (public)                                       dano 09/09 |
 |                                                                      |
 | The function is just a dummy. For the thermo elements, the           |
 | integration of the volume neumann (body forces) loads takes place    |
 | in the element. We need it there for the stabilisation terms!        |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Thermo::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  return DRT::ELEMENTS::TemperImplInterface::Impl(this)->EvaluateNeumann(
      this, params, discretization, lm, elevec1, elemat1);
}

FOUR_C_NAMESPACE_CLOSE
