/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluation of a scatra element that does not contain any physics

\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_calc.hpp"
#include "4C_scatra_ele_calc_no_physics.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | evaluate action                                        gebauer 06/19 |
 *----------------------------------------------------------------------*/
template <Core::FE::CellType distype, int probdim>
int Discret::ELEMENTS::ScaTraEleCalcNoPhysics<distype, probdim>::evaluate_action(
    Core::Elements::Element* ele, Teuchos::ParameterList& params,
    Core::FE::Discretization& discretization, const ScaTra::Action& action,
    Core::Elements::Element::LocationArray& la, Core::LinAlg::SerialDenseMatrix& elemat1_epetra,
    Core::LinAlg::SerialDenseMatrix& elemat2_epetra,
    Core::LinAlg::SerialDenseVector& elevec1_epetra,
    Core::LinAlg::SerialDenseVector& elevec2_epetra,
    Core::LinAlg::SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case ScaTra::Action::time_update_material:
      break;

    default:
    {
      FOUR_C_THROW(
          "Not acting on this action. This ImplType is designed to be a dummy element without "
          "any physics used in SSI simulations. At the moment, only the minimal set of actions "
          "are implemented, that are needed for simulating a one-way coupling from scatra to "
          "the structure, while reading the scatra results from a result file.");
      break;
    }
  }  // switch(action)

  return 0;
}

FOUR_C_NAMESPACE_CLOSE

// include forward declaration of template classes
#include "4C_scatra_ele_calc_no_physics_fwd.hpp"
