/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluation of a scatra element that does not contain any physics

\level 2


*/
/*----------------------------------------------------------------------*/

#include "baci_scatra_ele_calc.H"
#include "baci_scatra_ele.H"
#include "baci_scatra_ele_action.H"
#include "baci_scatra_ele_calc_no_physics.H"

/*----------------------------------------------------------------------*
 | evaluate action                                        gebauer 06/19 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcNoPhysics<distype, probdim>::EvaluateAction(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1_epetra,
    CORE::LINALG::SerialDenseMatrix& elemat2_epetra,
    CORE::LINALG::SerialDenseVector& elevec1_epetra,
    CORE::LINALG::SerialDenseVector& elevec2_epetra,
    CORE::LINALG::SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::Action::time_update_material:
      break;

    default:
    {
      dserror(
          "Not acting on this action. This ImplType is designed to be a dummy element without "
          "any physics used in SSI simulations. At the moment, only the minimal set of actions "
          "are implemented, that are needed for simulating a one-way coupling from scatra to "
          "the structure, while reading the scatra results from a result file.");
      break;
    }
  }  // switch(action)

  return 0;
}

// include forward declaration of template classes
#include "baci_scatra_ele_calc_no_physics_fwd.hpp"
