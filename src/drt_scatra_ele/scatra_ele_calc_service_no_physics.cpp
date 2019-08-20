/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluation of a scatra element that does not contain any physics

\level 2

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_calc.H"
#include "scatra_ele.H"
#include "scatra_ele_action.H"
#include "scatra_ele_calc_no_physics.H"

/*----------------------------------------------------------------------*
 | evaluate action                                        gebauer 06/19 |
 *----------------------------------------------------------------------*/
template <DRT::Element::DiscretizationType distype, int probdim>
int DRT::ELEMENTS::ScaTraEleCalcNoPhysics<distype, probdim>::EvaluateAction(DRT::Element* ele,
    Teuchos::ParameterList& params, DRT::Discretization& discretization,
    const SCATRA::Action& action, DRT::Element::LocationArray& la,
    Epetra_SerialDenseMatrix& elemat1_epetra, Epetra_SerialDenseMatrix& elemat2_epetra,
    Epetra_SerialDenseVector& elevec1_epetra, Epetra_SerialDenseVector& elevec2_epetra,
    Epetra_SerialDenseVector& elevec3_epetra)
{
  // determine and evaluate action
  switch (action)
  {
    case SCATRA::time_update_material:
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
#include "scatra_ele_calc_no_physics_fwd.hpp"