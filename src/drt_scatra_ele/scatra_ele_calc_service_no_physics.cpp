/*----------------------------------------------------------------------*/
/*!

\brief Evaluation of a scatra element that does not contain any physics. Currently only implements
 the minimal set of actions needed for reading the scatra results from a restart file and simulating
 a one-way coupling to the structure. This ImplType is currently not capable to be used in solving
 the scatra equations, as the needed actions are not implemented yet.

\level 2

\maintainer Amadeus Gebauer

*/
/*----------------------------------------------------------------------*/

#include "scatra_ele_calc.H"
#include "scatra_ele.H"
#include "scatra_ele_action.H"
#include "scatra_ele_calc_no_physics.h"

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

// template classes

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
// template classes
// 1D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line2, 1>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line2, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line2, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::line3, 1>;

// 2D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tri3, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tri3, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tri6, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::quad4, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::quad4, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::quad8>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::quad9, 2>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::nurbs9, 2>;

// 3D elements
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::hex8, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::hex20>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::hex27, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tet4, 3>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::tet10, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::wedge6>;
template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::pyramid5, 3>;
// template class DRT::ELEMENTS::ScaTraEleCalcNoPhysics<DRT::Element::nurbs27>;