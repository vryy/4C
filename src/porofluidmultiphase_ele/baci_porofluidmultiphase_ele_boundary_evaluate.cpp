/*----------------------------------------------------------------------*/
/*! \file
 \brief evaluation methods of the porofluidmultiphase boundary element

   \level 3

 *----------------------------------------------------------------------*/

#include "baci_lib_discret.hpp"
#include "baci_lib_element.hpp"
#include "baci_porofluidmultiphase_ele.hpp"
#include "baci_porofluidmultiphase_ele_action.hpp"
#include "baci_porofluidmultiphase_ele_boundary_factory.hpp"
#include "baci_porofluidmultiphase_ele_interface.hpp"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             vuong 08/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  dserror("not implemented. Use the Evaluate() method with Location Array instead!");
  return -1;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             vuong 08/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  // we assume here, that numdofpernode is equal for every node within
  // the element and does not change during the computations
  const int numdofpernode = NumDofPerNode(*(Nodes()[0]));

  // copy pointers to matrices and vectors into std::vector
  std::vector<CORE::LINALG::SerialDenseMatrix*> elemat(2);
  elemat[0] = &elemat1;
  elemat[1] = &elemat2;
  std::vector<CORE::LINALG::SerialDenseVector*> elevec(3);
  elevec[0] = &elevec1;
  elevec[1] = &elevec2;
  elevec[2] = &elevec3;

  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Transport
  // boundary element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return DRT::ELEMENTS::PoroFluidMultiPhaseBoundaryFactory::ProvideImpl(
      this, numdofpernode, discretization.Name())
      ->Evaluate(this, params, discretization, la, elemat, elevec);
}


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary condition on boundary element   vuong 08/16 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::PoroFluidMultiPhaseBoundary::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  // add Neumann boundary condition to parameter list
  params.set<DRT::Condition*>("condition", &condition);

  // build the location array
  LocationArray la(discretization.NumDofSets());
  DRT::Element::LocationVector(discretization, la, false);

  // evaluate boundary element
  return Evaluate(params, discretization, la, *elemat1, *elemat1, elevec1, elevec1, elevec1);
}

BACI_NAMESPACE_CLOSE
