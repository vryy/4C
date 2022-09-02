/*----------------------------------------------------------------------*/
/*! \file

\brief A C++ wrapper for the solid element

This file contains the element-specific evaluation routines such as
Evaluate(...), EvaluateNeumann(...), etc.

\level 1
*/

#include "../drt_lib/drt_elements_paramsinterface.H"
#include "../drt_structure_new/str_elements_paramsinterface.H"
#include "../drt_lib/drt_dserror.H"
#include "solid_ele.H"
#include "solid_ele_factory.H"
#include "../drt_patspec/patspec.H"
#include "solid_ele_interface.H"

int DRT::ELEMENTS::Solid::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  // get ptr to interface to time integration
  SetParamsInterfacePtr(params);

  // check for patient specific data
  // TODO: do we really need to do this here ???
  PATSPEC::GetILTDistance(Id(), params, discretization);
  PATSPEC::GetLocalRadius(Id(), params, discretization);
  PATSPEC::GetInnerRadius(Id(), params, discretization);

  // get the calculation class
  DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->Evaluate(
      this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);

  return 0;
}
int DRT::ELEMENTS::Solid::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  dserror("stop");
  return 0;
  //  return DRT::ELEMENTS::SolidFactory::ProvideImpl(Shape(),eletech_)->EvaluateNeumann(
  //      this,params,discretization,condition,lm,NULL,elevec1,elemat1);
}
