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
  if (!material_post_setup_)
  {
    DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->MaterialPostSetup(*this);
    material_post_setup_ = true;
  }

  // get ptr to interface to time integration
  SetParamsInterfacePtr(params);

  // get the calculation class
  DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->Evaluate(
      this, params, discretization, lm, elemat1, elemat2, elevec1, elevec2, elevec3);

  return 0;
}
int DRT::ELEMENTS::Solid::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  return DRT::ELEMENTS::SolidFactory::ProvideImpl(this)->EvaluateNeumann(
      this, params, discretization, condition, lm, nullptr, elevec1, elemat1);
}
