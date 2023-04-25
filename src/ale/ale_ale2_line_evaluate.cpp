/*----------------------------------------------------------------------------*/
/*! \file

\brief Evaluate ALE2 elements

\level 1

*/
/*----------------------------------------------------------------------------*/

#include "ale_ale2.H"
#include "lib_discret.H"
#include "lib_dserror.H"

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
int DRT::ELEMENTS::Ale2Line::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  dserror("Line Neumann condition not implemented for Ale2");

  return 0;
}
