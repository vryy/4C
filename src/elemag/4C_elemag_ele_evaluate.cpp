/*--------------------------------------------------------------------------*/
/*! \file
\brief electromagnetic element evaluation base routines

\level 2

*/
/*--------------------------------------------------------------------------*/

#include "4C_elemag_diff_ele.hpp"
#include "4C_elemag_ele.hpp"
#include "4C_elemag_ele_factory.hpp"
#include "4C_elemag_ele_interface.hpp"
#include "4C_inpar_elemag.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::ElemagType::PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix2,
    Teuchos::RCP<Epetra_Vector> systemvector1, Teuchos::RCP<Epetra_Vector> systemvector2,
    Teuchos::RCP<Epetra_Vector> systemvector3)
{
  return;
}

/*---------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Elemag::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  Teuchos::RCP<MAT::Material> mat = Material();

  if (dynamic_cast<const DRT::ELEMENTS::ElemagDiff*>(this))
    return DRT::ELEMENTS::ElemagFactory::ProvideImpl(Shape(), "diff")
        ->Evaluate(
            this, discretization, lm, params, mat, elemat1, elemat2, elevec1, elevec2, elevec3);
  else
    return DRT::ELEMENTS::ElemagFactory::ProvideImpl(Shape(), "std")
        ->Evaluate(
            this, discretization, lm, params, mat, elemat1, elemat2, elevec1, elevec2, elevec3);

  return -1;
}

/*---------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Elemag::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseMatrix* elemat1)
{
  return 0;
}

FOUR_C_NAMESPACE_CLOSE
