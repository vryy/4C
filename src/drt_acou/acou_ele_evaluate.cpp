/*--------------------------------------------------------------------------*/
/*!
\file acou_ele_evaluate.cpp
\brief acoustic element evaluation base routines

\level 2

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15244
*/
/*--------------------------------------------------------------------------*/

#include "acou_ele_factory.H"
#include "../drt_inpar/inpar_acou.H"
#include "acou_ele_interface.H"
#include "acou_ele.H"
#include "acou_sol_ele.H"

/*---------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::AcouType::PreEvaluate(DRT::Discretization& dis, Teuchos::ParameterList& p,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix1,
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix2, Teuchos::RCP<Epetra_Vector> systemvector1,
    Teuchos::RCP<Epetra_Vector> systemvector2, Teuchos::RCP<Epetra_Vector> systemvector3)
{
  return;
}

/*---------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Acou::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  Teuchos::RCP<MAT::Material> mat = Material();

  if (dynamic_cast<const DRT::ELEMENTS::AcouSol*>(this))
    return DRT::ELEMENTS::AcouFactory::ProvideImpl(Shape(), "sol")
        ->Evaluate(
            this, discretization, lm, params, mat, elemat1, elemat2, elevec1, elevec2, elevec3);
  else
    return DRT::ELEMENTS::AcouFactory::ProvideImpl(Shape(), "std")
        ->Evaluate(
            this, discretization, lm, params, mat, elemat1, elemat2, elevec1, elevec2, elevec3);

  return -1;
}

/*---------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Acou::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}
