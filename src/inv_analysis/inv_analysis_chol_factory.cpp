/*----------------------------------------------------------------------*/
/*! \file
\brief Factory for the Cholesky factorization

\level 3

*/
/*----------------------------------------------------------------------*/
#include "inv_analysis_chol_factory.H"

#include "inv_analysis_chol_factor.H"
#include "utils_exceptions.H"

#include <Epetra_Map.h>
#include <Epetra_CrsMatrix.h>


Teuchos::RCP<INVANA::CholFactorBase> INVANA::CholFactory::Create(
    Teuchos::RCP<Epetra_CrsMatrix> A, Teuchos::ParameterList& p)
{
  Teuchos::RCP<CholFactorBase> factor = Teuchos::null;

  // initialize cholesky factor
  factor = Teuchos::rcp(new CholFactor(A));
  factor->SetParameters(p);
  factor->Initialize();

  // compute factorization
  int err = factor->Compute();
  if (err != 0) dserror("Factorization computation delivered: %d", err);
  factor->Print(std::cout);

  return factor;
}
