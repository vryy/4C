/*----------------------------------------------------------------------*/
/*!
\brief Factory for the Cholesky factorization

\level 3

\maintainer Sebastian Brandstaeter
*/
/*----------------------------------------------------------------------*/
#include "chol_factory.H"

#include "chol_factor.H"
#include "../drt_lib/drt_dserror.H"

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"


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
