/*----------------------------------------------------------------------*/
/*! \file
\brief Wrapper for the scatra time integrator.
\level 1
 */
/*----------------------------------------------------------------------*/


#include "adapter_scatra_wrapper.H"
#include "scatra_timint_implicit.H"
#include "lib_dserror.H"
#include "linalg_sparseoperator.H"


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
ADAPTER::AdapterScatraWrapper::AdapterScatraWrapper(Teuchos::RCP<ScatraInterface> scatra)
    : scatra_timint_(scatra)
{
  Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntImpl>(scatra_timint_, true)
      ->SetModelEvaluatroPtr(this);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::AdapterScatraWrapper::EvaluateAdditionalSolutionDependingModels(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix, Teuchos::RCP<Epetra_Vector> rhs)
{
  // do nothing so far
  return;
}
