/*----------------------------------------------------------------------*/
/*! \file
\brief Wrapper for the scatra time integrator.
\level 1
 */
/*----------------------------------------------------------------------*/


#include "baci_adapter_scatra_wrapper.H"

#include "baci_linalg_sparseoperator.H"
#include "baci_scatra_timint_implicit.H"
#include "baci_utils_exceptions.H"

BACI_NAMESPACE_OPEN


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
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix, Teuchos::RCP<Epetra_Vector> rhs)
{
  // do nothing so far
  return;
}

BACI_NAMESPACE_CLOSE
