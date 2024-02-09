/*----------------------------------------------------------------------*/
/*! \file
\brief Wrapper for the scatra time integrator.
\level 1
 */
/*----------------------------------------------------------------------*/


#include "baci_adapter_scatra_wrapper.hpp"

#include "baci_linalg_sparseoperator.hpp"
#include "baci_scatra_timint_implicit.hpp"
#include "baci_utils_exceptions.hpp"

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
