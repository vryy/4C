/*----------------------------------------------------------------------*/
/*! \file
\brief Wrapper for the scatra time integrator.
\level 1
 */
/*----------------------------------------------------------------------*/


#include "4C_adapter_scatra_wrapper.hpp"

#include "4C_linalg_sparseoperator.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
ADAPTER::AdapterScatraWrapper::AdapterScatraWrapper(Teuchos::RCP<ScatraInterface> scatra)
    : scatra_timint_(scatra)
{
  Teuchos::rcp_dynamic_cast<SCATRA::ScaTraTimIntImpl>(scatra_timint_, true)
      ->set_model_evaluatro_ptr(this);
}


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ADAPTER::AdapterScatraWrapper::evaluate_additional_solution_depending_models(
    Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix, Teuchos::RCP<Epetra_Vector> rhs)
{
  // do nothing so far
  return;
}

FOUR_C_NAMESPACE_CLOSE
