/*----------------------------------------------------------------------*/
/*! \file
\brief Wrapper for the scatra time integrator.
\level 1
\maintainer Anh-Tu Vuong
 */
/*----------------------------------------------------------------------*/


#include "ad_scatra_wrapper.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_lib/drt_dserror.H"
#include "../linalg/linalg_sparseoperator.H"


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
