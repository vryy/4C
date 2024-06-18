/*----------------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the ALE time integration

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "4C_adapter_ale_wrapper.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::AleNOXCorrectionWrapper::prepare_time_step()
{
  AleWrapper::prepare_time_step();

  if (stepinc_ != Teuchos::null) stepinc_->PutScalar(0.0);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void Adapter::AleNOXCorrectionWrapper::evaluate(Teuchos::RCP<const Epetra_Vector> stepinc)
{
  if (stepinc != Teuchos::null)
  {
    // iteration increments
    Teuchos::RCP<Epetra_Vector> iterinc = Teuchos::rcp(new Epetra_Vector(*stepinc));
    if (stepinc_ != Teuchos::null)
    {
      iterinc->Update(-1.0, *stepinc_, 1.0);

      // update incremental displacement member to provided step increments
      // shortly: disinc_^<i> := disp^<i+1>
      stepinc_->Update(1.0, *stepinc, 0.0);
    }
    else
    {
      stepinc_ = Teuchos::rcp(new Epetra_Vector(*stepinc));
    }

    // do structural update with provided residual displacements - iteration increment
    AleWrapper::evaluate(iterinc);
  }
  else
  {
    AleWrapper::evaluate(Teuchos::null);
  }

  return;
}

FOUR_C_NAMESPACE_CLOSE
