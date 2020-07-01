/*----------------------------------------------------------------------------*/
/*! \file

\brief Wrapper for the ALE time integration

\level 1

*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/* header inclusions */
#include "ad_ale_wrapper.H"


/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleNOXCorrectionWrapper::PrepareTimeStep()
{
  AleWrapper::PrepareTimeStep();

  if (stepinc_ != Teuchos::null) stepinc_->PutScalar(0.0);

  return;
}

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
void ADAPTER::AleNOXCorrectionWrapper::Evaluate(Teuchos::RCP<const Epetra_Vector> stepinc)
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
    AleWrapper::Evaluate(iterinc);
  }
  else
  {
    AleWrapper::Evaluate(Teuchos::null);
  }

  return;
}
