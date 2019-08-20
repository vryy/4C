/*----------------------------------------------------------------------*/
/*!

\brief Wrapper for the field time integration

\level 2

\maintainer Christoph Ager

*/
/*-----------------------------------------------------------------------*/

#include "ad_field_wrapper.H"

/*-----------------------------------------------------------------------/
| start new time step                                                    |
/-----------------------------------------------------------------------*/
void ADAPTER::FieldWrapper::PrepareTimeStep()
{
  field_->PrepareTimeStep();
  if (NOXCorrection_) ResetStepinc();
}


void ADAPTER::FieldWrapper::UpdateStateIncrementally(Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  if (NOXCorrection_) GetIterinc(disiterinc);
  field_->UpdateStateIncrementally(disiterinc);
}

/*-----------------------------------------------------------------------/
| update dofs and evaluate elements                                      |
/-----------------------------------------------------------------------*/
void ADAPTER::FieldWrapper::Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  if (NOXCorrection_) GetIterinc(disiterinc);
  field_->Evaluate(disiterinc);
}

/*-----------------------------------------------------------------------/
| update dofs and evaluate elements                                      |
/-----------------------------------------------------------------------*/
void ADAPTER::FieldWrapper::Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc, bool firstiter)
{
  if (NOXCorrection_) GetIterinc(disiterinc);
  field_->Evaluate(disiterinc, firstiter);
}

/*-----------------------------------------------------------------------/
| Reset Step Increment                                                   |
/-----------------------------------------------------------------------*/
void ADAPTER::FieldWrapper::ResetStepinc()
{
  if (stepinc_ != Teuchos::null) stepinc_->PutScalar(0.);
}

/*-----------------------------------------------------------------------/
| Get Iteration Increment from Step Increment                            |
/-----------------------------------------------------------------------*/
void ADAPTER::FieldWrapper::GetIterinc(Teuchos::RCP<const Epetra_Vector>& stepinc)
{
  // The field solver always expects an iteration increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest iteration
  // increment only.
  // Naming:
  //
  // x^n+1_i+1 = x^n+1_i + iterinc  (sometimes referred to as residual increment), and
  //
  // x^n+1_i+1 = x^n     + stepinc
  if (stepinc != Teuchos::null)
  {
    // iteration increments
    Teuchos::RCP<Epetra_Vector> iterinc = Teuchos::rcp(new Epetra_Vector(*stepinc));
    if (stepinc_ != Teuchos::null)
    {
      iterinc->Update(-1.0, *stepinc_, 1.0);

      // update incremental dof member to provided step increments
      // shortly: disinc_^<i> := disp^<i+1>
      stepinc_->Update(1.0, *stepinc, 0.0);
    }
    else
    {
      stepinc_ = Teuchos::rcp(new Epetra_Vector(*stepinc));
    }
    // output is iterinc!
    stepinc = Teuchos::rcp(new const Epetra_Vector(*iterinc));
  }
}
