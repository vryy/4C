/*----------------------------------------------------------------------*/
/*!
\file ad_field_wrapper.cpp

\brief Wrapper for the field time integration

<pre>
Maintainer: Ager Christoph
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
</pre>
*/
/*----------------------------------------------------------------------*/

#include "ad_field_wrapper.H"

/*-----------------------------------------------------------------------/
| start new time step                                                    |
/*----------------------------------------------------------------------*/
void ADAPTER::FieldWrapper::PrepareTimeStep()
{
  field_->PrepareTimeStep();
  if (NOXCorrection_) ResetStepinc();
}

/*-----------------------------------------------------------------------/
| update dofs and evaluate elements                                      |
/*----------------------------------------------------------------------*/
void ADAPTER::FieldWrapper::Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc)
{
  if (!NOXCorrection_)
    field_->Evaluate(disiterinc);
  else
    field_->Evaluate(GetIterinc(disiterinc));
}

/*-----------------------------------------------------------------------/
| update dofs and evaluate elements                                      |
/*----------------------------------------------------------------------*/
void ADAPTER::FieldWrapper::Evaluate(Teuchos::RCP<const Epetra_Vector> disiterinc, bool firstiter)
{
  if (!NOXCorrection_)
    field_->Evaluate(disiterinc,firstiter);
  else
    field_->Evaluate(GetIterinc(disiterinc),firstiter);
}

/*-----------------------------------------------------------------------/
| Reset Step Increment                                                   |
/*----------------------------------------------------------------------*/
void ADAPTER::FieldWrapper::ResetStepinc()
{
  if (stepinc_!=Teuchos::null)
    stepinc_->PutScalar(0.);
}

/*-----------------------------------------------------------------------/
| Get Iteration Increment from Step Increment                            |
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>  ADAPTER::FieldWrapper::GetIterinc(Teuchos::RCP<const Epetra_Vector> stepinc)
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

  if (stepinc!=Teuchos::null)
  {
    // iteration increments
    Teuchos::RCP<Epetra_Vector> iterinc = Teuchos::rcp(new Epetra_Vector(*stepinc));
    if (stepinc_!=Teuchos::null)
    {
      iterinc->Update(-1.0,*stepinc_,1.0);

      // update incremental dof member to provided step increments
      // shortly: disinc_^<i> := disp^<i+1>
      stepinc_->Update(1.0,*stepinc,0.0);
    }
    else
    {
      stepinc_ = Teuchos::rcp(new Epetra_Vector(*stepinc));
    }

    // do field update with provided residual dofs - iteration increment
    return iterinc;
  }
  else
    return Teuchos::null;
}
