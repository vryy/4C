/*----------------------------------------------------------------------*/
/*!

\brief Wrapper for the structural time integration

\level 1

\maintainer Anh-Tu Vuong

*/
/*----------------------------------------------------------------------*/

#include "ad_str_wrapper.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureNOXCorrectionWrapper::PrepareTimeStep()
{
  StructureWrapper::PrepareTimeStep();
  if (disstepinc_ != Teuchos::null) disstepinc_->PutScalar(0.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureNOXCorrectionWrapper::Evaluate(Teuchos::RCP<const Epetra_Vector> disstepinc)
{
  // The field solver always expects an iteration increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest iteration
  // increment only.
  // Naming:
  //
  // x^n+1_i+1 = x^n+1_i + disiterinc  (sometimes referred to as residual increment), and
  //
  // x^n+1_i+1 = x^n     + disstepinc

  if (disstepinc != Teuchos::null)
  {
    // iteration increments
    Teuchos::RCP<Epetra_Vector> disiterinc = Teuchos::rcp(new Epetra_Vector(*disstepinc));
    if (disstepinc_ != Teuchos::null)
    {
      disiterinc->Update(-1.0, *disstepinc_, 1.0);

      // update incremental displacement member to provided step increments
      // shortly: disinc_^<i> := disp^<i+1>
      disstepinc_->Update(1.0, *disstepinc, 0.0);
    }
    else
    {
      disstepinc_ = Teuchos::rcp(new Epetra_Vector(*disstepinc));
    }

    // do structural update with provided residual displacements - iteration increment
    StructureWrapper::Evaluate(disiterinc);
  }
  else
  {
    StructureWrapper::Evaluate(Teuchos::null);
  }
}
