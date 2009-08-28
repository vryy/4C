
#ifdef CCADISCRET

#include "adapter_structure_wrapper.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureNOXCorrectionWrapper::PrepareTimeStep()
{
  StructureWrapper::PrepareTimeStep();
  if (disinc_!=Teuchos::null)
    disinc_->PutScalar(0.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ADAPTER::StructureNOXCorrectionWrapper::Evaluate(Teuchos::RCP<const Epetra_Vector> disp)
{
  // Yes, this is complicated. But we have to be very careful
  // here. The field solver always expects an increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest
  // increment only.
  if (disp!=Teuchos::null)
  {
    // residual displacements (or iteration increments or iteratively
    // incremental displacements)
    Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(new Epetra_Vector(*disp));
    if (disinc_!=Teuchos::null)
    {
      disi->Update(-1.0,*disinc_,1.0);

      // update incremental displacement member to provided step increments
      // shortly: disinc_^<i> := disp^<i+1>
      disinc_->Update(1.0,*disp,0.0);
    }
    else
    {
      disinc_ = Teuchos::rcp(new Epetra_Vector(*disp));
    }

    // do structural update with provided residual displacements
    StructureWrapper::Evaluate(disi);
  }
  else
  {
    StructureWrapper::Evaluate(Teuchos::null);
  }
}

#endif
