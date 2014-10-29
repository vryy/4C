/*----------------------------------------------------------------------*/
/*!
\file ad_str_invana.cpp

\brief Wrapper for the structural time integration which gives fine grained
       access in the time loop

<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-10361
</pre>
*/
/*----------------------------------------------------------------------*/

#include "ad_str_invana.H"
#include "../drt_lib/drt_utils_timintmstep.H"
#include "../drt_inpar/inpar_structure.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureInvana::StructureInvana(Teuchos::RCP<Structure> structure):
  StructureTimeLoop(structure),
stime_(Teuchos::null),
sdis_(Teuchos::null),
svel_(Teuchos::null),
singlesteponly_(false)
{
  stime_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<double>(0, 0, 0.0));
  // initialize to zero!! ResizeStorage only initializes the new vectors
  sdis_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
  svel_ = Teuchos::rcp(new DRT::UTILS::TimIntMStep<Epetra_Vector>(0, 0, DofRowMapView(), true));
}

/*----------------------------------------------------------------------*/
/* Resizing of multi-step quantities */
void ADAPTER::StructureInvana::ResizeStorage()
{
  // resize storage
  stime_->Resize(-(Step()-1), 0, 0.0);
  // -> time is already stored for the next time step
  sdis_->Resize(-(Step()-1), 0, DofRowMapView(), true);
  svel_->Resize(-(Step()-1), 0, DofRowMapView(), true);

  return;
}

void ADAPTER::StructureInvana::SetTimeStepStateOld(double time, int step,
    Teuchos::RCP<Epetra_Vector> disp,
    Teuchos::RCP<Epetra_Vector> vel)
{
  singlesteponly_=true;
  SetTimen(time);
  SetStepn(step);
  WriteAccessDispnp()->Update(1.0,*disp,0.0);
  WriteAccessVelnp()->Update(1.0,*vel,0.0);

  // to stop after this step:
  double timemax=time+Dt();
  SetTimeEnd(timemax);
}

void ADAPTER::StructureInvana::PrePredict()
{
  if (singlesteponly_)
  {
    // do what was not done in the update of the previous
    // loop during a standard integrate()-call since we only
    // do this timestep
    Update();
  }

  return;
}

void ADAPTER::StructureInvana::PostPredict()
{
  ResizeStorage();

  return;
}

void ADAPTER::StructureInvana::PostUpdate()
{
  // state and time was already updated so Dispn() and Veln()
  // (pointing to "(*state_)(0)") give the results of
  // this time step
  sdis_->UpdateSteps(*Dispn());
  svel_->UpdateSteps(*Veln());
  stime_->UpdateSteps(TimeOld());

  return;
}

