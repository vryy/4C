/*----------------------------------------------------------------------*/
/*!

\brief Wrapper for the structural time integration which gives fine grained
       access in the time loop

\level 2

\maintainer Harald Willmann

*/
/*----------------------------------------------------------------------*/

#include "ad_str_invana.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/inpar_structure.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::StructureInvana::StructureInvana(Teuchos::RCP<Structure> structure)
    : StructureTimeLoop(structure),
      stime_(Teuchos::null),
      sdis_(Teuchos::null),
      svel_(Teuchos::null),
      singlesteponly_(false)
{
  // initialize storage
  stime_ = Teuchos::rcp(new std::map<int, double>);
  sdis_ = Teuchos::rcp(new std::map<int, Epetra_Vector>);
  svel_ = Teuchos::rcp(new std::map<int, Epetra_Vector>);
}

/*----------------------------------------------------------------------*/
/* Resizing of multi-step quantities */
void ADAPTER::StructureInvana::ResizeStorage() { return; }

void ADAPTER::StructureInvana::SetTimeStepStateOld(
    double time, int step, Teuchos::RCP<Epetra_Vector> disp, Teuchos::RCP<Epetra_Vector> vel)
{
  singlesteponly_ = true;

  // ------- TIME STEP HACKS
  // set timen_ lower to pass NotFinished()
  SetTimen(time - Dt());

  // set (*time_)[0] to have the input target parameter "time"
  // after the Adaptivity's PrepareTimeStep-hack  which is
  // called during Predict()
  SetTime(time - Dt());
  // ------- TIME STEP HACKS

  SetStepn(step);

  // here one has to set into (*dis_)(0) since in Predict
  // disn_ is updated from (*dis_)(0) depending on predictortype
  WriteAccessDispn()->Update(1.0, *disp, 0.0);
  WriteAccessVeln()->Update(1.0, *vel, 0.0);

  // to stop after this step:
  SetTimeEnd(time);
}

void ADAPTER::StructureInvana::PrePredict()
{
  if (singlesteponly_)
  {
    // set the correct internal variables for timestep Step()
    Teuchos::ParameterList p;
    p.set("timestep", Step());  // Step() gives the target step
    p.set("action", "calc_struct_recover_istep");
    Discretization()->Evaluate(
        p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  }

  return;
}

void ADAPTER::StructureInvana::PostPredict() { return; }

void ADAPTER::StructureInvana::PreUpdate()
{
  // do the update of the internal variables history before
  // quantities are prepared for the next target time step
  // via Update()
  Teuchos::ParameterList p;
  p.set("timestep", Step());
  p.set("action", "calc_struct_store_istep");
  Discretization()->Evaluate(
      p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}

void ADAPTER::StructureInvana::PostUpdate()
{
  // state and time was already updated so Dispn() and Veln()
  // (pointing to "(*state_)(0)") give the results of
  // this time step
  sdis_->insert(std::make_pair(StepOld(), *Dispn()));
  svel_->insert(std::make_pair(StepOld(), *Veln()));
  (*stime_)[StepOld()] = TimeOld();


  return;
}
