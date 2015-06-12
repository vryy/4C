/*!------------------------------------------------------------------------------------------------*
 \file ssi_partitioned_1wc.cpp

 \brief one way coupled partitioned scalar structure interaction

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *------------------------------------------------------------------------------------------------*/

#include "ssi_partitioned_1wc.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"

#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_scatra/scatra_timint_implicit.H"

SSI::SSI_Part1WC::SSI_Part1WC(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams)
  : SSI_Part(comm, globaltimeparams,scatraparams,structparams),
    isscatrafromfile_(false)
{

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC::DoStructStep()
{

  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  structure_-> Solve();
  // calculate stresses, strains, energies
  structure_-> PrepareOutput();
  // update all single field solvers
  structure_-> Update();
  // write output to files
  structure_-> Output();
  // write output to screen
  structure_->PrintStep();
  // clean up
  structure_->Discretization()->ClearState(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC::DoScatraStep()
{

  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n TRANSPORT SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  if(isscatrafromfile_)
  {
    int diffsteps = structure_->Dt()/scatra_->ScaTraField()->Dt();
    if (scatra_->ScaTraField()->Step() % diffsteps ==0){
      scatra_->ScaTraField()->ReadRestart(scatra_->ScaTraField()->Step()); // read results from restart file
    }
  }
  else scatra_->ScaTraField()->Solve(); // really solve scatra problem


  // -------------------------------------------------------------------
  //                         update solution
  //        current solution becomes old solution of next timestep
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->Update();

  // -------------------------------------------------------------------
  // evaluate error for problems with analytical solution
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->EvaluateErrorComparedToAnalyticalSol();

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->Output();

  //cleanup
  scatra_->ScaTraField()->Discretization()->ClearState(true);
}

/*----------------------------------------------------------------------*/
//prepare time step
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC_SolidToScatra::PrepareTimeStep()
{
  IncrementTimeAndStep();
  //PrintHeader();

  structure_-> PrepareTimeStep();

  const int diffsteps = scatra_->ScaTraField()->Dt()/structure_->Dt();

  if (structure_->Step() % diffsteps == 0)
  {
    SetStructSolution(structure_->Dispn(),structure_->Veln());
    scatra_->ScaTraField()->PrepareTimeStep();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Part1WC_SolidToScatra::SSI_Part1WC_SolidToScatra(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams)
  : SSI_Part1WC(comm, globaltimeparams, scatraparams, structparams)
{
  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<DRT::DofSet> structdofset
    = structure_->Discretization()->GetDofSetProxy();

  //do some checks
  {
    INPAR::SCATRA::ConvForm convform
    = DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(scatraparams,"CONVFORM");
    if ( convform != INPAR::SCATRA::convform_conservative )
      dserror("If the scalar tranport problem is solved on the deforming domain, the conservative form must be"
          "used to include volume changes! Set 'CONVFORM' to 'conservative' in the SCALAR TRANSPORT DYNAMIC section!");
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC_SolidToScatra::Timeloop()
{
  //InitialCalculations();

  if (structure_->Dt() > scatra_->ScaTraField()->Dt())
    dserror("Timestepsize of scatra should be equal or bigger than solid timestep in solid to scatra interaction");

  const int diffsteps = scatra_->ScaTraField()->Dt()/structure_->Dt();

  while (NotFinished())
  {
    PrepareTimeStep();
    DoStructStep();  // It has its own time and timestep variables, and it increments them by itself.
    if (structure_->Step() % diffsteps == 0)
    {
      SetStructSolution(structure_->Dispnp(),structure_->Velnp());
      DoScatraStep();  // It has its own time and timestep variables, and it increments them by itself.
    }
  }

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Part1WC_ScatraToSolid::SSI_Part1WC_ScatraToSolid(const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& scatraparams,
    const Teuchos::ParameterList& structparams)
  : SSI_Part1WC(comm, globaltimeparams, scatraparams, structparams)
{
  // Flag for reading scatra result from restart file instead of computing it
  isscatrafromfile_ = DRT::INPUT::IntegralValue<bool>(DRT::Problem::Instance()->SSIControlParams(),"SCATRA_FROM_RESTART_FILE");

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC_ScatraToSolid::Timeloop()
{

  if (structure_->Dt() < scatra_->ScaTraField()->Dt())
  {
    dserror("Timestepsize of solid should be equal or bigger than scatra timestep in scatra to solid interaction");
  }

  const int diffsteps = structure_->Dt()/scatra_->ScaTraField()->Dt();
  while (NotFinished())
  {
    PrepareTimeStep();
    DoScatraStep();  // It has its own time and timestep variables, and it increments them by itself.
    if (scatra_->ScaTraField()->Step()  % diffsteps ==0)
    {
      SetScatraSolution(scatra_->ScaTraField()->Phinp());
      // PrepareTimeStep() is called after solving the scalar transport, because then the predictor will include the
      // new scalar solution
      structure_-> PrepareTimeStep();
      DoStructStep();  // It has its own time and timestep variables, and it increments them by itself.
    }
  }

}

/*----------------------------------------------------------------------*/
//prepare time step
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC_ScatraToSolid::PrepareTimeStep()
{
  IncrementTimeAndStep();
  //PrintHeader();

  scatra_->ScaTraField()->PrepareTimeStep();
  //PrepareTimeStep of structure field is called later
}
