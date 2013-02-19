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

SSI::SSI_Part1WC::SSI_Part1WC(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams)
  : SSI_Part(comm, timeparams)
{
  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<DRT::DofSet> structdofset
    = structure_->Discretization()->GetDofSetProxy();
  // build a proxy of the temperature discretization for the structure field
  Teuchos::RCP<DRT::DofSet> scatradofset
    = scatra_->ScaTraField().Discretization()->GetDofSetProxy();

  // check if scatra field has 2 discretizations, so that coupling is possible
  if (scatra_->ScaTraField().Discretization()->AddDofSet(structdofset)!=1)
    dserror("unexpected dof sets in scatra field");
  if (structure_->Discretization()->AddDofSet(scatradofset)!=1)
    dserror("unexpected dof sets in structure field");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC::Timeloop()
{
  //InitialCalculations();

  while (NotFinished())
  {
    PrepareTimeStep();

    DoStructStep(); // It has its own time and timestep variables, and it increments them by itself.
    SetStructSolution();
    DoScatraStep(); // It has its own time and timestep variables, and it increments them by itself.
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC::DoStructStep()
{

  if (Comm().MyPID() == 0)
  {
    cout
        << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // solve the step problem. Methods obtained from poroelast->TimeLoop(sdynparams); --> sdynparams
  structure_-> PrepareTimeStep();
  // Newton-Raphson iteration
  structure_-> Solve();
  // calculate stresses, strains, energies
  structure_-> PrepareOutput();
  // update all single field solvers
  structure_-> Update();
  // write output to screen and files
  structure_-> Output();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC::DoScatraStep()
{

  if (Comm().MyPID() == 0)
  {
    cout
        << "\n***********************\n TRANSPORT SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_->ScaTraField().PrepareTimeStep();

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_->ScaTraField().Solve();

  // -------------------------------------------------------------------
  //                         update solution
  //        current solution becomes old solution of next timestep
  // -------------------------------------------------------------------
  scatra_->ScaTraField().Update();

  // -------------------------------------------------------------------
  // evaluate error for problems with analytical solution
  // -------------------------------------------------------------------
  scatra_->ScaTraField().EvaluateErrorComparedToAnalyticalSol();

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  scatra_->ScaTraField().Output();
}

/*----------------------------------------------------------------------*/
//prepare time step
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC::PrepareTimeStep()
{
  IncrementTimeAndStep();
  PrintHeader();

  //PrepareTimeStep of single fields is called in DoStructStep and DoScatraStep
}
