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

  // solve the step problem
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

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
SSI::SSI_Part1WC_SolidToScatra::SSI_Part1WC_SolidToScatra(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams)
  : SSI_Part1WC(comm, timeparams)
{
  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<DRT::DofSet> structdofset
    = structure_->Discretization()->GetDofSetProxy();

  // check if scatra field has 2 discretizations, so that coupling is possible
  if (scatra_->ScaTraField().Discretization()->AddDofSet(structdofset)!=1)
    dserror("unexpected dof sets in scatra field");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC_SolidToScatra::Timeloop()
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
SSI::SSI_Part1WC_ScatraToSolid::SSI_Part1WC_ScatraToSolid(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams)
  : SSI_Part1WC(comm, timeparams)
{
  // build a proxy of the scatra discretization for the structure field
  Teuchos::RCP<DRT::DofSet> scatradofset
    = scatra_->ScaTraField().Discretization()->GetDofSetProxy();

  // check if structure field has 2 discretizations, so that coupling is possible
  if (structure_->Discretization()->AddDofSet(scatradofset)!=1)
    dserror("unexpected dof sets in structure field");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void SSI::SSI_Part1WC_ScatraToSolid::Timeloop()
{
  //InitialCalculations();

  while (NotFinished())
  {
    PrepareTimeStep();

    DoScatraStep(); // It has its own time and timestep variables, and it increments them by itself.
    SetScatraSolution();
    DoStructStep(); // It has its own time and timestep variables, and it increments them by itself.
  }
}
