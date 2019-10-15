/*----------------------------------------------------------------------*/
/*! \file

 \brief control routine of poroelasticity problems

\level 2

\maintainer Johannes Kremheller
            kremheller@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249
 *------------------------------------------------------------------------------------------------*/


/*----------------------------------------------------------------------*
 |  headers                                                             |
 *----------------------------------------------------------------------*/
#include "poro_dyn.H"

#include <Teuchos_TimeMonitor.hpp>

#include "poro_base.H"
#include "poro_scatra_base.H"
#include "poroelast_utils.H"
#include "poroelast_utils_setup.H"
#include "poro_utils_clonestrategy.H"

#include "../drt_inpar/inpar_poroelast.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"


/*------------------------------------------------------------------------------------------------*
 | main control routine for poroelasticity problems                                   vuong 01/12 |
 *------------------------------------------------------------------------------------------------*/
void poroelast_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // create a communicator
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  // print Logo to screen
  if (comm.MyPID() == 0) POROELAST::PrintLogo();

  // setup of the discretizations, including clone strategy
  POROELAST::UTILS::SetupPoro<POROELAST::UTILS::PoroelastCloneStrategy>();

  // access the problem-specific parameter list
  const Teuchos::ParameterList& poroelastdyn = problem->PoroelastDynamicParams();

  // choose algorithm depending on solution type
  Teuchos::RCP<POROELAST::PoroBase> poroalgo =
      POROELAST::UTILS::CreatePoroAlgorithm(poroelastdyn, comm);

  // read the restart information, set vectors and variables
  const int restart = problem->Restart();
  poroalgo->ReadRestart(restart);

  // now do the coupling setup and create the combined dofmap
  poroalgo->SetupSystem();

  // solve the whole problem
  poroalgo->TimeLoop();

  // summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test
  poroalgo->TestResults(comm);

  return;
}  // poroelast_drt()

/*------------------------------------------------------------------------------------------------*
 | main control routine for poro scatra problems                                   vuong 01/12 |
 *------------------------------------------------------------------------------------------------*/
void poro_scatra_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // 1.- Initialization
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  // 2.- Parameter reading
  const Teuchos::ParameterList& poroscatradynparams = problem->PoroScatraControlParams();

  POROELAST::UTILS::SetupPoroScatraDiscretizations<POROELAST::UTILS::PoroelastCloneStrategy,
      POROELAST::UTILS::PoroScatraCloneStrategy>();

  // 3.- Creation of Poroelastic + Scalar_Transport problem. (Discretization called inside)
  Teuchos::RCP<POROELAST::PoroScatraBase> poro_scatra =
      POROELAST::UTILS::CreatePoroScatraAlgorithm(poroscatradynparams, comm);

  // 3.1- Read restart if needed. (Discretization called inside)
  const int restart = problem->Restart();
  poro_scatra->ReadRestart(restart);

  // 4.- Run of the actual problem.

  // 4.1.- Some setup needed for the poroelastic subproblem.
  poro_scatra->SetupSystem();

  // 4.2.- Solve the whole problem
  poro_scatra->Timeloop();

  // 4.3.- Summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // 5. - perform the result test
  poro_scatra->TestResults(comm);
}
