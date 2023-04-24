/*----------------------------------------------------------------------*/
/*! \file

 \brief control routine of poroelasticity coupled with scalar transport problems

\level 2

 *------------------------------------------------------------------------------------------------*/

#include "poroelast_scatra_dyn.H"

#include <Teuchos_TimeMonitor.hpp>


#include "poroelast_scatra_base.H"
#include "poroelast_scatra_utils.H"
#include "poroelast_scatra_utils_setup.H"
#include "poroelast_scatra_utils_clonestrategy.H"
#include "poroelast_utils_clonestrategy.H"


void poro_scatra_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  // 1.- Initialization
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  // 2.- Parameter reading
  const Teuchos::ParameterList& poroscatradynparams = problem->PoroScatraControlParams();

  POROELASTSCATRA::UTILS::SetupPoroScatraDiscretizations<
      POROELASTSCATRA::UTILS::PoroelastCloneStrategyforScatraElements,
      POROELASTSCATRA::UTILS::PoroScatraCloneStrategy>();

  // 3.- Creation of Poroelastic + Scalar_Transport problem. (Discretization called inside)
  Teuchos::RCP<POROELASTSCATRA::PoroScatraBase> poro_scatra =
      POROELASTSCATRA::UTILS::CreatePoroScatraAlgorithm(poroscatradynparams, comm);

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
