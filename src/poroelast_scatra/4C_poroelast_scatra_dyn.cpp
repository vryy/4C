/*----------------------------------------------------------------------*/
/*! \file

 \brief control routine of poroelasticity coupled with scalar transport problems

\level 2

 *------------------------------------------------------------------------------------------------*/

#include "4C_poroelast_scatra_dyn.hpp"

#include "4C_poroelast_scatra_base.hpp"
#include "4C_poroelast_scatra_utils.hpp"
#include "4C_poroelast_scatra_utils_clonestrategy.hpp"
#include "4C_poroelast_scatra_utils_setup.hpp"
#include "4C_poroelast_utils_clonestrategy.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


void poro_scatra_drt()
{
  Global::Problem* problem = Global::Problem::Instance();

  // 1.- Initialization
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  // 2.- Parameter reading
  const Teuchos::ParameterList& poroscatradynparams = problem->poro_scatra_control_params();

  PoroElastScaTra::UTILS::SetupPoroScatraDiscretizations<
      PoroElastScaTra::UTILS::PoroelastCloneStrategyforScatraElements,
      PoroElastScaTra::UTILS::PoroScatraCloneStrategy>();

  // 3.- Creation of Poroelastic + Scalar_Transport problem. (discretization called inside)
  Teuchos::RCP<PoroElastScaTra::PoroScatraBase> poro_scatra =
      PoroElastScaTra::UTILS::CreatePoroScatraAlgorithm(poroscatradynparams, comm);

  // 3.1- Read restart if needed. (discretization called inside)
  const int restart = problem->Restart();
  poro_scatra->read_restart(restart);

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

FOUR_C_NAMESPACE_CLOSE
