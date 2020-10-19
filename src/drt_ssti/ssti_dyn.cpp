/*----------------------------------------------------------------------*/
/*! \file
 \brief control routine for scalar structure thermo interaction

 \level 2

 *------------------------------------------------------------------------------------------------*/

#include "ssti_dyn.H"

#include "ssti_algorithm.H"
#include "ssti_monolithic.H"
#include "ssti_utils.H"

#include "../drt_inpar/inpar_ssti.H"

#include "../drt_io/io_control.H"

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ssti_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  SSTI::PrintSSTILogo(problem->GetDis("structure")->Comm().MyPID());

  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  auto ssti = SSTI::BuildSSTI(Teuchos::getIntegralValue<INPAR::SSTI::SolutionScheme>(
                                  problem->SSTIControlParams(), "COUPALGO"),
      comm, problem->SSTIControlParams());

  ssti->Init(comm, problem->SSTIControlParams(), problem->ScalarTransportDynamicParams(),
      problem->SSTIControlParams().sublist("THERMO"), problem->StructuralDynamicParams());

  ssti->Setup();

  const int restart = problem->Restart();
  const double restarttime = problem->RestartTime();

  if (restarttime > 0.0)
    dserror("Restart from time is not supported.");
  else if (restart)
    ssti->ReadRestart(restart);

  ssti->SetupSystem();

  ssti->Timeloop();

  Teuchos::TimeMonitor::summarize();

  ssti->TestResults(comm);
}
