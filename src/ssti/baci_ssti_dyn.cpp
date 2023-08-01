/*----------------------------------------------------------------------*/
/*! \file
 \brief control routine for scalar structure thermo interaction

 \level 2

 *------------------------------------------------------------------------------------------------*/

#include "baci_ssti_dyn.H"

#include "baci_ssti_algorithm.H"
#include "baci_ssti_monolithic.H"
#include "baci_ssti_utils.H"

#include "baci_inpar_ssti.H"

#include "baci_io_control.H"

#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ssti_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

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
