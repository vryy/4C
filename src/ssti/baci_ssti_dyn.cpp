/*----------------------------------------------------------------------*/
/*! \file
 \brief control routine for scalar structure thermo interaction

 \level 2

 *------------------------------------------------------------------------------------------------*/

#include "baci_ssti_dyn.hpp"

#include "baci_global_data.hpp"
#include "baci_inpar_ssti.hpp"
#include "baci_io_control.hpp"
#include "baci_lib_discret.hpp"
#include "baci_ssti_algorithm.hpp"
#include "baci_ssti_monolithic.hpp"
#include "baci_ssti_utils.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ssti_drt()
{
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  auto ssti = SSTI::BuildSSTI(Teuchos::getIntegralValue<INPAR::SSTI::SolutionScheme>(
                                  problem->SSTIControlParams(), "COUPALGO"),
      comm, problem->SSTIControlParams());

  ssti->Init(comm, problem->SSTIControlParams(), problem->ScalarTransportDynamicParams(),
      problem->SSTIControlParams().sublist("THERMO"), problem->StructuralDynamicParams());

  ssti->Setup();

  const int restart = problem->Restart();
  if (restart) ssti->ReadRestart(restart);

  ssti->SetupSystem();

  ssti->Timeloop();

  Teuchos::TimeMonitor::summarize();

  ssti->TestResults(comm);
}

FOUR_C_NAMESPACE_CLOSE
