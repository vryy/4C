/*----------------------------------------------------------------------*/
/*! \file
 \brief control routine for scalar structure thermo interaction

 \level 2

 *------------------------------------------------------------------------------------------------*/

#include "4C_ssti_dyn.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_ssti.hpp"
#include "4C_io_control.hpp"
#include "4C_ssti_algorithm.hpp"
#include "4C_ssti_monolithic.hpp"
#include "4C_ssti_utils.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ssti_drt()
{
  Global::Problem* problem = Global::Problem::Instance();

  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

  auto ssti = SSTI::BuildSSTI(Teuchos::getIntegralValue<Inpar::SSTI::SolutionScheme>(
                                  problem->SSTIControlParams(), "COUPALGO"),
      comm, problem->SSTIControlParams());

  ssti->init(comm, problem->SSTIControlParams(), problem->scalar_transport_dynamic_params(),
      problem->SSTIControlParams().sublist("THERMO"), problem->structural_dynamic_params());

  ssti->setup();

  const int restart = problem->Restart();
  if (restart) ssti->read_restart(restart);

  ssti->SetupSystem();

  ssti->Timeloop();

  Teuchos::TimeMonitor::summarize();

  ssti->TestResults(comm);
}

FOUR_C_NAMESPACE_CLOSE
