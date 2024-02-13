/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for thermo-structure-interaction problems.


\level 2
*/


/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/
#include "baci_tsi_dyn.hpp"

#include "baci_adapter_str_structure.hpp"
#include "baci_adapter_thermo.hpp"
#include "baci_global_data.hpp"
#include "baci_inpar_tsi.hpp"
#include "baci_lib_discret.hpp"
#include "baci_tsi_algorithm.hpp"
#include "baci_tsi_monolithic.hpp"
#include "baci_tsi_partitioned.hpp"
#include "baci_tsi_utils.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_TimeMonitor.hpp>

BACI_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | entry point for TSI in DRT                                dano 12/09 |
 *----------------------------------------------------------------------*/
void tsi_dyn_drt()
{
  // create a communicator
  const Epetra_Comm& comm = GLOBAL::Problem::Instance()->GetDis("structure")->Comm();

  // print TSI-Logo to screen
  if (comm.MyPID() == 0) TSI::printlogo();

  // setup of the discretizations, including clone strategy
  TSI::UTILS::SetupTSI(comm);

  // access the problem-specific parameter list
  const Teuchos::ParameterList& tsidyn = GLOBAL::Problem::Instance()->TSIDynamicParams();
  // access the problem-specific parameter list
  const Teuchos::ParameterList& sdynparams = GLOBAL::Problem::Instance()->StructuralDynamicParams();
  const INPAR::TSI::SolutionSchemeOverFields coupling =
      INPUT::IntegralValue<INPAR::TSI::SolutionSchemeOverFields>(tsidyn, "COUPALGO");

  // create an empty TSI::Algorithm instance
  Teuchos::RCP<TSI::Algorithm> tsi;

  // choose algorithm depending on solution type
  switch (coupling)
  {
    case INPAR::TSI::Monolithic:
    {
      // create an TSI::Monolithic instance
      tsi = Teuchos::rcp(new TSI::Monolithic(comm, sdynparams));
      break;
    }
    case INPAR::TSI::OneWay:
    case INPAR::TSI::SequStagg:
    case INPAR::TSI::IterStagg:
    case INPAR::TSI::IterStaggAitken:
    case INPAR::TSI::IterStaggAitkenIrons:
    case INPAR::TSI::IterStaggFixedRel:
    {
      // Any partitioned algorithm. Stable of working horses.
      // create an TSI::Algorithm instance
      tsi = Teuchos::rcp(new TSI::Partitioned(comm));
      break;
    }
    default:
      dserror("Unknown solutiontype for thermo-structure interaction: %d", coupling);
      break;
  }  // end switch

  const int restart = GLOBAL::Problem::Instance()->Restart();
  if (restart)
  {
    // read the restart information, set vectors and variables
    tsi->ReadRestart(restart);
  }

  // now do the coupling setup and create the combined dofmap
  tsi->SetupSystem();

  // solve the whole tsi problem
  tsi->TimeLoop();

  // summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test
  GLOBAL::Problem::Instance()->AddFieldTest(tsi->StructureField()->CreateFieldTest());
  GLOBAL::Problem::Instance()->AddFieldTest(tsi->ThermoField()->CreateFieldTest());
  GLOBAL::Problem::Instance()->TestAll(comm);

  return;
}  // tsi_dyn_drt()


/*----------------------------------------------------------------------*/

BACI_NAMESPACE_CLOSE
