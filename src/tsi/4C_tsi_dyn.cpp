/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for thermo-structure-interaction problems.


\level 2
*/


/*----------------------------------------------------------------------*
 | definitions                                               dano 12/09 |
 *----------------------------------------------------------------------*/
#include "4C_tsi_dyn.hpp"

#include "4C_adapter_str_structure.hpp"
#include "4C_adapter_thermo.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_tsi.hpp"
#include "4C_tsi_algorithm.hpp"
#include "4C_tsi_monolithic.hpp"
#include "4C_tsi_partitioned.hpp"
#include "4C_tsi_utils.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | entry point for TSI in DRT                                dano 12/09 |
 *----------------------------------------------------------------------*/
void tsi_dyn_drt()
{
  // create a communicator
  const Epetra_Comm& comm = Global::Problem::instance()->get_dis("structure")->get_comm();

  // print TSI-Logo to screen
  if (comm.MyPID() == 0) TSI::printlogo();

  // setup of the discretizations, including clone strategy
  TSI::UTILS::SetupTSI(comm);

  // access the problem-specific parameter list
  const Teuchos::ParameterList& tsidyn = Global::Problem::instance()->tsi_dynamic_params();
  // access the problem-specific parameter list
  const Teuchos::ParameterList& sdynparams =
      Global::Problem::instance()->structural_dynamic_params();
  const Inpar::TSI::SolutionSchemeOverFields coupling =
      Core::UTILS::IntegralValue<Inpar::TSI::SolutionSchemeOverFields>(tsidyn, "COUPALGO");

  // create an empty TSI::Algorithm instance
  Teuchos::RCP<TSI::Algorithm> tsi;

  // choose algorithm depending on solution type
  switch (coupling)
  {
    case Inpar::TSI::Monolithic:
    {
      // create an TSI::Monolithic instance
      tsi = Teuchos::rcp(new TSI::Monolithic(comm, sdynparams));
      break;
    }
    case Inpar::TSI::OneWay:
    case Inpar::TSI::SequStagg:
    case Inpar::TSI::IterStagg:
    case Inpar::TSI::IterStaggAitken:
    case Inpar::TSI::IterStaggAitkenIrons:
    case Inpar::TSI::IterStaggFixedRel:
    {
      // Any partitioned algorithm. Stable of working horses.
      // create an TSI::Algorithm instance
      tsi = Teuchos::rcp(new TSI::Partitioned(comm));
      break;
    }
    default:
      FOUR_C_THROW("Unknown solutiontype for thermo-structure interaction: %d", coupling);
      break;
  }  // end switch

  const int restart = Global::Problem::instance()->restart();
  if (restart)
  {
    // read the restart information, set vectors and variables
    tsi->read_restart(restart);
  }

  // now do the coupling setup and create the combined dofmap
  tsi->setup_system();

  // solve the whole tsi problem
  tsi->time_loop();

  // summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test
  Global::Problem::instance()->add_field_test(tsi->structure_field()->create_field_test());
  Global::Problem::instance()->add_field_test(tsi->thermo_field()->create_field_test());
  Global::Problem::instance()->test_all(comm);

  return;
}  // tsi_dyn_drt()


/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
