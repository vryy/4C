/*----------------------------------------------------------------------*/
/*! \file

 \brief control routine of poroelasticity problems

\level 2

 *------------------------------------------------------------------------------------------------*/

#include "4C_poroelast_dyn.hpp"

#include "4C_poroelast_base.hpp"
#include "4C_poroelast_utils.hpp"
#include "4C_poroelast_utils_clonestrategy.hpp"
#include "4C_poroelast_utils_setup.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


void poroelast_drt()
{
  Global::Problem* problem = Global::Problem::instance();

  // create a communicator
  const Epetra_Comm& comm = problem->get_dis("structure")->get_comm();

  // print Logo to screen
  if (comm.MyPID() == 0) PoroElast::PrintLogo();

  // setup of the discretizations, including clone strategy
  PoroElast::UTILS::SetupPoro<PoroElast::UTILS::PoroelastCloneStrategy>();

  // access the problem-specific parameter list
  const Teuchos::ParameterList& poroelastdyn = problem->poroelast_dynamic_params();

  // choose algorithm depending on solution type
  Teuchos::RCP<PoroElast::PoroBase> poroalgo =
      PoroElast::UTILS::CreatePoroAlgorithm(poroelastdyn, comm);

  // read the restart information, set vectors and variables
  const int restart = problem->restart();
  poroalgo->read_restart(restart);

  // now do the coupling setup and create the combined dofmap
  poroalgo->setup_system();

  // solve the whole problem
  poroalgo->time_loop();

  // summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test
  poroalgo->test_results(comm);
}

FOUR_C_NAMESPACE_CLOSE
