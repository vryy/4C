/*----------------------------------------------------------------------*/
/*! \file

 \brief control routine of poroelasticity problems

\level 2

 *------------------------------------------------------------------------------------------------*/

#include "baci_poroelast_dyn.H"

#include "baci_poroelast_base.H"
#include "baci_poroelast_utils.H"
#include "baci_poroelast_utils_clonestrategy.H"
#include "baci_poroelast_utils_setup.H"

#include <Teuchos_TimeMonitor.hpp>


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
}
