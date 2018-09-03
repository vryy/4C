/*!------------------------------------------------------------------------------------------------*
\file topopt_dyn.cpp

\brief control routine of topology optimization for fluid domains

\level 2

<pre>
\maintainer Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/


#include <Epetra_MpiComm.h>
#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>

#include "topopt_dyn.H"
#include "topopt_algorithm.H"
#include "topopt_optimizer.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "topopt_utils.H"



/*------------------------------------------------------------------------------------------------*
 | main control routine for fluid topology optimization                          winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void fluid_topopt_dyn()
{
  // create a communicator
#ifdef PARALLEL
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();
#else
  Epetra_SerialComm comm;
#endif

  //------------------------------------------------------------------------------------------------
  // print Logo on screen
  //------------------------------------------------------------------------------------------------
  if (comm.MyPID() == 0) TOPOPT::printTopOptLogo();

  //------------------------------------------------------------------------------------------------
  // create optimization discretization by copying the fluid discretization (fill with opti
  // elements)
  //------------------------------------------------------------------------------------------------
  // get discretization ids

  DRT::Problem* problem = DRT::Problem::Instance();

  // access fluid discretization
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  if (!fluiddis->Filled()) fluiddis->FillComplete(false, false, false);
  if (fluiddis->NumGlobalNodes() == 0) dserror("No fluid discretization found!");

  // access optimization discretization (it should be empty if it will be cloned)
  Teuchos::RCP<DRT::Discretization> optidis = problem->GetDis("opti");
  if (!optidis->Filled()) optidis->FillComplete(false, false, false);

  if (optidis->NumGlobalNodes() == 0)
  {
    DRT::UTILS::CloneDiscretization<TOPOPT::TopoptFluidCloneStrategy>(fluiddis, optidis);
    optidis->FillComplete();
  }
  else
    dserror("Optimization discretization is not empty as it should be!");
  // TODO this shall be ok later if optimization has different discretization than fluid
  // (winklmaier) therefore later a fillcomplete has to be called here also!

  //------------------------------------------------------------------------------------------------
  // create a topology optimization algorithm
  //------------------------------------------------------------------------------------------------
  // get the topology optimization parameter list
  Teuchos::ParameterList topoptdyn = problem->OptimizationControlParams();
  // create a Algorithm instance
  Teuchos::RCP<TOPOPT::Algorithm> topopt_ = Teuchos::rcp(new TOPOPT::Algorithm(comm, topoptdyn));

  //------------------------------------------------------------------------------------------------
  // restart
  //------------------------------------------------------------------------------------------------
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    // check where we restart
    const INPAR::TOPOPT::Restart restartaction =
        DRT::INPUT::IntegralValue<INPAR::TOPOPT::Restart>(topoptdyn, "RESTART_ACTION");

    // read the restart information, set vectors and variables
    topopt_->Restart(restart, restartaction);
  }

  //------------------------------------------------------------------------------------------------
  // call the optimization routine
  //------------------------------------------------------------------------------------------------
  topopt_->OptimizationLoop();

  //------------------------------------------------------------------------------------------------
  // validate the results
  //------------------------------------------------------------------------------------------------
  // summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test
  problem->AddFieldTest(topopt_->FluidField()->CreateFieldTest());
  problem->AddFieldTest(topopt_->AdjointFluidField()->CreateFieldTest());
  problem->AddFieldTest(topopt_->Optimizer()->CreateFieldTest());
  problem->TestAll(comm);

  return;

}  // fluid_topopt_dyn()
