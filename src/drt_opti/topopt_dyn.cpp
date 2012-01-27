/*!------------------------------------------------------------------------------------------------*
\file topopt_dyn.cpp

\brief control routine of topology optimization for fluid domains

<pre>
Maintainer: Martin Winklmaier
            winklmaier@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifdef CCADISCRET

#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include "Epetra_Time.h"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_scatra/scatra_utils.H"
#include "topopt_algorithm.H"
#include "topopt_dyn.H"
#include "topopt_utils.H"
#include <Teuchos_TimeMonitor.hpp>



/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;



/*------------------------------------------------------------------------------------------------*
 | main control routine for fluid topology optimization                          winklmaier 12/11 |
 *------------------------------------------------------------------------------------------------*/
void fluid_topopt_dyn()
{
  // create a communicator
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  //------------------------------------------------------------------------------------------------
  // print COMBUST-Logo on screen
  //------------------------------------------------------------------------------------------------
  if (comm.MyPID()==0) TOPOPT::printTopOptLogo();

  //------------------------------------------------------------------------------------------------
  // create G-function discretization by copying the fluid discretization (fill with scatra elements)
  //------------------------------------------------------------------------------------------------
  // get discretization ids
  int disnumff = genprob.numff; // discretization number fluid; typically 0

  // access fluid discretization
  RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(disnumff,0);
  if (!fluiddis->Filled()) fluiddis->FillComplete(false,false,false);
  if (fluiddis->NumGlobalNodes()==0)
    dserror("No fluid discretization found!");

  // currently the optimization field has no own discretization
  // it uses the fluid node maps for creating the required vector(s)

//  //------------------------------------------------------------------------------------------------
//  // create a topology optimization algorithm
//  //------------------------------------------------------------------------------------------------
  // get the topology optimization parameter list
  Teuchos::ParameterList topoptdyn = DRT::Problem::Instance()->OptimizationControlParams();
  // create a COMBUST::Algorithm instance
  Teuchos::RCP<TOPOPT::Algorithm> topopt_ = Teuchos::rcp(new TOPOPT::Algorithm(comm,topoptdyn));

  //------------------------------------------------------------------------------------------------
  // restart
  //------------------------------------------------------------------------------------------------
  if (genprob.restart)
  {
    // check where we restart
    const int restartaction = DRT::INPUT::IntegralValue<int>(topoptdyn,"RESTART_ACTION");
cout << "test restart action: 0=fluid,1=adjoint,2=grad,3=opti-step: " << restartaction << endl;
    // read the restart information, set vectors and variables
    topopt_->Restart(genprob.restart,restartaction);
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
  DRT::Problem::Instance()->AddFieldTest(topopt_->FluidField().CreateFieldTest());
  DRT::Problem::Instance()->TestAll(comm);

  return;

} // fluid_topopt_dyn()

#endif  // #ifdef CCADISCRET
