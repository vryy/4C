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


#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif
#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>

#include "topopt_dyn.H"
#include "topopt_algorithm.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "topopt_utils.H"



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
  const Epetra_Comm& comm = DRT::Problem::Instance()->Dis(genprob.numff,0)->Comm();
#else
  Epetra_SerialComm comm;
#endif

  //------------------------------------------------------------------------------------------------
  // print COMBUST-Logo on screen
  //------------------------------------------------------------------------------------------------
  if (comm.MyPID()==0) TOPOPT::printTopOptLogo();

  //------------------------------------------------------------------------------------------------
  // create optimization discretization by copying the fluid discretization (fill with opti elements)
  //------------------------------------------------------------------------------------------------
  // get discretization ids
  int disnumff = genprob.numff; // discretization number fluid; typically 0
  int disnumof = genprob.numof; // discretization number optimization field; typically 1

  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  // access fluid discretization
  RCP<DRT::Discretization> fluiddis = problem->Dis(disnumff,0);
  if (!fluiddis->Filled()) fluiddis->FillComplete(false,false,false);
  if (fluiddis->NumGlobalNodes()==0)
    dserror("No fluid discretization found!");

  // access optimization discretization (it should be empty if it will be cloned)
  RCP<DRT::Discretization> optidis = problem->Dis(disnumof,0);
  if (!optidis->Filled()) optidis->FillComplete(false,false,false);

  if (optidis->NumGlobalNodes()==0)
  {
    Epetra_Time time(comm);

    //get the material map for cloning
    std::map<std::pair<string,string>,std::map<int,int> > clonefieldmatmap = problem->ClonedMaterialMap();

    // create the fluid discretization
    {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<TOPOPT::TopoptFluidCloneStrategy> > clonewizard
      = Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<TOPOPT::TopoptFluidCloneStrategy>() );

      std::pair<string,string> key("fluid","scatra");
      std::map<int,int> matmap = clonefieldmatmap[key];
      clonewizard->CreateMatchingDiscretization(fluiddis,optidis,matmap);
    }

    if (comm.MyPID()==0)
      cout<<"Created optimization discretization from fluid discretization in...."
          <<time.ElapsedTime() << " secs\n\n";
  }
  else
    dserror("Optimization discretization is not empty as it should be!");
  // TODO this shall be ok later if optimization has different discretization than fluid (winklmaier)
  // therefore later a fillcomplete has to be called here also!

  //------------------------------------------------------------------------------------------------
  // create a topology optimization algorithm
  //------------------------------------------------------------------------------------------------
  // get the topology optimization parameter list
  Teuchos::ParameterList topoptdyn = DRT::Problem::Instance()->OptimizationControlParams();
  // create a COMBUST::Algorithm instance
  Teuchos::RCP<TOPOPT::Algorithm> topopt_ = Teuchos::rcp(new TOPOPT::Algorithm(comm,topoptdyn));

  //------------------------------------------------------------------------------------------------
  // restart
  //------------------------------------------------------------------------------------------------
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    // check where we restart
    const int restartaction = DRT::INPUT::IntegralValue<int>(topoptdyn,"RESTART_ACTION");
cout << "test restart action: 0=fluid,1=adjoint,2=grad,3=opti-step: " << restartaction << endl;
    // read the restart information, set vectors and variables
    topopt_->Restart(restart,restartaction);
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
  problem->AddFieldTest(topopt_->FluidField().CreateFieldTest());
  problem->AddFieldTest(topopt_->AdjointFluidField()->CreateFieldTest());

  problem->TestAll(comm);

  return;

} // fluid_topopt_dyn()

