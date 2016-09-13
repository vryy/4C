/*!---------------------------------------------------------------------*
\file combust_dyn.cpp

\brief control of dynamic combustion analysis

\level 2

<pre>
\maintainer Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>
*----------------------------------------------------------------------*/
#include "../drt_comm/comm_utils.H"

#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_scatra/scatra_utils_clonestrategy.H"

#include "../drt_scatra_ele/scatra_ele.H"

#include <Epetra_Time.h>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>

#include "combust_algorithm.H"
#include "combust_dyn.H"


/*------------------------------------------------------------------------------------------------*
 | main control routine for dynamic combustion analysis                               henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void combust_dyn()
{
  // create a communicator
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();

  //------------------------------------------------------------------------------------------------
  // print COMBUST-Logo on screen
  //------------------------------------------------------------------------------------------------
  if (comm.MyPID()==0)
  {
//      IO::cout << "     ___            ___    \n"
//               << "    /   \\          /   \\ \n"
//               << "    \\_   \\        /  __/ \n"
//               << "     _\\   \\      /  /__  " << " Das ist               \n"
//               << "     \\___  \\____/   __/  " << " das Verbrennungsmodul \n"
//               << "         \\_       _/      " << " in BACI               \n"
//               << "           | @ @  \\_      " << "                       \n"
//               << "           |               " << " Der Elch wird bald    \n"
//               << "         _/     /\\        " << " ein feuerspeiender    \n"
//               << "        /o)  (o/\\ \\_     " << " Drache sein!          \n"
//               << "        \\_____/ /         \n"
//               << "          \\____/          \n"
//               << "                           " << IO::endl;
      IO::cout << "\n ENTER COMBUSTION AND TWO-PHASE FLOW MODULE ...\n" << IO::endl;
  }

  //------------------------------------------------------------------------------------------------
  // create G-function discretization by copying the fluid discretization (fill with scatra elements)
  //------------------------------------------------------------------------------------------------

  // access fluid discretization
  Teuchos::RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis("fluid");
  if (!fluiddis->Filled()) fluiddis->FillComplete(false,false,false);
  if (fluiddis->NumGlobalNodes()==0)
    dserror("No fluid discretization found!");

  // access G-function discretization (it should be empty)
  Teuchos::RCP<DRT::Discretization> gfuncdis = DRT::Problem::Instance()->GetDis("scatra");
  /* remark: here the ScaTra discretization is renamed according to its physical
   * meaning in the combustion context, namely the G-function.*/
  if (!gfuncdis->Filled()) gfuncdis->FillComplete(false,false,false);

  if (gfuncdis->NumGlobalNodes()==0)
  {
    // fill scatra discretization by cloning fluid discretization
    DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(fluiddis,gfuncdis);

    // set implementation type of cloned scatra elements to levelset
    for(int i=0; i<gfuncdis->NumMyColElements(); ++i)
    {
      DRT::ELEMENTS::Transport* element = dynamic_cast<DRT::ELEMENTS::Transport*>(gfuncdis->lColElement(i));
      if(element == NULL)
        dserror("Invalid element type!");
      else
        element->SetImplType(INPAR::SCATRA::impltype_levelset);
    }
  }

  else
    dserror("G-function discretization is not empty. Fluid and G-function already present!");

  //------------------------------------------------------------------------------------------------
  // create a combustion algorithm
  //------------------------------------------------------------------------------------------------
  // get the combustion parameter list
  Teuchos::ParameterList combustdyn = DRT::Problem::Instance()->CombustionDynamicParams();
  const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
  combustdyn.set<bool>("GMSH_OUTPUT",DRT::INPUT::IntegralValue<bool>(fdyn,"GMSH_OUTPUT"));

  // get linear solver id from SCALAR TRANSPORT DYNAMIC
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
  const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror("no linear solver defined for ELCH problem. Please set LINEAR_SOLVER in SCALAR TRANSPORT DYNAMIC to a valid number!");

  // create a COMBUST::Algorithm instance
  Teuchos::RCP<COMBUST::Algorithm> combust = Teuchos::rcp(new COMBUST::Algorithm(comm, combustdyn, DRT::Problem::Instance()->SolverParams(linsolvernumber)));
  combust->Init(
      combustdyn,
      DRT::Problem::Instance()->ScalarTransportDynamicParams(),
      DRT::Problem::Instance()->SolverParams(linsolvernumber) );
  combust->Setup();

  //------------------------------------------------------------------------------------------------
  // restart
  //------------------------------------------------------------------------------------------------
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    // check if we restart from standard fluid problem
    const bool restartfromfluid = (bool)DRT::INPUT::IntegralValue<int>(combustdyn,"RESTART_FROM_FLUID");
    // turn on/off read scatra restart from input file
    const bool restartscatrainput = (bool)DRT::INPUT::IntegralValue<int>(combustdyn,"RESTART_SCATRA_INPUT");

    // read the restart information, set vectors and variables
    combust->Restart(restart, restartscatrainput, restartfromfluid);
  }
  //------------------------------------------------------------------------------------------------
  // call one of the available time integration schemes
  //------------------------------------------------------------------------------------------------
  INPAR::FLUID::TimeIntegrationScheme timeintscheme = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(fdyn,"TIMEINTEGR");

  if (timeintscheme == INPAR::FLUID::timeint_one_step_theta or
      timeintscheme == INPAR::FLUID::timeint_afgenalpha or
      timeintscheme == INPAR::FLUID::timeint_bdf2)
  {
    // solve a dynamic combustion problem
    combust->TimeLoop();
  }
  else if (timeintscheme == INPAR::FLUID::timeint_stationary)
  {
    // solve a static combustion problem
    combust->SolveStationaryProblem();
  }
  else
  {
    // every time integration scheme must be either static or dynamic
    dserror("the combustion module can not handle this type of time integration scheme");
  }

  //------------------------------------------------------------------------------------------------
  // validate the results
  //------------------------------------------------------------------------------------------------
  // summarize the performance measurements
  Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = COMM_UTILS::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), IO::cout.cout_replacement(), false, true, false);
  IO::cout << IO::flush;

  // perform the result test
  combust->TestResults();

  return;

} // combust_dyn()
