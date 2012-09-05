/*!---------------------------------------------------------------------*
\file combust_dyn.cpp

\brief control of dynamic combustion analysis

<pre>
Maintainer: Florian Henke
            henke@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15265
</pre>
*----------------------------------------------------------------------*/

#include "combust_dyn.H"
#include "combust_utils.H"
#include "combust_algorithm.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_comm/comm_utils.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_Time.h>


/*------------------------------------------------------------------------------------------------*
 | main  control routine for dynamic combustion analysis                              henke 06/08 |
 *------------------------------------------------------------------------------------------------*/
void combust_dyn()
{
  // create a communicator
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();

  //------------------------------------------------------------------------------------------------
  // print COMBUST-Logo on screen
  //------------------------------------------------------------------------------------------------
  if (comm.MyPID()==0) COMBUST::printCombustLogo();

  //------------------------------------------------------------------------------------------------
  // create G-function discretization by copying the fluid discretization (fill with scatra elements)
  //------------------------------------------------------------------------------------------------

  // access fluid discretization
  RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->GetDis("fluid");
  if (!fluiddis->Filled()) fluiddis->FillComplete(false,false,false);
  if (fluiddis->NumGlobalNodes()==0)
    dserror("No fluid discretization found!");

  // access G-function discretization (it should be empty)
  RCP<DRT::Discretization> gfuncdis = DRT::Problem::Instance()->GetDis("scatra");
  /* remark: here the ScaTra discretization is renamed according to its physical
   * meaning in the combustion context, namely the G-function.*/
  if (!gfuncdis->Filled()) gfuncdis->FillComplete(false,false,false);

  if (gfuncdis->NumGlobalNodes()==0)
  {
    Epetra_Time time(comm);

    // access the scalar transport parameter list
    const Teuchos::ParameterList& scatracontrol = DRT::Problem::Instance()->ScalarTransportDynamicParams();
    // fetch the desired material id for the transport elements
    const int matid = scatracontrol.get<int>("MATID");

    // create the scatra discretization
    {
    Teuchos::RCP<DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy> > clonewizard =
          Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy>() );

    clonewizard->CreateMatchingDiscretization(fluiddis,gfuncdis,matid);
    }
    if (comm.MyPID()==0)
      cout<<"Created G-function discretization from fluid discretization in...."
          <<time.ElapsedTime() << " secs\n\n";
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
  Teuchos::RCP<COMBUST::Algorithm> combust_ = Teuchos::rcp(new COMBUST::Algorithm(comm, combustdyn, DRT::Problem::Instance()->SolverParams(linsolvernumber)));

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
    combust_->Restart(restart, restartscatrainput, restartfromfluid);
  }
  //------------------------------------------------------------------------------------------------
  // call one of the available time integration schemes
  //------------------------------------------------------------------------------------------------
  INPAR::FLUID::TimeIntegrationScheme timeintscheme = DRT::INPUT::IntegralValue<INPAR::FLUID::TimeIntegrationScheme>(combustdyn,"TIMEINT");

  if (timeintscheme == INPAR::FLUID::timeint_one_step_theta or
      timeintscheme == INPAR::FLUID::timeint_afgenalpha)
  {
    // solve a dynamic combustion problem
    combust_->TimeLoop();
    /* remark: Hier kann auch z.B. Genalpha mit combust->TimeLoop() gerufen werden, weil der
     * combustion Algorithmus ja schon weiss welche Zeitintegration er hat. Es muss dann eine Klasse
     * "GenalphaTimeInt" existieren, die eine Funktion TimeLoop() hat. Dann muss allerdings auch
     * das ADAPTER::FluidCombust ein entsprechendes Object GenalphaTimeInt haben. Momentan hat ein
     * FluidCombust immer ein CombustFluidImplicitTimeInt!
     */
  }
  else if (timeintscheme == INPAR::FLUID::timeint_stationary)
  {
    // solve a static combustion problem
    combust_->SolveStationaryProblem();
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
#ifdef TRILINOS_DEV
  Teuchos::RCP<const Teuchos::Comm<int> > TeuchosComm = COMM_UTILS::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, false);
#else
  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);
#endif

  // perform the result test
  DRT::Problem::Instance()->AddFieldTest(combust_->FluidField().CreateFieldTest());
  DRT::Problem::Instance()->AddFieldTest(combust_->CreateScaTraFieldTest());
  DRT::Problem::Instance()->TestAll(comm);

  return;

} // combust_dyn()
