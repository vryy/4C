/*!----------------------------------------------------------------------
\brief Main control routine for all arterial network  solvers,

     including instationary solvers based on

     o

\maintainer Johannes Kremheller

\level 3

*----------------------------------------------------------------------*/

#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "art_net_dyn_drt.H"
#include "artery_resulttest.H"
#include "../drt_adapter/ad_art_net.H"
#include "art_net_utils.H"

#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/drt_validparameters.H"

#include "../drt_lib/drt_utils_createdis.H"

/*----------------------------------------------------------------------*
 * Main control routine for arterial network including various solvers:
 *
 *        o
 *
 *----------------------------------------------------------------------*/
void dyn_art_net_drt() { dyn_art_net_drt(false); }

Teuchos::RCP<ADAPTER::ArtNet> dyn_art_net_drt(bool CoupledTo3D)
{
  if (DRT::Problem::Instance()->DoesExistDis("artery") == false)
  {
#if 0
    if (actdis->Comm().MyPID()==0)
    {
      std::cout<<"+--------------------- WARNING ---------------------+"<<std::endl;
      std::cout<<"|                                                   |"<<std::endl;
      std::cout<<"| One-dimesional arterial network is compiled, but  |"<<std::endl;
      std::cout<<"| no artery elements are defined!                   |"<<std::endl;
      std::cout<<"|                                                   |"<<std::endl;
      std::cout<<"+---------------------------------------------------+"<<std::endl;
    }
#endif
    return Teuchos::null;
  }

  // define the discretization names
  const std::string artery_disname = "artery";
  const std::string scatra_disname = "artery_scatra";

  // access the problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;

  actdis = problem->GetDis(artery_disname);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // If discretization is empty, then return empty time integration
  // -------------------------------------------------------------------
  if (actdis->NumGlobalElements() < 1)
  {
    return Teuchos::null;
  }

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();
  output->WriteMesh(0, 0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  const Teuchos::ParameterList& artdyn = problem->ArterialDynamicParams();

  if (actdis->Comm().MyPID() == 0) DRT::INPUT::PrintDefaultParameters(IO::cout, artdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // get the solver number
  const int linsolvernumber = artdyn.get<int>("LINEAR_SOLVER");
  // check if the solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined. Please set LINEAR_SOLVER in ARTERIAL DYNAMIC to a valid "
        "number!");

  // solution output
  if (artdyn.get<std::string>("SOLVESCATRA") == "yes")
  {
    if (actdis->Comm().MyPID() == 0)
    {
      std::cout << "<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cout << "<  ARTERY:  ScaTra coupling present  >" << std::endl;
      std::cout << "<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>" << std::endl;
    }
    Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis(scatra_disname);
    // fill scatra discretization by cloning artery discretization
    DRT::UTILS::CloneDiscretization<ART::ArteryScatraCloneStrategy>(actdis, scatradis);
    scatradis->FillComplete();

    // the problem is one way coupled, scatra needs only artery

    // build a proxy of the structure discretization for the scatra field
    Teuchos::RCP<DRT::DofSetInterface> arterydofset = actdis->GetDofSetProxy();

    // check if ScatraField has 2 discretizations, so that coupling is possible
    if (scatradis->AddDofSet(arterydofset) != 1) dserror("unexpected dof sets in scatra field");

    scatradis->FillComplete(true, false, false);
  }
  else
  {
    if (actdis->Comm().MyPID() == 0)
    {
      std::cout << "<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>" << std::endl;
      std::cout << "<  ARTERY: no ScaTra coupling present  >" << std::endl;
      std::cout << "<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>" << std::endl;
    }
  }

  // flag for writing the hemodynamic physiological results
  // arterytimeparams.set ("write stresses"
  // ,DRT::INPUT::IntegralValue<int>(ioflags,"HEMO_PHYS_RESULTS"));
  //---------------------- A method to initialize the flow inside the
  //                       arteries.
  //  int init = DRT::INPUT::IntegralValue<int> (artdyn,"INITIALFIELD");

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  INPAR::ARTDYN::TimeIntegrationScheme timintscheme =
      DRT::INPUT::IntegralValue<INPAR::ARTDYN::TimeIntegrationScheme>(artdyn, "DYNAMICTYP");

  // build art net time integrator
  Teuchos::RCP<ADAPTER::ArtNet> artnettimint = ART::UTILS::CreateAlgorithm(timintscheme, actdis,
      linsolvernumber, artdyn, artdyn, DRT::Problem::Instance()->ErrorFile()->Handle(), output);

  // initialize
  artnettimint->Init(artdyn, artdyn, scatra_disname);

  // Initialize state save vectors
  if (CoupledTo3D)
  {
    artnettimint->InitSaveState();
  }

  // initial field from restart or calculated by given function
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart && !CoupledTo3D)
  {
    // read the restart information, set vectors and variables
    artnettimint->ReadRestart(restart);
  }
  else
  {
    // artnetexplicit.SetInitialData(init,startfuncno);
  }

  // assign materials
  // note: to be done after potential restart, as in ReadRestart()
  //       the secondary material is destroyed
  if (artdyn.get<std::string>("SOLVESCATRA") == "yes")
    ART::UTILS::AssignMaterialPointers(artery_disname, scatra_disname);

  if (!CoupledTo3D)
  {
    // call time-integration (or stationary) scheme
    Teuchos::RCP<Teuchos::ParameterList> param_temp;
    artnettimint->Integrate(CoupledTo3D, param_temp);

    // result test
    artnettimint->TestResults();

    return artnettimint;
    //    return  Teuchos::null;
  }
  else
  {
    return artnettimint;
  }

}  // end of dyn_art_net_drt()
