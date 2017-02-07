/*!----------------------------------------------------------------------
\file art_net_dyn_drt.cpp
\brief Main control routine for all arterial network  solvers,

     including instationary solvers based on

     o

\maintainer Lena Yoshihara

\level 3

*----------------------------------------------------------------------*/
#ifdef D_ARTNET

#include <ctime>
#include <cstdlib>
#include <iostream>

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include "art_net_dyn_drt.H"
#include "artery_resulttest.H"

#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_inpar/drt_validparameters.H"



/*----------------------------------------------------------------------*
 * Main control routine for arterial network including various solvers:
 *
 *        o
 *
 *----------------------------------------------------------------------*/
void dyn_art_net_drt()
{
  dyn_art_net_drt(false);
}

Teuchos::RCP<ART::ArtNetExplicitTimeInt> dyn_art_net_drt(bool CoupledTo3D)
{
  if(DRT::Problem::Instance()->DoesExistDis("artery")==false)
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

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;

  actdis = DRT::Problem::Instance()->GetDis("artery");

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // If discretization is empty, then return empty time integration
  // -------------------------------------------------------------------
  if (actdis->NumGlobalElements()<1)
  {
    return Teuchos::null;
  }

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter>  output = actdis->Writer();
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  // const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  //  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& artdyn   = DRT::Problem::Instance()->ArterialDynamicParams();

  if (actdis->Comm().MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(IO::cout, artdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  // get the solver number
  const int linsolvernumber = artdyn.get<int>("LINEAR_SOLVER");
  // check if the solver has a valid solver number
  if (linsolvernumber == (-1))
    dserror("no linear solver defined. Please set LINEAR_SOLVER in ARTERIAL DYNAMIC to a valid number!");
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp( new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                                                       actdis->Comm(),
                                                       DRT::Problem::Instance()->ErrorFile()->Handle()),false );
  actdis->ComputeNullSpaceIfNecessary(solver->Params());

  // -------------------------------------------------------------------
  // set parameters in list required for all schemes
  // -------------------------------------------------------------------
  Teuchos::ParameterList arterytimeparams;

  // -------------------------------------- number of degrees of freedom
  // number of degrees of freedom
  const int ndim = DRT::Problem::Instance()->NDim();
  arterytimeparams.set<int>              ("number of degrees of freedom" ,2*ndim);

  // -------------------------------------------------- time integration
  // the default time step size
  arterytimeparams.set<double>           ("time step size"           ,artdyn.get<double>("TIMESTEP"));
  // maximum number of timesteps
  arterytimeparams.set<int>              ("max number timesteps"     ,artdyn.get<int>("NUMSTEP"));

  // ----------------------------------------------- restart and output
  // restart
  arterytimeparams.set                  ("write restart every"       ,artdyn.get<int>("RESTARTEVRY"));
  // solution output
  arterytimeparams.set                  ("write solution every"      ,artdyn.get<int>("RESULTSEVRY"));

  // solution output
  if (artdyn.get<std::string>("SOLVESCATRA")=="yes")
  {
    std::cout<<"Scatra will be solved"<<std::endl;
    arterytimeparams.set ("solve scatra", true);
  }
  else
  {
    std::cout<<"Scatra will not be solved"<<std::endl;
    arterytimeparams.set ("solve scatra", false);
  }

  // flag for writing the hemodynamic physiological results
  //arterytimeparams.set ("write stresses"  ,DRT::INPUT::IntegralValue<int>(ioflags,"HEMO_PHYS_RESULTS"));
  //---------------------- A method to initialize the flow inside the
  //                       arteries.
  //  int init = DRT::INPUT::IntegralValue<int> (artdyn,"INITIALFIELD");

  //------------------------------------------------------------------
  // create all vectors and variables associated with the time
  // integration (call the constructor);
  // the only parameter from the list required here is the number of
  // velocity degrees of freedom
  //------------------------------------------------------------------
  Teuchos::RCP<ART::ArtNetExplicitTimeInt> artnetexplicit
    =
    Teuchos::rcp(new ART::ArtNetExplicitTimeInt(actdis,*solver,arterytimeparams,*output),false);

  // Initialize state save vectors
  if (CoupledTo3D)
  {
    artnetexplicit->InitSaveState();
  }

  // initial field from restart or calculated by given function
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart && !CoupledTo3D)
  {
    // read the restart information, set vectors and variables
    artnetexplicit->ReadRestart(restart);
  }
  else
  {
    // artnetexplicit.SetInitialData(init,startfuncno);
  }

  arterytimeparams.set<FILE*>("err file",DRT::Problem::Instance()->ErrorFile()->Handle());

  if (!CoupledTo3D)
  {
    // call time-integration (or stationary) scheme
    Teuchos::RCP<Teuchos::ParameterList> param_temp;
    artnetexplicit->Integrate(CoupledTo3D,param_temp);

    /*
    Teuchos::RCP<DRT::ResultTest> resulttest
      = Teuchos::rcp(new ART::ArteryResultTest(*artnetexplicit));
    DRT::Problem::Instance()->AddFieldTest(resulttest);
    DRT::Problem::Instance()->TestAll(actdis->Comm());
    */
    artnetexplicit->TestResults();

    return artnetexplicit;
    //    return  Teuchos::null;
  }
  else
  {
    return artnetexplicit;
  }

} // end of dyn_art_net_drt()


#endif //#ifdef D_ARTNET
