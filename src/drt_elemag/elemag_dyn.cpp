/*----------------------------------------------------------------------*/
/*!
\file elemag_dyn.cpp

\brief Main control routine for electromagnetic simulations

<pre>
\level 3

\maintainer Volker Gravemeier
            gravemeier@lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/

#include <iostream>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "elemag_dyn.H"
#include "elemag_timeint.H"
//#include "elemag_ele.H"
//#include "elemag_impl.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/inpar_elemag.H"
#include "../drt_lib/drt_discret_hdg.H"
#include "../drt_lib/drt_dofset_independent.H"
#include "../linalg/linalg_solver.H"
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_lib/drt_dofset_predefineddofnumber.H"

void electromagnetics_drt()
{
  // declare abbreviation
  DRT::Problem* problem = DRT::Problem::Instance();

  // declare problem-specific parameter list for electromagnetics
  const Teuchos::ParameterList& elemagparams = problem->ElectromagneticParams();

  // declare discretization and check their existence
  Teuchos::RCP<DRT::DiscretizationHDG> elemagdishdg =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationHDG>(problem->GetDis("elemag"));
  if (elemagdishdg == Teuchos::null)
    dserror("Failed to cast DRT::Discretization to DRT::DiscretizationHDG.");

  // declare communicator and print module information to screen
  const Epetra_Comm& comm = elemagdishdg->Comm();
  if (comm.MyPID() == 0)
  {
    std::cout << "---------------------------------------------------------------------------------"
              << std::endl;
    std::cout << "---------- You are now about to enter the module for electromagnetics! ----------"
              << std::endl;
    std::cout << "---------------------------------------------------------------------------------"
              << std::endl;
  }

  // call fill complete on discretization
  if (not elemagdishdg->Filled() || not elemagdishdg->HaveDofs()) elemagdishdg->FillComplete();

  // build map
  const Teuchos::RCP<Epetra_IntVector> eledofs =
      Teuchos::rcp(new Epetra_IntVector(*elemagdishdg->ElementColMap()));
  // for(int i=0; i<elemagdishdg->NumMyColElements(); ++i)
  //{
  //  (*eledofs)[i] =
  //  dynamic_cast<DRT::ELEMENTS::Elemag*>(elemagdishdg->lColElement(i))->NumDofPerElementAuxiliary();
  //}

  Teuchos::RCP<DRT::DofSetInterface> dofsetaux =
      Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(0, eledofs, 0, false));
  elemagdishdg->AddDofSet(dofsetaux);

  // call fill complete on discretization
  elemagdishdg->FillComplete();

  // create solver
  const int linsolvernumber_elemag = elemagparams.get<int>("LINEAR_SOLVER");
  if (linsolvernumber_elemag == (-1))
    dserror(
        "There is not any linear solver defined for electromagnetic problem. Please set "
        "LINEAR_SOLVER in ELECTROMAGNETIC DYNAMIC to a valid number!");
  Teuchos::RCP<LINALG::Solver> solver = Teuchos::rcp(new LINALG::Solver(
      problem->SolverParams(linsolvernumber_elemag), comm, problem->ErrorFile()->Handle()));

  // declare output writer
  Teuchos::RCP<IO::DiscretizationWriter> output = elemagdishdg->Writer();

  // declare electromagnetic parameter list
  Teuchos::RCP<Teuchos::ParameterList> params =
      Teuchos::rcp(new Teuchos::ParameterList(elemagparams));
  // params->set<bool>("writeelemagoutput",true);

  // set restart step if required
  int restart = problem->Restart();
  params->set<int>("restart", restart);

  // create algorithm depending on time-integration scheme
  INPAR::ELEMAG::DynamicType elemagdyna =
      DRT::INPUT::IntegralValue<INPAR::ELEMAG::DynamicType>(elemagparams, "TIMEINT");
  Teuchos::RCP<ELEMAG::ElemagTimeInt> elemagalgo;
  switch (elemagdyna)
  {
    case INPAR::ELEMAG::elemag_ost:
    {
      // elemagalgo = Teuchos::rcp(new ELEMAG::TimIntOST(elemagdishdg,solver,params,output));
      elemagalgo = Teuchos::rcp(new ELEMAG::ElemagTimeInt(elemagdishdg, solver, params, output));
      break;
    }
    default:
      dserror("Unknown time-integration scheme for problem type electromagnetics");
      break;
  }

  // print information to screen
  elemagalgo->PrintInformationToScreen();

  // set initial field
  // if (restart) elemagalgo->ReadRestart(restart);
  // else
  //{
  int startfuncno = elemagparams.get<int>("STARTFUNCNO");
  elemagalgo->SetInitialField(startfuncno);
  //}

  // call time-integration scheme
  elemagalgo->Integrate();

  // print computing time
  Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm = COMM_UTILS::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, true);

  // do result test if required
  // problem->AddFieldTest(elemagalgo->CreateFieldTest());
  problem->TestAll(comm);

  return;
}
