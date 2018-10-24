/*----------------------------------------------------------------------*/
/*!
\file elemag_dyn.cpp

\brief Main control routine for electromagnetic simulations

<pre>
\level 3

\maintainer Luca Berardocco
            berardocco@lnm.mw.tum.de
            089 - 289-15244
</pre>
*/
/*----------------------------------------------------------------------*/

#include <iostream>

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "elemag_dyn.H"
#include "elemag_timeint.H"
#include "elemag_ele.H"
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

  // The function NumDofPerElementAuxiliary() of the electromagnetic elements return nsd_*2. This
  // does not assure that the code will work in any case (more spatial dimensions might give
  // problems)
  if (problem->NDim() != 3)
  {
    dserror(
        "The implementation of electromagnetic propagation only supports 3D problems.\n"
        "It is necessary to change the spatial dimension of your problem.");
  }

  // declare problem-specific parameter list for electromagnetics
  const Teuchos::ParameterList& elemagparams = problem->ElectromagneticParams();

  // declare discretization and check their existence
  Teuchos::RCP<DRT::DiscretizationHDG> elemagdishdg =
      Teuchos::rcp_dynamic_cast<DRT::DiscretizationHDG>(problem->GetDis("elemag"));
  if (elemagdishdg == Teuchos::null)
    dserror("Failed to cast DRT::Discretization to DRT::DiscretizationHDG.");

#ifdef DEBUG
  elemagdishdg->PrintFaces(std::cout);
#endif

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
  // Asking the discretization how many internal DOF the elements have and creating the additional
  // DofSet
  int eledofs = dynamic_cast<DRT::ELEMENTS::Elemag*>(elemagdishdg->lColElement(0))
                    ->NumDofPerElementAuxiliary();
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
      dserror("One step theta not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::TimIntOST(elemagdishdg,solver,params,output));
      break;
    }
    case INPAR::ELEMAG::elemag_implicit_euler:
    {
      elemagalgo = Teuchos::rcp(new ELEMAG::ElemagTimeInt(elemagdishdg, solver, params, output));
      break;
    }
    case INPAR::ELEMAG::elemag_explicit_euler:
    {
      dserror("Explicit euler method not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::TimeIntExplEuler(elemagdishdg,solver,params,output));
      break;
    }
    case INPAR::ELEMAG::elemag_rk:
    {
      dserror("Runge-Kutta methods not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::TimeIntRK(elemagdishdg,solver,params,output));
      break;
    }
    case INPAR::ELEMAG::elemag_cn:
    {
      dserror("Crank-Nicolson method not yet implemented.");
      // elemagalgo = Teuchos::rcp(new ELEMAG::TimeIntCN(elemagdishdg,solver,params,output));
      break;
    }
    default:
      dserror("Unknown time-integration scheme for problem type electromagnetics");
      break;
  }

  // Initialize the evolution algorithm
  elemagalgo->Init();

  // print information to screen
  elemagalgo->PrintInformationToScreen();

  // set initial field
  if (restart)
    elemagalgo->ReadRestart(restart);
  else
  {
    int startfuncno = elemagparams.get<int>("STARTFUNCNO");
    INPAR::ELEMAG::InitialField init =
        DRT::INPUT::IntegralValue<INPAR::ELEMAG::InitialField>(elemagparams, "INITIALFIELD");
    elemagalgo->SetInitialField(init, startfuncno);
  }

  // call time-integration scheme
  elemagalgo->Integrate();

  // Computing the error at the las time step (the conditional stateme nt is inside for now)
  if (DRT::INPUT::IntegralValue<bool>(elemagparams, "CALCERR"))
  {
    Teuchos::RCP<Epetra_SerialDenseVector> errors = elemagalgo->ComputeError();
    if (comm.MyPID() == 0)
    {
      std::cout
          << "-----------------------------------------------------------------------------------"
          << std::endl;
      std::cout << "---------------------------- Result error wrt FUNCT "
                << elemagparams.get<int>("ERRORFUNCNO") << " -----------------------------"
                << std::endl;
      std::cout << "Electric L2-error: " << sqrt((*errors)[0]) << std::endl;
      std::cout << "Magnetic L2-error: " << sqrt((*errors)[2]) << std::endl;
      std::cout
          << "-----------------------------------------------------------------------------------"
          << std::endl;
    }
  }

  // print computing time
  Teuchos::RCP<const Teuchos::Comm<int>> TeuchosComm = COMM_UTILS::toTeuchosComm<int>(comm);
  Teuchos::TimeMonitor::summarize(TeuchosComm.ptr(), std::cout, false, true, true);

  // do result test if required
  problem->AddFieldTest(elemagalgo->CreateFieldTest());
  problem->TestAll(comm);

  return;
}
