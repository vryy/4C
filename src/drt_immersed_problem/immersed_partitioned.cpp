/*----------------------------------------------------------------------*/
/*! \file

\brief base class for all partitioned immersed algorithms

\level 2

\maintainer Jonas Eichinger
*----------------------------------------------------------------------*/
#include "immersed_base.H"
#include "immersed_partitioned.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"

#include "../drt_inpar/inpar_cell.H"

#include "../drt_fsi/fsi_nox_aitken.H"
#include "../drt_fsi/fsi_nox_fixpoint.H"

#include <Epetra_Time.h>
#include <Teuchos_TimeMonitor.hpp>


IMMERSED::ImmersedPartitioned::ImmersedPartitioned(const Epetra_Comm& comm)
    : ImmersedBase(),
      ADAPTER::AlgorithmBase(comm, DRT::Problem::Instance()->CellMigrationParams()),
      counter_(7)
{
  // keep constructor empty
}  // ImmersedPartitioned constructor

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitioned::Timeloop(
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
  const Teuchos::ParameterList& immerseddyn = DRT::Problem::Instance()->ImmersedMethodParams();

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = noxparameterlist_;

  // sublists
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method", "Newton"));
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");

  // Create printing utilities
  utils_ = Teuchos::rcp(new NOX::Utils(printParams));

  // ==================================================================

  // log solver iterations

  Teuchos::RCP<std::ofstream> log;
  if (Comm().MyPID() == 0)
  {
    std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
    s.append(".iteration");
    log = Teuchos::rcp(new std::ofstream(s.c_str()));
    (*log) << "# num procs      = " << Comm().NumProc() << "\n"
           << "# Method         = " << nlParams.sublist("Direction").get("Method", "Newton") << "\n"
           << "# Jacobian       = " << nlParams.get("Jacobian", "None") << "\n"
           << "# Preconditioner = " << nlParams.get("Preconditioner", "None") << "\n"
           << "# Line Search    = " << nlParams.sublist("Line Search").get("Method", "Aitken")
           << "\n"
           << "# Predictor      = '"
           << immerseddyn.sublist("PARTITIONED SOLVER").get<std::string>("PREDICTOR") << "'\n"
           << "#\n"
           << "# step | time | time/step | #nliter  |R|  #liter  Residual  Jac  Prec  FD_Res  "
              "MF_Res  MF_Jac  User\n";
  }

  // construct the timer
  Teuchos::Time timer("time step timer");

  // ==================================================================

  while (NotFinished())
  {
    // set time step size for this step
    SetFieldDt();

    // Increment all field counters and predict field values whenever
    // appropriate.
    PrepareTimeStep();

    // reset all counters
    std::fill(counter_.begin(), counter_.end(), 0);
    lsParams.sublist("Output").set("Total Number of Linear Iterations", 0);
    linsolvcount_.resize(0);

    // start time measurement
    Teuchos::RCP<Teuchos::TimeMonitor> timemonitor =
        Teuchos::rcp(new Teuchos::TimeMonitor(timer, true));

    /*----------------- CSD - predictor for itnum==0 --------------------*/

    // Begin Nonlinear Solver ************************************

    // Get initial guess.
    Teuchos::RCP<Epetra_Vector> soln = InitialGuess();

    NOX::Epetra::Vector noxSoln(soln, NOX::Epetra::Vector::CreateView);

    // Create the linear system
    Teuchos::RCP<NOX::Epetra::LinearSystem> linSys =
        CreateLinearSystem(nlParams, interface, noxSoln, utils_);

    // Create the Group
    Teuchos::RCP<NOX::Epetra::Group> grp =
        Teuchos::rcp(new NOX::Epetra::Group(printParams, interface, noxSoln, linSys));

    // Convergence Tests
    Teuchos::RCP<NOX::StatusTest::Combo> combo = CreateStatusTest(nlParams, grp);

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(
        grp, combo, Teuchos::RCP<Teuchos::ParameterList>(&nlParams, false));

    // solve the whole thing
    NOX::StatusTest::StatusType status = solver->solve();

    if (status != NOX::StatusTest::Converged) dserror("Nonlinear solver failed to converge!");

    // End Nonlinear Solver **************************************

    // Output the parameter list
    if (utils_->isPrintType(NOX::Utils::Parameters))
      if (Step() == 1 and Comm().MyPID() == 0)
      {
        utils_->out() << std::endl
                      << "Final Parameters" << std::endl
                      << "****************" << std::endl;
        solver->getList().print(utils_->out());
        utils_->out() << std::endl;
      }

    // ==================================================================

    // stop time measurement
    timemonitor = Teuchos::null;

    if (Comm().MyPID() == 0)
    {
      (*log) << Step() << "\t" << Time() << "\t" << timer.totalElapsedTime() << "\t"
             << nlParams.sublist("Output").get("Nonlinear Iterations", 0) << "\t"
             << nlParams.sublist("Output").get("2-Norm of Residual", 0.) << "\t"
             << lsParams.sublist("Output").get("Total Number of Linear Iterations", 0);
      for (std::vector<int>::size_type i = 0; i < counter_.size(); ++i)
      {
        (*log) << " " << counter_[i];
      }
      (*log) << std::endl;
      log->flush();
    }

    // ==================================================================


    // calculate stresses, strains, energies
    PrepareOutput();

    // prepare field variables for new time step
    Update();

    // write current solution
    Output();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitioned::DoStep(
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface,
    const Teuchos::RCP<Epetra_Vector>& coupling_info)
{
  PrintStepInfo();

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = noxparameterlist_;

  // sublists
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method", "Newton"));
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");

  // Create printing utilities
  utils_ = Teuchos::rcp(new NOX::Utils(printParams));

  // construct the timer
  Teuchos::Time timer("time step timer");

  // set time step size for this step
  SetFieldDt();

  // increment counters
  IncrementTimeAndStep();

  // print info
  PrintHeader();

  // reset all counters
  std::fill(counter_.begin(), counter_.end(), 0);
  lsParams.sublist("Output").set("Total Number of Linear Iterations", 0);
  linsolvcount_.resize(0);

  // start time measurement
  Teuchos::RCP<Teuchos::TimeMonitor> timemonitor =
      Teuchos::rcp(new Teuchos::TimeMonitor(timer, true));

  /*----------------- CSD - predictor for itnum==0 --------------------*/

  // Begin Nonlinear Solver ************************************

  // Get initial guess.
  Teuchos::RCP<Epetra_Vector> soln = InitialGuess();

  NOX::Epetra::Vector noxSoln(soln, NOX::Epetra::Vector::CreateView);

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys =
      CreateLinearSystem(nlParams, interface, noxSoln, utils_);

  // Create the Group
  Teuchos::RCP<NOX::Epetra::Group> grp =
      Teuchos::rcp(new NOX::Epetra::Group(printParams, interface, noxSoln, linSys));

  // Convergence Tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo = CreateStatusTest(nlParams, grp);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver =
      NOX::Solver::buildSolver(grp, combo, Teuchos::RCP<Teuchos::ParameterList>(&nlParams, false));

  // solve the whole thing
  NOX::StatusTest::StatusType status = solver->solve();

  if (status != NOX::StatusTest::Converged) dserror("Nonlinear solver failed to converge!");

  // End Nonlinear Solver **************************************

  // Output the parameter list
  if (utils_->isPrintType(NOX::Utils::Parameters))
    if (Step() == 1 and Comm().MyPID() == 0)
    {
      utils_->out() << std::endl
                    << "Final Parameters" << std::endl
                    << "****************" << std::endl;
      solver->getList().print(utils_->out());
      utils_->out() << std::endl;
    }

  // return the inter-module coupling information
  return ReturnCouplingInfo();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool IMMERSED::ImmersedPartitioned::computeF(
    const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  const char* flags[] = {"Residual", "Jac", "Prec", "FD_Res", "MF_Res", "MF_Jac", "User", NULL};

  Epetra_Time timer(x.Comm());
  const double startTime = timer.WallTime();

  if (Comm().MyPID() == 0)
  {
    utils_->out() << "\n "
                  << "Global residual calculation"
                  << ".\n";
    if (fillFlag != Residual) utils_->out() << " fillFlag = " << flags[fillFlag] << "\n";
  }

  // we count the number of times the residuum is build
  counter_[fillFlag] += 1;

  if (!x.Map().UniqueGIDs()) dserror("source map not unique");


  // Do the coupling step. The real work is in here.
  CouplingOp(x, F, fillFlag);


  const double endTime = timer.WallTime();
  if (Comm().MyPID() == 0)
    utils_->out() << "\nTime for residual calculation: " << endTime - startTime << " secs\n\n";
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem> IMMERSED::ImmersedPartitioned::CreateLinearSystem(
    Teuchos::ParameterList& nlParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface, NOX::Epetra::Vector& noxSoln,
    Teuchos::RCP<NOX::Utils> utils)
{
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method", "Aitken"));
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");


  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac;

  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys;

  // ==================================================================
  // decide on Jacobian and preconditioner
  // We might want to use no preconditioner at all. Some kind of
  // Jacobian has to be provided, otherwise the linear system uses
  // plain finite differences.

  const std::string jacobian = nlParams.get("Jacobian", "None");
  std::string preconditioner = nlParams.get("Preconditioner", "None");

  // No Jacobian at all. Do a fix point iteration.
  if (jacobian == "None")
  {
    preconditioner = "None";
  }
  else
  {
    dserror("unsupported Jacobian '%s'", jacobian.c_str());
  }

  // ==================================================================

  // No preconditioning at all.
  if (preconditioner == "None")
  {
    if (Teuchos::is_null(iJac))
    {
      // if no Jacobian has been set this better be the fix point
      // method.
      if (dirParams.get("Method", "Newton") != "User Defined")
      {
        if (Comm().MyPID() == 0)
          utils->out() << "Warning: No Jacobian for solver " << dirParams.get("Method", "Newton")
                       << "\n";
      }
      linSys = Teuchos::rcp(
          new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, interface, noxSoln));
    }
    else
    {
      dserror("only fix point interation supported");
    }
  }
  else
  {
    dserror("unsupported preconditioner '%s'", preconditioner.c_str());
  }

  return linSys;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo> IMMERSED::ImmersedPartitioned::CreateStatusTest(
    Teuchos::ParameterList& nlParams, Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // setup the real tests
  CreateStatusTest(nlParams, grp, converged);

  return combo;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitioned::CreateStatusTest(Teuchos::ParameterList& nlParams,
    Teuchos::RCP<NOX::Epetra::Group> grp, Teuchos::RCP<NOX::StatusTest::Combo> converged)
{
  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
      Teuchos::rcp(new NOX::StatusTest::NormF(nlParams.get("Norm abs F", 1.0e-6)));
  converged->addStatusTest(absresid);

  if (nlParams.isParameter("Norm Update"))
  {
    Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
        Teuchos::rcp(new NOX::StatusTest::NormUpdate(nlParams.get("Norm Update", 1.0e-5)));
    converged->addStatusTest(update);
  }

  if (nlParams.isParameter("Norm rel F"))
  {
    Teuchos::RCP<NOX::StatusTest::NormF> relresid =
        Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), nlParams.get("Norm rel F", 1.0e-2)));
    converged->addStatusTest(relresid);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitioned::PrepareOutput()
{
  dserror("not implemented in this class. my be overridden by derived class.");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitioned::CouplingOp(
    const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  dserror("not implemented in this class. my be overridden by derived class.");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitioned::SetDefaultParameters(
    const Teuchos::ParameterList& immersedpart, Teuchos::ParameterList& list)
{
  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = list;

  nlParams.set("Nonlinear Solver", "Line Search Based");
  nlParams.set("Preconditioner", "None");
  nlParams.set("Norm abs F", immersedpart.get<double>("CONVTOL"));
  nlParams.set("Max Iterations", immersedpart.get<int>("ITEMAX"));

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& lineSearchParams = nlParams.sublist("Line Search");

  //
  // Set parameters for NOX to chose the solver direction and line
  // search step.
  //
  switch (DRT::INPUT::IntegralValue<int>(immersedpart, "COUPALGO"))
  {
    case INPAR::CELL::cell_iter_stagg_fixed_rel_param:
    {
      // fixed-point solver with fixed relaxation parameter
      SetMethod("ITERATIVE STAGGERED SCHEME WITH FIXED RELAXATION PARAMETER");

      nlParams.set("Jacobian", "None");

      dirParams.set("Method", "User Defined");
      Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
          Teuchos::rcp(new NOX::FSI::FixPointFactory());
      dirParams.set("User Defined Direction Factory", fixpointfactory);

      lineSearchParams.set("Method", "Full Step");
      lineSearchParams.sublist("Full Step").set("Full Step", immersedpart.get<double>("RELAX"));
      break;
    }
    case INPAR::CELL::cell_iter_stagg_AITKEN_rel_param:
    {
      // fixed-point solver with Aitken relaxation parameter
      SetMethod("ITERATIVE STAGGERED SCHEME WITH RELAXATION PARAMETER VIA AITKEN ITERATION");

      nlParams.set("Jacobian", "None");

      dirParams.set("Method", "User Defined");
      Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
          Teuchos::rcp(new NOX::FSI::FixPointFactory());
      dirParams.set("User Defined Direction Factory", fixpointfactory);

      Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> aitkenfactory =
          Teuchos::rcp(new NOX::FSI::AitkenFactory());
      lineSearchParams.set("Method", "User Defined");
      lineSearchParams.set("User Defined Line Search Factory", aitkenfactory);

      lineSearchParams.sublist("Aitken").set("max step size", immersedpart.get<double>("MAXOMEGA"));
      break;
    }
    case INPAR::CELL::cell_basic_sequ_stagg:
    {
      // sequential coupling (no iteration!)
      SetMethod("BASIC SEQUENTIAL STAGGERED SCHEME");

      nlParams.set("Jacobian", "None");
      nlParams.set("Max Iterations", 1);

      dirParams.set("Method", "User Defined");
      Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
          Teuchos::rcp(new NOX::FSI::FixPointFactory());
      dirParams.set("User Defined Direction Factory", fixpointfactory);

      lineSearchParams.set("Method", "Full Step");
      lineSearchParams.sublist("Full Step").set("Full Step", 1.0);
      break;
    }
    default:
    {
      dserror("coupling method type '%s' unsupported",
          immersedpart.get<std::string>("COUPALGO").c_str());
      break;
    }
  }

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm().MyPID());

  // set default output flag to no output
  // The field solver will output a lot, anyway.
  printParams.get("Output Information",
      ::NOX::Utils::Warning | ::NOX::Utils::OuterIteration | ::NOX::Utils::OuterIterationStatusTest
      // ::NOX::Utils::Parameters
  );

  Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
  solverOptions.set<std::string>("Status Test Check Type", "Complete");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void IMMERSED::ImmersedPartitioned::BackgroundOp(
    Teuchos::RCP<Epetra_Vector> backgrd_dirichlet_values, const FillType fillFlag)
{
  if (Comm().MyPID() == 0 and utils_->isPrintType(NOX::Utils::OuterIteration))
    utils_->out() << std::endl << "Background operator" << std::endl;
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> IMMERSED::ImmersedPartitioned::ImmersedOp(
    Teuchos::RCP<Epetra_Vector> bdry_traction, const FillType fillFlag)
{
  if (Comm().MyPID() == 0 and utils_->isPrintType(NOX::Utils::OuterIteration))
    utils_->out() << std::endl << "Immersed operator" << std::endl;
  return Teuchos::null;
}
