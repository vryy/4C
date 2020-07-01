/*----------------------------------------------------------------------*/
/*! \file
\brief Class for performing UQ with redairway problems


\level 2
*/
/*----------------------------------------------------------------------*/

#ifdef HAVE_FFTW

#include "../drt_io/io.H"
#include "../drt_lib/drt_globalproblem.H"
#include "drt_uq_redairways.H"
#include "../drt_io/io_pstream.H"
#include "randomfield.H"
#include "randomfield_fourier.H"
#include "randomfield_spectral.H"
#include "mc_mat_par_manager.H"
#include "../drt_io/io_control.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils_sparse_algebra_math.H"
#include "../drt_comm/comm_utils.H"
#include "../drt_red_airways/airwayimplicitintegration.H"
#include "../drt_red_airways/red_airway_resulttest.H"
#include "../drt_mat/maxwell_0d_acinus_Ogden.H"
#include "../drt_lib/drt_colors.H"

/*----------------------------------------------------------------------*/
/* standard constructor */
UQ::UQ_REDAIRWAYS::UQ_REDAIRWAYS(Teuchos::RCP<DRT::Discretization> dis) : discret_(dis)
{
  // set degrees of freedom in the discretization
  if (not discret_->Filled() || not discret_->HaveDofs()) discret_->FillComplete();

  filename_ = DRT::Problem::Instance()->OutputControlFile()->FileName();

  // input parameters multi level monte carlo
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

  // get starting random seed
  start_random_seed_ = mlmcp.get<int>("INITRANDOMSEED");

  // controlling parameter
  start_run_ = mlmcp.get<int>("START_RUN");
  int numruns = mlmcp.get<int>("NUMRUNS");
  DRT::Problem* problem = DRT::Problem::Instance();
  Teuchos::RCP<Epetra_Comm> lcomm = problem->GetNPGroup()->LocalComm();
  int NNestedGroups = problem->GetNPGroup()->NumGroups();
  int i = problem->GetNPGroup()->GroupId();

  numruns_pergroup_ = int(ceil(numruns / NNestedGroups));
  start_run_ += (i)*numruns_pergroup_;

  numb_run_ = start_run_;  // counter of how many runs were made monte carlo

  reduced_output_ = DRT::INPUT::IntegralValue<int>(mlmcp, "REDUCED_OUTPUT");

  // set up managers to create stochasticity
  my_matpar_manager_ = Teuchos::rcp(new UQ::MCMatParManager(discret_));
  my_matpar_manager_->SetupRandomFields(2);
}

/*----------------------------------------------------------------------*/
/* analyse */
void UQ::UQ_REDAIRWAYS::Integrate()
{
  const int myrank = discret_->Comm().MyPID();
  ///------
  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  const Teuchos::ParameterList& rawdyn = DRT::Problem::Instance()->ReducedDAirwayDynamicParams();
  const Teuchos::ParameterList& mlmcp = DRT::Problem::Instance()->MultiLevelMonteCarloParams();

  // nested par
  int numruns = numruns_pergroup_ + start_run_;

  // get initial random seed from inputfile
  unsigned int random_seed = mlmcp.get<int>("INITRANDOMSEED");

  const double t0 = Teuchos::Time::wallTime();
  do
  {
    if (myrank == 0)
    {
      IO::cout << "================================================================================"
               << IO::endl;
      IO::cout << "                            UQ USING MONTE CARLO                              "
               << IO::endl;
      IO::cout << "                              RUN: " << numb_run_ << "  of  " << numruns
               << IO::endl;
      IO::cout << "================================================================================"
               << IO::endl;
    }
    const double t1 = Teuchos::Time::wallTime();

    my_matpar_manager_->SetUpStochMats((random_seed + (unsigned int)numb_run_), 1.0, false);
    const double t2 = Teuchos::Time::wallTime();
    discret_->Comm().Barrier();

    discret_->Writer()->NewResultFile(filename_, (numb_run_));
    // -------------------------------------------------------------------
    // context for output and restart
    // -------------------------------------------------------------------
    Teuchos::RCP<IO::DiscretizationWriter> output = discret_->Writer();
    output->WriteMesh(0, 0.0);

    // store time
    double t3 = 0;
    double t4 = 0;

    // -------------------------------------------------------------------
    // create a solver
    // -------------------------------------------------------------------
    // get the solver number
    const int linsolvernumber = rawdyn.get<int>("LINEAR_SOLVER");
    // check if the present solver has a valid solver number
    if (linsolvernumber == (-1))
      dserror(
          "no linear solver defined. Please set LINEAR_SOLVER in REDUCED DIMENSIONAL AIRWAYS "
          "DYNAMIC to a valid number!");
    Teuchos::RCP<LINALG::Solver> solver =
        Teuchos::rcp(new LINALG::Solver(DRT::Problem::Instance()->SolverParams(linsolvernumber),
                         discret_->Comm(), DRT::Problem::Instance()->ErrorFile()->Handle()),
            false);
    discret_->ComputeNullSpaceIfNecessary(solver->Params());

    // -------------------------------------------------------------------
    // set parameters in list required for all schemes
    // -------------------------------------------------------------------
    Teuchos::ParameterList airwaystimeparams;

    // -------------------------------------- number of degrees of freedom
    // number of degrees of freedom
    const int ndim = DRT::Problem::Instance()->NDim();
    airwaystimeparams.set<int>("number of degrees of freedom", 1 * ndim);

    // -------------------------------------------------- time integration
    // the default time step size
    airwaystimeparams.set<double>("time step size", rawdyn.get<double>("TIMESTEP"));
    // maximum number of timesteps
    airwaystimeparams.set<int>("max number timesteps", rawdyn.get<int>("NUMSTEP"));

    // ----------------------------------------------- restart and output
    // restart
    airwaystimeparams.set("write restart every", rawdyn.get<int>("RESTARTEVRY"));
    // solution output
    airwaystimeparams.set("write solution every", rawdyn.get<int>("RESULTSEVRY"));

    // ----------------------------------------------- solver parameters
    // solver type
    airwaystimeparams.set("solver type", rawdyn.get<std::string>("SOLVERTYPE"));
    // tolerance
    airwaystimeparams.set("tolerance", rawdyn.get<double>("TOLERANCE"));
    // maximum number of iterations
    airwaystimeparams.set("maximum iteration steps", rawdyn.get<int>("MAXITERATIONS"));
    // solveScatra
    if (rawdyn.get<std::string>("SOLVESCATRA") == "yes")
      airwaystimeparams.set("SolveScatra", true);
    else
      airwaystimeparams.set("SolveScatra", false);

    // compute Interdependency
    if (rawdyn.get<std::string>("COMPAWACINTER") == "yes")
      airwaystimeparams.set("CompAwAcInter", true);
    else
      airwaystimeparams.set("CompAwAcInter", false);

    // prestress parameters
    // Adjust acini volume with pre-stress condition
    if (rawdyn.get<std::string>("CALCV0PRESTRESS") == "yes")
      airwaystimeparams.set("CalcV0PreStress", true);
    else
      airwaystimeparams.set("CalcV0PreStress", false);

    airwaystimeparams.set("transpulmpress", rawdyn.get<double>("TRANSPULMPRESS"));

    //------------------------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor);
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    //------------------------------------------------------------------

    Teuchos::RCP<AIRWAY::RedAirwayImplicitTimeInt> airwayimplicit = Teuchos::rcp(
        new AIRWAY::RedAirwayImplicitTimeInt(discret_, *solver, airwaystimeparams, *output));

    // initial field from restart or calculated by given function
    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      dserror("Restart not possible in UQ");
      // read the restart information, set vectors and variables
      airwayimplicit->ReadRestart(restart);
    }

    airwaystimeparams.set<FILE*>("err file", DRT::Problem::Instance()->ErrorFile()->Handle());

    t3 = Teuchos::Time::wallTime();

    // call time-integration (or stationary) scheme
    Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
    params->set<int>("run", numb_run_);
    airwayimplicit->Integrate(false, params);

    t4 = Teuchos::Time::wallTime();

    if (!reduced_output_)
    {
      Teuchos::RCP<DRT::ResultTest> resulttest =
          Teuchos::rcp(new AIRWAY::RedAirwayResultTest(*airwayimplicit));

      DRT::Problem::Instance()->AddFieldTest(resulttest);
      DRT::Problem::Instance()->TestAll(discret_->Comm());
    }

    numb_run_++;
    const double t5 = Teuchos::Time::wallTime();
    if (!reduced_output_)
    {
      // plot some time measurements
      if (!discret_->Comm().MyPID())
      {
        IO::cout << "\n=================Time  Measurement================" << IO::endl;
        IO::cout << "Setup Stochastic Material:\t" << std::setprecision(4) << t2 - t1 << "\ts"
                 << IO::endl;
        IO::cout << "Forward Solve  :\t" << std::setprecision(4) << t4 - t3 << "\ts" << IO::endl;
        IO::cout << "Total Wall Time After " << numb_run_ << " runs:\t" << std::setprecision(4)
                 << t5 - t0 << "\ts" << IO::endl;
        IO::cout << "==================================================" << IO::endl;
      }
    }
  } while (numb_run_ < numruns);
  const double t6 = Teuchos::Time::wallTime();
  IO::cout << "\n=================Time  Measurement================" << IO::endl;
  IO::cout << "Total runtime:\t" << std::setprecision(4) << t6 - t0 << "\ts" << IO::endl;
  IO::cout << "==================================================" << IO::endl;
  return;
}

#endif  // HAVE_FFTW
