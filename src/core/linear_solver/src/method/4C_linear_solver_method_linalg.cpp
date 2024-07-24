/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation

\level 0

*-----------------------------------------------------------------------*/
#include "4C_linear_solver_method_linalg.hpp"

#include "4C_io_pstream.hpp"
#include "4C_linalg_sparsematrix.hpp"
#include "4C_linear_solver_method_direct.hpp"
#include "4C_linear_solver_method_iterative.hpp"
#include "4C_utils_parameter_list.hpp"

#include <BelosTypes.hpp>  // for Belos verbosity codes
#include <Epetra_LinearProblem.h>
#include <Epetra_MpiComm.h>
#include <ml_MultiLevelPreconditioner.h>  // includes for ML parameter list validation
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::LinAlg::Solver::Solver(const Teuchos::ParameterList& inparams, const Epetra_Comm& comm,
    const std::function<const Teuchos::ParameterList&(int)>& get_solver_params,
    Core::IO::Verbositylevel verbosity, const bool translate_params_to_belos)
    : comm_(comm), params_(Teuchos::rcp(new Teuchos::ParameterList()))
{
  if (translate_params_to_belos)
    *params_ = translate_solver_parameters(inparams, get_solver_params, verbosity);
  else
    *params_ = inparams;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Core::LinAlg::Solver::~Solver() { reset(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::Solver::reset() { solver_ = Teuchos::null; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::LinAlg::Solver::get_num_iters() const { return solver_->get_num_iters(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::Solver::adapt_tolerance(
    const double desirednlnres, const double currentnlnres, const double better)
{
  if (!params().isSublist("Belos Parameters")) FOUR_C_THROW("Adaptive tolerance only for Belos.");

  Teuchos::ParameterList& solver_params = params().sublist("Belos Parameters");

  if (!solver_params.isParameter("Convergence Tolerance"))
    FOUR_C_THROW("No iterative solver tolerance in ParameterList");

  const bool do_output = solver_params.get<int>("Output Frequency", 1) and !get_comm().MyPID();

  const std::string conv_test_strategy = solver_params.get<std::string>(
      "Implicit Residual Scaling", Belos::convertScaleTypeToString(Belos::ScaleType::None));

  if (conv_test_strategy != Belos::convertScaleTypeToString(Belos::ScaleType::NormOfInitRes))
  {
    FOUR_C_THROW(
        "You are using an adaptive tolerance for the linear solver. Therefore, the iterative "
        "solver needs to work with a relative residual norm. This can be achieved by setting "
        "'AZCONV' to 'AZ_r0' in the input file.");
  }

  // save original value of convergence
  const bool have_saved_value = solver_params.isParameter("Convergence Tolerance Saved");
  if (!have_saved_value)
  {
    solver_params.set<double>(
        "Convergence Tolerance Saved", solver_params.get<double>("Convergence Tolerance"));
  }

  const double input_tolerance = solver_params.get<double>("Convergence Tolerance Saved");

  if (do_output)
    std::cout << "                --- Solver input relative tolerance " << input_tolerance << "\n";
  if (currentnlnres * input_tolerance < desirednlnres)
  {
    double adapted_tolerance = desirednlnres * better / currentnlnres;
    if (adapted_tolerance > 1.0)
    {
      adapted_tolerance = 1.0;
      if (do_output)
      {
        std::cout << "WARNING:  Computed adapted relative tolerance bigger than 1\n";
        std::cout << "          Value constrained to 1, but consider adapting Parameter "
                     "ADAPTCONV_BETTER\n";
      }
    }
    if (adapted_tolerance < input_tolerance) adapted_tolerance = input_tolerance;
    if (do_output && adapted_tolerance > input_tolerance)
    {
      std::cout << "                *** Solver adapted relative tolerance " << adapted_tolerance
                << "\n";
    }

    solver_params.set<double>("Convergence Tolerance", adapted_tolerance);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::Solver::set_tolerance(const double tolerance)
{
  if (!params().isSublist("Belos Parameters"))
    FOUR_C_THROW("Set tolerance of linear solver only for Belos solver.");

  Teuchos::ParameterList& solver_params = params().sublist("Belos Parameters");

  const bool have_saved_value = solver_params.isParameter("Convergence Tolerance Saved");
  if (!have_saved_value)
  {
    solver_params.set<double>(
        "Convergence Tolerance Saved", solver_params.get<double>("Convergence Tolerance", 1.0e-8));
  }

  solver_params.set<double>("Convergence Tolerance", tolerance);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::Solver::reset_tolerance()
{
  if (!params().isSublist("Belos Parameters")) return;

  Teuchos::ParameterList& solver_params = params().sublist("Belos Parameters");

  const bool have_saved_value = solver_params.isParameter("Convergence Tolerance Saved");
  if (!have_saved_value) return;

  solver_params.set<double>(
      "Convergence Tolerance", solver_params.get<double>("Convergence Tolerance Saved"));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Core::LinAlg::Solver::setup(Teuchos::RCP<Epetra_Operator> matrix,
    Teuchos::RCP<Epetra_MultiVector> x, Teuchos::RCP<Epetra_MultiVector> b,
    const SolverParams& params)
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::Solver:  1)   Setup");

  FOUR_C_ASSERT(!(params.lin_tol_better > -1.0 and params.tolerance > 0.0),
      "Do not set tolerance and adaptive tolerance to the linear solver.");

  if (params.lin_tol_better > -1.0)
  {
    adapt_tolerance(params.nonlin_tolerance, params.nonlin_residual, params.lin_tol_better);
  }

  if (params.tolerance > 0.0)
  {
    set_tolerance(params.tolerance);
  }

  // reset data flags on demand
  bool refactor = params.refactor;
  if (params.reset)
  {
    reset();
    refactor = true;
  }

  if (solver_ == Teuchos::null)
  {
    // decide what solver to use
    std::string solvertype = Solver::params().get("solver", "none");

    if ("belos" == solvertype)
    {
      solver_ =
          Teuchos::rcp(new Core::LinearSolver::IterativeSolver<Epetra_Operator, Epetra_MultiVector>(
              comm_, Solver::params()));
    }
    else if ("umfpack" == solvertype or "superlu" == solvertype)
    {
      solver_ = Teuchos::rcp(
          new Core::LinearSolver::DirectSolver<Epetra_Operator, Epetra_MultiVector>(solvertype));
    }
    else
      FOUR_C_THROW("Unknown type of solver");
  }

  solver_->setup(matrix, x, b, refactor, params.reset, params.projector);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int Core::LinAlg::Solver::solve(Teuchos::RCP<Epetra_Operator> matrix,
    Teuchos::RCP<Epetra_MultiVector> x, Teuchos::RCP<Epetra_MultiVector> b,
    const SolverParams& params)
{
  setup(matrix, x, b, params);

  int error_value = 0;
  {
    TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::Solver:  2)   Solve");
    error_value = solver_->solve();
  }

  return error_value;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
Teuchos::ParameterList translate_four_c_to_ifpack(const Teuchos::ParameterList& inparams)
{
  Teuchos::ParameterList ifpacklist;

  ifpacklist.set("fact: level-of-fill", inparams.get<int>("IFPACKGFILL"));
  ifpacklist.set("partitioner: overlap", inparams.get<int>("IFPACKOVERLAP"));
  ifpacklist.set("schwarz: combine mode",
      inparams.get<std::string>("IFPACKCOMBINE"));    // can be "Zero", "Add", "Insert"
  ifpacklist.set("schwarz: reordering type", "rcm");  // "rcm" or "metis" or "amd"

  return ifpacklist;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
Teuchos::ParameterList translate_four_c_to_ml(const Teuchos::ParameterList& inparams,
    const std::function<const Teuchos::ParameterList&(int)>& get_solver_params,
    Core::IO::Verbositylevel verbosity)
{
  Teuchos::ParameterList mllist;

  ML_Epetra::SetDefaults("SA", mllist);
  const auto prectyp =
      Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(inparams, "AZPREC");

  switch (prectyp)
  {
    case Core::LinearSolver::PreconditionerType::multigrid_ml:
      break;
    case Core::LinearSolver::PreconditionerType::multigrid_ml_fluid2:
      mllist.set("energy minimization: enable", true);
      mllist.set("energy minimization: type", 3);
      mllist.set("aggregation: block scaling", false);
      break;
    default:
      FOUR_C_THROW("Unknown type of ml preconditioner");
      break;
  }

  // set repartitioning parameters
  // En-/Disable ML repartitioning. Note: ML requires parameter to be set as integer.
  bool doRepart = Core::UTILS::IntegralValue<bool>(inparams, "ML_REBALANCE");
  if (doRepart)
  {
    mllist.set("repartition: enable", 1);

    // these are the hard-coded ML repartitioning settings
    mllist.set("repartition: partitioner", "ParMETIS");
    mllist.set("repartition: max min ratio", 1.3);
    mllist.set("repartition: min per proc", 3000);
  }
  else
  {
    mllist.set("repartition: enable", 0);
  }

  mllist.set("ML output", inparams.get<int>("ML_PRINT"));
  if (inparams.get<int>("ML_PRINT") == 10)
    mllist.set("print unused", 1);
  else
    mllist.set("print unused", -2);
  mllist.set("increasing or decreasing", "increasing");
  mllist.set("coarse: max size", inparams.get<int>("ML_MAXCOARSESIZE"));
  mllist.set("coarse: pre or post", "pre");
  mllist.set("max levels", inparams.get<int>("ML_MAXLEVEL"));
  mllist.set("smoother: pre or post", "both");
  mllist.set("aggregation: threshold", inparams.get<double>("ML_PROLONG_THRES"));
  mllist.set("aggregation: damping factor", inparams.get<double>("ML_PROLONG_SMO"));
  mllist.set("aggregation: nodes per aggregate", inparams.get<int>("ML_AGG_SIZE"));
  // override the default sweeps=2 with a default sweeps=1
  // individual level sweeps are set below
  mllist.set("smoother: sweeps", 1);
  // save memory if this is an issue, make ML use single precision
  // mllist.set("low memory usage",true);
  switch (Core::UTILS::IntegralValue<int>(inparams, "ML_COARSEN"))
  {
    case 0:
      mllist.set("aggregation: type", "Uncoupled");
      break;
    case 1:
      mllist.set("aggregation: type", "METIS");
      break;
    case 2:
      mllist.set("aggregation: type", "VBMETIS");
      break;
    case 3:
      mllist.set("aggregation: type", "MIS");
      break;
    default:
      FOUR_C_THROW("Unknown type of coarsening for ML");
      break;
  }

  // set ml smoothers
  const int mlmaxlevel = inparams.get<int>("ML_MAXLEVEL");

  // create vector of integers containing smoothing steps/polynomial order of level
  std::vector<int> mlsmotimessteps;
  {
    std::istringstream mlsmotimes(Teuchos::getNumericStringParameter(inparams, "ML_SMOTIMES"));
    std::string word;
    while (mlsmotimes >> word) mlsmotimessteps.push_back(std::atoi(word.c_str()));
  }

  if ((int)mlsmotimessteps.size() < mlmaxlevel)
    FOUR_C_THROW(
        "Not enough smoothing steps ML_SMOTIMES=%d, must be larger/equal than ML_MAXLEVEL=%d\n",
        mlsmotimessteps.size(), mlmaxlevel);

  for (int i = 0; i < mlmaxlevel - 1; ++i)
  {
    char levelstr[19];
    sprintf(levelstr, "(level %d)", i);
    Teuchos::ParameterList& smolevelsublist =
        mllist.sublist("smoother: list " + std::string(levelstr));
    int type = 0;
    double damp = 0.0;
    if (i == 0)
    {
      type = Core::UTILS::IntegralValue<int>(inparams, "ML_SMOOTHERFINE");
      damp = inparams.get<double>("ML_DAMPFINE");
    }
    else if (i < mlmaxlevel - 1)
    {
      type = Core::UTILS::IntegralValue<int>(inparams, "ML_SMOOTHERMED");
      damp = inparams.get<double>("ML_DAMPMED");
    }

    switch (type)
    {
      case 0:  // SGS
        smolevelsublist.set("smoother: type", "symmetric Gauss-Seidel");
        smolevelsublist.set("smoother: sweeps", mlsmotimessteps[i]);
        smolevelsublist.set("smoother: damping factor", damp);
        break;
      case 7:  // GS
        smolevelsublist.set("smoother: type", "Gauss-Seidel");
        smolevelsublist.set("smoother: sweeps", mlsmotimessteps[i]);
        smolevelsublist.set("smoother: damping factor", damp);
        break;
      case 2:  // Chebychev
        smolevelsublist.set("smoother: type", "MLS");
        smolevelsublist.set("smoother: sweeps", mlsmotimessteps[i]);
        break;
      case 3:  // MLS
        smolevelsublist.set("smoother: type", "MLS");
        smolevelsublist.set("smoother: MLS polynomial order", -mlsmotimessteps[i]);
        break;
      case 4:  // Ifpack's ILU
      {
        smolevelsublist.set("smoother: type", "IFPACK");
        smolevelsublist.set("smoother: ifpack type", "ILU");
        smolevelsublist.set("smoother: ifpack overlap", inparams.get<int>("IFPACKOVERLAP"));
        smolevelsublist.set<double>("smoother: ifpack level-of-fill",
            (double)mlsmotimessteps[i]);  // 12.01.2012: TW fixed double->int
        Teuchos::ParameterList& ifpacklist = mllist.sublist("smoother: ifpack list");
        ifpacklist.set("schwarz: reordering type", "rcm");  // "rcm" or "metis" or "amd" or "true"
        ifpacklist.set("schwarz: combine mode",
            inparams.get<std::string>("IFPACKCOMBINE"));  // can be "Zero", "Insert", "Add"
        ifpacklist.set("partitioner: overlap", inparams.get<int>("IFPACKOVERLAP"));
      }
      break;
      case 5:  // Amesos' KLU
        smolevelsublist.set("smoother: type", "Amesos-KLU");
        break;
      case 9:  // Amesos' Umfpack
        smolevelsublist.set("smoother: type", "Amesos-UMFPACK");
        break;
      case 6:  // Amesos' SuperLU_Dist
        smolevelsublist.set("smoother: type", "Amesos-Superludist");
        break;
      case 10:  // Braess-Sarazin smoother (only for MueLu with BlockedOperators)
      {
        smolevelsublist.set("smoother: type", "Braess-Sarazin");
        smolevelsublist.set("smoother: sweeps", mlsmotimessteps[i]);
        smolevelsublist.set("smoother: damping factor", damp);
        Teuchos::ParameterList& SchurCompList = smolevelsublist.sublist("smoother: SchurComp list");
        SchurCompList = Core::LinAlg::Solver::translate_solver_parameters(
            get_solver_params(inparams.get<int>("SUB_SOLVER1")), get_solver_params, verbosity);
      }
      break;
      case 11:  // SIMPLE smoother  (only for MueLu with BlockedOperators)
      case 12:  // SIMPLEC smoother (only for MueLu with BlockedOperators)
      {
        if (type == 11)
          smolevelsublist.set("smoother: type", "SIMPLE");
        else if (type == 12)
          smolevelsublist.set("smoother: type", "SIMPLEC");
        smolevelsublist.set("smoother: type", "SIMPLE");
        smolevelsublist.set("smoother: sweeps", mlsmotimessteps[i]);
        smolevelsublist.set("smoother: damping factor", damp);
        Teuchos::ParameterList& predictList = smolevelsublist.sublist("smoother: Predictor list");
        predictList = Core::LinAlg::Solver::translate_solver_parameters(
            get_solver_params(inparams.get<int>("SUB_SOLVER1")), get_solver_params, verbosity);
        Teuchos::ParameterList& SchurCompList = smolevelsublist.sublist("smoother: SchurComp list");
        SchurCompList = Core::LinAlg::Solver::translate_solver_parameters(
            get_solver_params(inparams.get<int>("SUB_SOLVER2")), get_solver_params, verbosity);
      }
      break;
      case 13:  // IBD: indefinite block diagonal preconditioner
      {
        smolevelsublist.set("smoother: type", "IBD");
        smolevelsublist.set("smoother: sweeps", mlsmotimessteps[i]);
        smolevelsublist.set("smoother: damping factor", damp);
        Teuchos::ParameterList& predictList = smolevelsublist.sublist("smoother: Predictor list");
        predictList = Core::LinAlg::Solver::translate_solver_parameters(
            get_solver_params(inparams.get<int>("SUB_SOLVER1")), get_solver_params, verbosity);
        Teuchos::ParameterList& SchurCompList = smolevelsublist.sublist("smoother: SchurComp list");
        SchurCompList = Core::LinAlg::Solver::translate_solver_parameters(
            get_solver_params(inparams.get<int>("SUB_SOLVER2")), get_solver_params, verbosity);
      }
      break;
      case 14:  // Uzawa: inexact Uzawa smoother
      {
        smolevelsublist.set("smoother: type", "Uzawa");
        smolevelsublist.set("smoother: sweeps", mlsmotimessteps[i]);
        smolevelsublist.set("smoother: damping factor", damp);
        Teuchos::ParameterList& predictList = smolevelsublist.sublist("smoother: Predictor list");
        predictList = Core::LinAlg::Solver::translate_solver_parameters(
            get_solver_params(inparams.get<int>("SUB_SOLVER1")), get_solver_params, verbosity);
        Teuchos::ParameterList& SchurCompList = smolevelsublist.sublist("smoother: SchurComp list");
        SchurCompList = Core::LinAlg::Solver::translate_solver_parameters(
            get_solver_params(inparams.get<int>("SUB_SOLVER2")), get_solver_params, verbosity);
      }
      break;
      default:
        FOUR_C_THROW("Unknown type of smoother for ML: tuple %d", type);
        break;
    }  // switch (type)
  }    // for (int i=0; i<azvar->mlmaxlevel-1; ++i)

  // set coarse grid solver
  const int coarse = mlmaxlevel - 1;
  switch (Core::UTILS::IntegralValue<int>(inparams, "ML_SMOOTHERCOARSE"))
  {
    case 0:
      mllist.set("coarse: type", "symmetric Gauss-Seidel");
      mllist.set("coarse: sweeps", mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor", inparams.get<double>("ML_DAMPCOARSE"));
      break;
    case 7:
      mllist.set("coarse: type", "Gauss-Seidel");
      mllist.set("coarse: sweeps", mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor", inparams.get<double>("ML_DAMPCOARSE"));
      break;
    case 2:  // Chebychev
      mllist.set("smoother: type", "MLS");
      mllist.set("smoother: sweeps", mlsmotimessteps[coarse]);
      break;
    case 3:
      mllist.set("coarse: type", "MLS");
      mllist.set("coarse: MLS polynomial order", -mlsmotimessteps[coarse]);
      break;
    case 4:
    {
      mllist.set("coarse: type", "IFPACK");
      mllist.set("coarse: ifpack type", "ILU");
      mllist.set("coarse: ifpack overlap", inparams.get<int>("IFPACKOVERLAP"));
      mllist.set<double>("coarse: ifpack level-of-fill",
          (double)mlsmotimessteps[coarse]);  // 12.01.2012: TW fixed double -> int
      Teuchos::ParameterList& ifpacklist = mllist.sublist("smoother: ifpack list");
      ifpacklist.set<int>("fact: level-of-fill", (int)mlsmotimessteps[coarse]);
      ifpacklist.set("schwarz: reordering type", "rcm");
      ifpacklist.set("schwarz: combine mode",
          inparams.get<std::string>("IFPACKCOMBINE"));  // can be "Zero", "Insert", "Add"
      ifpacklist.set("partitioner: overlap", inparams.get<int>("IFPACKOVERLAP"));
    }
    break;
    case 5:
      mllist.set("coarse: type", "Amesos-KLU");
      break;
    case 9:
      mllist.set("coarse: type", "Amesos-UMFPACK");
      break;
    case 6:
      mllist.set("coarse: type", "Amesos-Superludist");
      break;
    case 10:  // Braess-Sarazin smoother (only for MueLu with BlockedOperators)
    {
      mllist.set("coarse: type", "Braess-Sarazin");
      mllist.set("coarse: sweeps", mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor", inparams.get<double>("ML_DAMPCOARSE"));
      Teuchos::ParameterList& SchurCompList = mllist.sublist("coarse: SchurComp list");
      SchurCompList = Core::LinAlg::Solver::translate_solver_parameters(
          get_solver_params(inparams.get<int>("SUB_SOLVER2")), get_solver_params, verbosity);
    }
    break;
    case 11:  // SIMPLE smoother  (only for MueLu with BlockedOperators)
    case 12:  // SIMPLEC smoother (only for MueLu with BlockedOperators)
    {
      int type = Core::UTILS::IntegralValue<int>(inparams, "ML_SMOOTHERCOARSE");
      if (type == 11)
        mllist.set("coarse: type", "SIMPLE");
      else if (type == 12)
        mllist.set("coarse: type", "SIMPLEC");
      mllist.set("coarse: sweeps", mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor", inparams.get<double>("ML_DAMPCOARSE"));
      Teuchos::ParameterList& predictList = mllist.sublist("coarse: Predictor list");
      predictList = Core::LinAlg::Solver::translate_solver_parameters(
          get_solver_params(inparams.get<int>("SUB_SOLVER1")), get_solver_params, verbosity);
      Teuchos::ParameterList& SchurCompList = mllist.sublist("coarse: SchurComp list");
      SchurCompList = Core::LinAlg::Solver::translate_solver_parameters(
          get_solver_params(inparams.get<int>("SUB_SOLVER2")), get_solver_params, verbosity);
    }
    break;
    case 13:  // IBD: indefinite block diagonal preconditioner
    {
      mllist.set("coarse: type", "IBD");
      mllist.set("coarse: sweeps", mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor", inparams.get<double>("ML_DAMPCOARSE"));
      Teuchos::ParameterList& predictList = mllist.sublist("coarse: Predictor list");
      predictList = Core::LinAlg::Solver::translate_solver_parameters(
          get_solver_params(inparams.get<int>("SUB_SOLVER1")), get_solver_params, verbosity);
      Teuchos::ParameterList& SchurCompList = mllist.sublist("coarse: SchurComp list");
      SchurCompList = Core::LinAlg::Solver::translate_solver_parameters(
          get_solver_params(inparams.get<int>("SUB_SOLVER2")), get_solver_params, verbosity);
    }
    break;
    case 14:  // Uzawa: inexact Uzawa smoother
    {
      mllist.set("coarse: type", "Uzawa");
      mllist.set("coarse: sweeps", mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor", inparams.get<double>("ML_DAMPCOARSE"));
      Teuchos::ParameterList& predictList = mllist.sublist("coarse: Predictor list");
      predictList = Core::LinAlg::Solver::translate_solver_parameters(
          get_solver_params(inparams.get<int>("SUB_SOLVER1")), get_solver_params, verbosity);
      Teuchos::ParameterList& SchurCompList = mllist.sublist("coarse: SchurComp list");
      SchurCompList = Core::LinAlg::Solver::translate_solver_parameters(
          get_solver_params(inparams.get<int>("SUB_SOLVER2")), get_solver_params, verbosity);
    }
    break;
    default:
      FOUR_C_THROW("Unknown type of coarse solver for ML");
      break;
  }  // switch (azvar->mlsmotype_coarse)
  // default values for nullspace
  mllist.set("PDE equations", 1);
  mllist.set("null space: dimension", 1);
  mllist.set("null space: type", "pre-computed");
  mllist.set("null space: add default vectors", false);
  mllist.set<double*>("null space: vectors", nullptr);

  return mllist;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
Teuchos::ParameterList translate_four_c_to_muelu(
    const Teuchos::ParameterList& inparams, Teuchos::ParameterList* azlist)
{
  Teuchos::ParameterList muelulist;

  std::string xmlfile = inparams.get<std::string>("MUELU_XML_FILE");
  if (xmlfile != "none") muelulist.set("MUELU_XML_FILE", xmlfile);

  return muelulist;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
Teuchos::ParameterList translate_four_c_to_belos(const Teuchos::ParameterList& inparams,
    const std::function<const Teuchos::ParameterList&(int)>& get_solver_params,
    Core::IO::Verbositylevel verbosity)
{
  Teuchos::ParameterList outparams;
  outparams.set("solver", "belos");
  Teuchos::ParameterList& beloslist = outparams.sublist("Belos Parameters");

  beloslist.set("reuse", inparams.get<int>("AZREUSE"));
  beloslist.set("ncall", 0);

  // try to get an xml file if possible
  std::string xmlfile = inparams.get<std::string>("SOLVER_XML_FILE");
  if (xmlfile != "none")
  {
    beloslist.set("SOLVER_XML_FILE", xmlfile);
  }
  else
  {
    switch (verbosity)
    {
      case Core::IO::minimal:
        beloslist.set("Output Style", Belos::OutputType::Brief);
        beloslist.set("Verbosity", Belos::MsgType::Warnings);
        break;
      case Core::IO::standard:
        beloslist.set("Output Style", Belos::OutputType::Brief);
        beloslist.set("Verbosity", Belos::MsgType::Warnings + Belos::MsgType::StatusTestDetails);
        break;
      case Core::IO::verbose:
        beloslist.set("Output Style", Belos::OutputType::General);
        beloslist.set("Verbosity", Belos::MsgType::Warnings + Belos::MsgType::StatusTestDetails +
                                       Belos::MsgType::FinalSummary);
        break;
      case Core::IO::debug:
        beloslist.set("Output Style", Belos::OutputType::General);
        beloslist.set("Verbosity", Belos::MsgType::Debug);
      default:
        break;
    }
    beloslist.set("Output Frequency", inparams.get<int>("AZOUTPUT"));

    // set tolerances and iterations
    beloslist.set("Maximum Iterations", inparams.get<int>("AZITER"));
    beloslist.set("Convergence Tolerance", inparams.get<double>("AZTOL"));
    beloslist.set("Implicit Residual Scaling",
        Belos::convertScaleTypeToString(
            Teuchos::getIntegralValue<Belos::ScaleType>(inparams, "AZCONV")));

    // set type of solver
    switch (Teuchos::getIntegralValue<Core::LinearSolver::IterativeSolverType>(inparams, "AZSOLVE"))
    {
      case Core::LinearSolver::IterativeSolverType::cg:
        beloslist.set("Solver Type", "CG");
        break;
      case Core::LinearSolver::IterativeSolverType::bicgstab:
        beloslist.set("Solver Type", "BiCGSTAB");
        break;
      case Core::LinearSolver::IterativeSolverType::gmres:
        beloslist.set("Solver Type", "GMRES");
        beloslist.set("Num Blocks", inparams.get<int>("AZSUB"));
        break;
      default:
      {
        FOUR_C_THROW("Flag '%s'! \nUnknown solver for Belos.",
            Teuchos::getIntegralValue<Core::LinearSolver::IterativeSolverType>(
                inparams, "AZSOLVE"));
        break;
      }
    }
  }

  // set type of preconditioner
  const auto azprectyp =
      Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(inparams, "AZPREC");

  switch (azprectyp)
  {
    case Core::LinearSolver::PreconditionerType::ilu:
      beloslist.set("Preconditioner Type", "ILU");
      break;
    case Core::LinearSolver::PreconditionerType::multigrid_ml:
    case Core::LinearSolver::PreconditionerType::multigrid_ml_fluid2:
    case Core::LinearSolver::PreconditionerType::multigrid_muelu:
      beloslist.set("Preconditioner Type", "ML");
      break;
    case Core::LinearSolver::PreconditionerType::multigrid_muelu_contactsp:
      beloslist.set("Preconditioner Type", "ContactSP");
      break;
    case Core::LinearSolver::PreconditionerType::multigrid_muelu_beamsolid:
      beloslist.set("Preconditioner Type", "BeamSolid");
      break;
    case Core::LinearSolver::PreconditionerType::multigrid_nxn:
      beloslist.set("Preconditioner Type", "AMGnxn");
      break;
    case Core::LinearSolver::PreconditionerType::block_gauss_seidel_2x2:
      beloslist.set("Preconditioner Type", "ML");
      break;
    case Core::LinearSolver::PreconditionerType::cheap_simple:
      beloslist.set("Preconditioner Type", "CheapSIMPLE");
      break;
    default:
      FOUR_C_THROW("Unknown preconditioner for Belos");
      break;
  }

  // set parameters for Ifpack if used
  if (azprectyp == Core::LinearSolver::PreconditionerType::ilu)
  {
    Teuchos::ParameterList& ifpacklist = outparams.sublist("IFPACK Parameters");
    ifpacklist = translate_four_c_to_ifpack(inparams);
  }

  // set parameters for CheapSIMPLE if used
  if (azprectyp == Core::LinearSolver::PreconditionerType::cheap_simple)
  {
    Teuchos::ParameterList& simplelist = outparams.sublist("CheapSIMPLE Parameters");
    simplelist.set("Prec Type", "CheapSIMPLE");  // not used
    Teuchos::ParameterList& predictList = simplelist.sublist("Inverse1");
    predictList = Core::LinAlg::Solver::translate_solver_parameters(
        get_solver_params(inparams.get<int>("SUB_SOLVER1")), get_solver_params, verbosity);
    Teuchos::ParameterList& schurList = simplelist.sublist("Inverse2");
    schurList = Core::LinAlg::Solver::translate_solver_parameters(
        get_solver_params(inparams.get<int>("SUB_SOLVER2")), get_solver_params, verbosity);
  }

  // set parameters for ML if used
  if (azprectyp == Core::LinearSolver::PreconditionerType::multigrid_ml ||
      azprectyp == Core::LinearSolver::PreconditionerType::multigrid_ml_fluid2)
  {
    Teuchos::ParameterList& mllist = outparams.sublist("ML Parameters");
    mllist = translate_four_c_to_ml(inparams, get_solver_params, verbosity);
  }
  if (azprectyp == Core::LinearSolver::PreconditionerType::multigrid_muelu)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu Parameters");
    muelulist = translate_four_c_to_muelu(inparams, &beloslist);
  }
  if (azprectyp == Core::LinearSolver::PreconditionerType::multigrid_muelu_contactsp)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu (Contact) Parameters");
    muelulist = translate_four_c_to_muelu(inparams, &beloslist);
  }
  if (azprectyp == Core::LinearSolver::PreconditionerType::multigrid_muelu_beamsolid)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu (BeamSolid) Parameters");
    muelulist = translate_four_c_to_muelu(inparams, &beloslist);
  }
  if (azprectyp == Core::LinearSolver::PreconditionerType::block_gauss_seidel_2x2)
  {
    Teuchos::ParameterList& bgslist = outparams.sublist("BGS Parameters");
    bgslist.set("numblocks", 2);

    // currently, the number of Gauss-Seidel iterations and the relaxation
    // parameter on the global level are set to 1 and 1.0, respectively
    bgslist.set("global_iter", 1);
    bgslist.set("global_omega", inparams.get<double>("BGS2X2_GLOBAL_DAMPING"));

    // the order of blocks in the given EpetraOperator can be changed in the
    // Gauss-Seidel procedure,
    // default: fliporder == 0, i.e., solve block1 --> block2
    std::string fliporder = inparams.get<std::string>("BGS2X2_FLIPORDER");
    bgslist.set("fliporder", (fliporder == "block1_block0_order") ? true : false);

    // currently, the number of Richardson iteratios and the relaxation
    // parameter on the individual block level are set to 1 and 1.0, respectively
    bgslist.set("block1_iter", 1);
    bgslist.set("block1_omega", inparams.get<double>("BGS2X2_BLOCK1_DAMPING"));
    bgslist.set("block2_iter", 1);
    bgslist.set("block2_omega", inparams.get<double>("BGS2X2_BLOCK2_DAMPING"));
  }
  if (azprectyp == Core::LinearSolver::PreconditionerType::multigrid_nxn)
  {
    Teuchos::ParameterList& amgnxnlist = outparams.sublist("AMGnxn Parameters");
    std::string amgnxn_xml = inparams.get<std::string>("AMGNXN_XML_FILE");
    amgnxnlist.set<std::string>("AMGNXN_XML_FILE", amgnxn_xml);
    std::string amgnxn_type = inparams.get<std::string>("AMGNXN_TYPE");
    amgnxnlist.set<std::string>("AMGNXN_TYPE", amgnxn_type);
  }

  return outparams;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::ParameterList Core::LinAlg::Solver::translate_solver_parameters(
    const Teuchos::ParameterList& inparams,
    const std::function<const Teuchos::ParameterList&(int)>& get_solver_params,
    Core::IO::Verbositylevel verbosity)
{
  TEUCHOS_FUNC_TIME_MONITOR("Core::LinAlg::Solver:  0)   translate_solver_parameters");

  Teuchos::ParameterList outparams;
  if (inparams.isParameter("NAME"))
    outparams.set<std::string>("name", inparams.get<std::string>("NAME"));

  switch (Teuchos::getIntegralValue<Core::LinearSolver::SolverType>(inparams, "SOLVER"))
  {
    case Core::LinearSolver::SolverType::undefined:
      std::cout << "undefined solver! Set " << inparams.name() << "  in your dat file!"
                << std::endl;
      FOUR_C_THROW("fix your dat file");
      break;
    case Core::LinearSolver::SolverType::umfpack:
      outparams.set("solver", "umfpack");
      break;
    case Core::LinearSolver::SolverType::superlu:
      outparams.set("solver", "superlu");
      break;
    case Core::LinearSolver::SolverType::belos:
      outparams = translate_four_c_to_belos(inparams, get_solver_params, verbosity);
      break;
    default:
      FOUR_C_THROW("Unsupported type of solver");
      break;
  }

  return outparams;
}

FOUR_C_NAMESPACE_CLOSE
