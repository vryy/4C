/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation

\level 0

*-----------------------------------------------------------------------*/
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <Epetra_MpiComm.h>
#include <Epetra_LinearProblem.h>

#include "baci_inpar_solver.H"
#include "baci_io_pstream.H"
#include "baci_lib_globalproblem.H"  // access global problem. can we avoid this?
#include "baci_linear_solver_method_linalg.H"
#include "baci_linalg_sparsematrix.H"

// some more Trilinos headers
#include <BelosTypes.hpp>                 // for Belos verbosity codes
#include <ml_MultiLevelPreconditioner.h>  // includes for ML parameter list validation

// solver interfaces
#include "baci_linear_solver_method_direct.H"
#include "baci_linear_solver_method_belos.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::Solver::Solver(
    Teuchos::RCP<Teuchos::ParameterList>& params, const Epetra_Comm& comm, FILE* outfile)
    : comm_(comm), params_(params), outfile_(outfile)
{
  Setup();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::Solver::Solver(const Epetra_Comm& comm, FILE* outfile)
    : comm_(comm), params_(Teuchos::rcp(new Teuchos::ParameterList())), outfile_(outfile)
{
  // set the default solver to umfpack
  Params().set("solver", "umfpack");
  Setup();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::Solver::Solver(
    const Teuchos::ParameterList& inparams, const Epetra_Comm& comm, FILE* outfile)
    : comm_(comm), params_(Teuchos::rcp(new Teuchos::ParameterList())), outfile_(outfile)
{
  *params_ = TranslateSolverParameters(inparams);
  Setup();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::Solver::Setup() { solver_ = Teuchos::null; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
CORE::LINALG::Solver::~Solver() { Reset(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::Solver::Reset() { solver_ = Teuchos::null; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::LINALG::Solver::getNumIters() const { return solver_->getNumIters(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::Solver::AdaptTolerance(
    const double desirednlnres, const double currentnlnres, const double better)
{
  if (!Params().isSublist("Belos Parameters")) return;

  const int myrank = Comm().MyPID();

  Teuchos::ParameterList& solverParams = Params().sublist("Belos Parameters");
  int output = solverParams.get<int>("Output Frequency", 1);

  std::string convtest = solverParams.get<std::string>(
      "Implicit Residual Scaling", Belos::convertScaleTypeToString(Belos::ScaleType::None));

  if (convtest != Belos::convertScaleTypeToString(Belos::ScaleType::NormOfInitRes))
    dserror("Using convergence adaptivity: Use AZ_r0 in input file");

  bool havesavedvalue = solverParams.isParameter("Convergence Tolerance Save");
  if (!havesavedvalue)
  {
    if (!solverParams.isParameter("Convergence Tolerance"))
    {
      std::cout << solverParams;
      dserror("No iterative solver tolerance in ParameterList");
    }
    solverParams.set<double>(
        "Convergence Tolerance Save", solverParams.get<double>("Convergence Tolerance", 1.e-8));
  }

  double tol = solverParams.get<double>("Convergence Tolerance Save", 1.e-8);

  if (!myrank && output)
    printf("                --- Solver input relative tolerance %10.3E\n", tol);
  if (currentnlnres * tol < desirednlnres)
  {
    double tolnew = desirednlnres * better / currentnlnres;
    if (tolnew > 1.0)
    {
      tolnew = 1.0;
      if (!myrank && output)
      {
        printf(
            "WARNING:  Computed adapted relative tolerance bigger than 1\n"
            "          Value constrained to 1, but consider adapting Parameter ADAPTCONV_BETTER\n");
      }
    }
    if (tolnew < tol) tolnew = tol;
    if (!myrank && output && tolnew > tol)
      printf("                *** Solver adapted relative tolerance %10.3E\n", tolnew);

    solverParams.set<double>("Convergence Tolerance", tolnew);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::Solver::ResetTolerance()
{
  if (!Params().isSublist("Belos Parameters")) return;

  Teuchos::ParameterList& solverParams = Params().sublist("Belos Parameters");

  bool havesavedvalue = solverParams.isParameter("Convergence Tolerance Save");
  if (!havesavedvalue) return;

  solverParams.set<double>(
      "Convergence Tolerance", solverParams.get<double>("Convergence Tolerance Save", 1.e-8));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void CORE::LINALG::Solver::Setup(Teuchos::RCP<Epetra_Operator> matrix,
    Teuchos::RCP<Epetra_MultiVector> x, Teuchos::RCP<Epetra_MultiVector> b, bool refactor,
    bool reset, Teuchos::RCP<CORE::LINALG::KrylovProjector> projector)
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::Solver:  1)   Setup");

  // reset data flags on demand
  if (reset)
  {
    Reset();
    refactor = true;
  }

  if (solver_ == Teuchos::null)
  {
    // decide what solver to use
    std::string solvertype = Params().get("solver", "none");

    if ("belos" == solvertype)
    {
      solver_ =
          Teuchos::rcp(new CORE::LINEAR_SOLVER::BelosSolver<Epetra_Operator, Epetra_MultiVector>(
              comm_, Params(), outfile_));
    }
    else if ("umfpack" == solvertype or "superlu" == solvertype)
    {
      solver_ = Teuchos::rcp(
          new CORE::LINEAR_SOLVER::DirectSolver<Epetra_Operator, Epetra_MultiVector>(solvertype));
    }
    else
      dserror("Unknown type of solver");
  }

  solver_->Setup(matrix, x, b, refactor, reset, projector);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::LINALG::Solver::Solve()
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::Solver:  2)   Solve");
  return solver_->Solve();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::LINALG::Solver::Solve(Teuchos::RCP<Epetra_Operator> matrix,
    Teuchos::RCP<Epetra_MultiVector> x, Teuchos::RCP<Epetra_MultiVector> b, bool refactor,
    bool reset, Teuchos::RCP<CORE::LINALG::KrylovProjector> projector)
{
  Setup(matrix, x, b, refactor, reset, projector);
  return Solve();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int CORE::LINALG::Solver::NoxSolve(Epetra_LinearProblem& linProblem, bool refactor, bool reset,
    Teuchos::RCP<CORE::LINALG::KrylovProjector> projector)
{
  Teuchos::RCP<Epetra_Operator> matrix = Teuchos::rcp(linProblem.GetOperator(), false);
  Teuchos::RCP<Epetra_MultiVector> x = Teuchos::rcp(linProblem.GetLHS(), false);
  Teuchos::RCP<Epetra_MultiVector> b = Teuchos::rcp(linProblem.GetRHS(), false);

  return this->Solve(matrix, x, b, refactor, reset, projector);
}


/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList CORE::LINALG::Solver::TranslateBACIToIfpack(
    const Teuchos::ParameterList& inparams)
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
const Teuchos::ParameterList CORE::LINALG::Solver::TranslateBACIToML(
    const Teuchos::ParameterList& inparams, Teuchos::ParameterList* azlist)
{
  Teuchos::ParameterList mllist;

  ML_Epetra::SetDefaults("SA", mllist);
  const auto prectyp =
      Teuchos::getIntegralValue<INPAR::SOLVER::PreconditionerType>(inparams, "AZPREC");

  switch (prectyp)
  {
    case INPAR::SOLVER::PreconditionerType::multigrid_ml:
      break;
    case INPAR::SOLVER::PreconditionerType::multigrid_ml_fluid:
      mllist.set("aggregation: use tentative restriction", true);
      break;
    case INPAR::SOLVER::PreconditionerType::multigrid_ml_fluid2:
      mllist.set("energy minimization: enable", true);
      mllist.set("energy minimization: type", 3);
      mllist.set("aggregation: block scaling", false);
      break;
    default:
      dserror("Unknown type of ml preconditioner");
      break;
  }

  // set repartitioning parameters
  // En-/Disable ML repartitioning. Note: ML requires parameter to be set as integer.
  bool doRepart = ::DRT::INPUT::IntegralValue<bool>(inparams, "ML_REBALANCE");
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
  switch (::DRT::INPUT::IntegralValue<int>(inparams, "ML_COARSEN"))
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
      dserror("Unknown type of coarsening for ML");
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
    dserror("Not enough smoothing steps ML_SMOTIMES=%d, must be larger/equal than ML_MAXLEVEL=%d\n",
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
      type = ::DRT::INPUT::IntegralValue<int>(inparams, "ML_SMOOTHERFINE");
      damp = inparams.get<double>("ML_DAMPFINE");
    }
    else if (i < mlmaxlevel - 1)
    {
      type = ::DRT::INPUT::IntegralValue<int>(inparams, "ML_SMOOTHERMED");
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
        SchurCompList = TranslateSolverParameters(
            ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
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
        predictList = TranslateSolverParameters(
            ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
        Teuchos::ParameterList& SchurCompList = smolevelsublist.sublist("smoother: SchurComp list");
        SchurCompList = TranslateSolverParameters(
            ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
      }
      break;
      case 13:  // IBD: indefinite block diagonal preconditioner
      {
        smolevelsublist.set("smoother: type", "IBD");
        smolevelsublist.set("smoother: sweeps", mlsmotimessteps[i]);
        smolevelsublist.set("smoother: damping factor", damp);
        Teuchos::ParameterList& predictList = smolevelsublist.sublist("smoother: Predictor list");
        predictList = TranslateSolverParameters(
            ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
        Teuchos::ParameterList& SchurCompList = smolevelsublist.sublist("smoother: SchurComp list");
        SchurCompList = TranslateSolverParameters(
            ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
      }
      break;
      case 14:  // Uzawa: inexact Uzawa smoother
      {
        smolevelsublist.set("smoother: type", "Uzawa");
        smolevelsublist.set("smoother: sweeps", mlsmotimessteps[i]);
        smolevelsublist.set("smoother: damping factor", damp);
        Teuchos::ParameterList& predictList = smolevelsublist.sublist("smoother: Predictor list");
        predictList = TranslateSolverParameters(
            ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
        Teuchos::ParameterList& SchurCompList = smolevelsublist.sublist("smoother: SchurComp list");
        SchurCompList = TranslateSolverParameters(
            ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
      }
      break;
      default:
        dserror("Unknown type of smoother for ML: tuple %d", type);
        break;
    }  // switch (type)
  }    // for (int i=0; i<azvar->mlmaxlevel-1; ++i)

  // set coarse grid solver
  const int coarse = mlmaxlevel - 1;
  switch (::DRT::INPUT::IntegralValue<int>(inparams, "ML_SMOOTHERCOARSE"))
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
      SchurCompList = TranslateSolverParameters(
          ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
    }
    break;
    case 11:  // SIMPLE smoother  (only for MueLu with BlockedOperators)
    case 12:  // SIMPLEC smoother (only for MueLu with BlockedOperators)
    {
      int type = ::DRT::INPUT::IntegralValue<int>(inparams, "ML_SMOOTHERCOARSE");
      if (type == 11)
        mllist.set("coarse: type", "SIMPLE");
      else if (type == 12)
        mllist.set("coarse: type", "SIMPLEC");
      mllist.set("coarse: sweeps", mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor", inparams.get<double>("ML_DAMPCOARSE"));
      Teuchos::ParameterList& predictList = mllist.sublist("coarse: Predictor list");
      predictList = TranslateSolverParameters(
          ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
      Teuchos::ParameterList& SchurCompList = mllist.sublist("coarse: SchurComp list");
      SchurCompList = TranslateSolverParameters(
          ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
    }
    break;
    case 13:  // IBD: indefinite block diagonal preconditioner
    {
      mllist.set("coarse: type", "IBD");
      mllist.set("coarse: sweeps", mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor", inparams.get<double>("ML_DAMPCOARSE"));
      Teuchos::ParameterList& predictList = mllist.sublist("coarse: Predictor list");
      predictList = TranslateSolverParameters(
          ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
      Teuchos::ParameterList& SchurCompList = mllist.sublist("coarse: SchurComp list");
      SchurCompList = TranslateSolverParameters(
          ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
    }
    break;
    case 14:  // Uzawa: inexact Uzawa smoother
    {
      mllist.set("coarse: type", "Uzawa");
      mllist.set("coarse: sweeps", mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor", inparams.get<double>("ML_DAMPCOARSE"));
      Teuchos::ParameterList& predictList = mllist.sublist("coarse: Predictor list");
      predictList = TranslateSolverParameters(
          ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
      Teuchos::ParameterList& SchurCompList = mllist.sublist("coarse: SchurComp list");
      SchurCompList = TranslateSolverParameters(
          ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
    }
    break;
    default:
      dserror("Unknown type of coarse solver for ML");
      break;
  }  // switch (azvar->mlsmotype_coarse)
  // default values for nullspace
  mllist.set("PDE equations", 1);
  mllist.set("null space: dimension", 1);
  mllist.set("null space: type", "pre-computed");
  mllist.set("null space: add default vectors", false);
  mllist.set<double*>("null space: vectors", NULL);

  return mllist;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList CORE::LINALG::Solver::TranslateBACIToMuelu(
    const Teuchos::ParameterList& inparams, Teuchos::ParameterList* azlist)
{
  Teuchos::ParameterList muelulist;

  std::string xmlfile = inparams.get<std::string>("MUELU_XML_FILE");
  if (xmlfile != "none") muelulist.set("MUELU_XML_FILE", xmlfile);

  muelulist.set<bool>(
      "MUELU_XML_ENFORCE", ::DRT::INPUT::IntegralValue<bool>(inparams, "MUELU_XML_ENFORCE"));
  muelulist.set<bool>("CORE::LINALG::MueLu_Preconditioner", true);

  return muelulist;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList CORE::LINALG::Solver::TranslateBACIToBelos(
    const Teuchos::ParameterList& inparams)
{
  Teuchos::ParameterList outparams;
  outparams.set("solver", "belos");
  Teuchos::ParameterList& beloslist = outparams.sublist("Belos Parameters");

  // set verbosity
  auto verbosityLevel = ::DRT::INPUT::IntegralValue<IO::verbositylevel>(
      ::DRT::Problem::Instance()->IOParams(), "VERBOSITY");

  switch (verbosityLevel)
  {
    case IO::minimal:
      beloslist.set("Output Style", Belos::OutputType::Brief);
      beloslist.set("Verbosity", Belos::MsgType::Warnings);
      break;
    case IO::standard:
      beloslist.set("Output Style", Belos::OutputType::Brief);
      beloslist.set("Verbosity", Belos::MsgType::Warnings + Belos::MsgType::StatusTestDetails);
      break;
    case IO::verbose:
      beloslist.set("Output Style", Belos::OutputType::General);
      beloslist.set("Verbosity", Belos::MsgType::Warnings + Belos::MsgType::StatusTestDetails +
                                     Belos::MsgType::FinalSummary);
      break;
    case IO::debug:
      beloslist.set("Output Style", Belos::OutputType::General);
      beloslist.set("Verbosity", Belos::MsgType::Debug);
    default:
      break;
  }
  beloslist.set("Output Frequency", inparams.get<int>("AZOUTPUT"));

  // set tolerances and iterations
  beloslist.set("Maximum Iterations", inparams.get<int>("AZITER"));
  beloslist.set("Convergence Tolerance", inparams.get<double>("AZTOL"));
  beloslist.set("reuse", inparams.get<int>("AZREUSE"));
  beloslist.set("ncall", 0);
  beloslist.set("Implicit Residual Scaling",
      Belos::convertScaleTypeToString(
          Teuchos::getIntegralValue<Belos::ScaleType>(inparams, "AZCONV")));

  // set type of solver
  switch (Teuchos::getIntegralValue<INPAR::SOLVER::IterativeSolverType>(inparams, "AZSOLVE"))
  {
    case INPAR::SOLVER::IterativeSolverType::cg:
      beloslist.set("Solver Type", "CG");
      break;
    case INPAR::SOLVER::IterativeSolverType::bicgstab:
      beloslist.set("Solver Type", "BiCGSTAB");
      break;
    case INPAR::SOLVER::IterativeSolverType::gmres:
      beloslist.set("Solver Type", "GMRES");
      beloslist.set("Num Blocks", inparams.get<int>("AZSUB"));
      break;
    default:
    {
      dserror("Flag '%s'! \nUnknown solver for Belos.",
          Teuchos::getIntegralValue<INPAR::SOLVER::IterativeSolverType>(inparams, "AZSOLVE"));
      break;
    }
  }

  // set type of preconditioner
  const auto azprectyp =
      Teuchos::getIntegralValue<INPAR::SOLVER::PreconditionerType>(inparams, "AZPREC");

  switch (azprectyp)
  {
    case INPAR::SOLVER::PreconditionerType::ilu:
      beloslist.set("Preconditioner Type", "ILU");
      break;
    case INPAR::SOLVER::PreconditionerType::icc:
      beloslist.set("Preconditioner Type", "IC");
      break;
    case INPAR::SOLVER::PreconditionerType::multigrid_ml:
    case INPAR::SOLVER::PreconditionerType::multigrid_ml_fluid:
    case INPAR::SOLVER::PreconditionerType::multigrid_ml_fluid2:
    case INPAR::SOLVER::PreconditionerType::multigrid_muelu:
      beloslist.set("Preconditioner Type", "ML");
      break;
    case INPAR::SOLVER::PreconditionerType::multigrid_muelu_fluid:
      beloslist.set("Preconditioner Type", "Fluid");
      break;
    case INPAR::SOLVER::PreconditionerType::multigrid_muelu_tsi:
      beloslist.set("Preconditioner Type", "TSI");
      break;
    case INPAR::SOLVER::PreconditionerType::multigrid_muelu_contactsp:
      beloslist.set("Preconditioner Type", "ContactSP");
      break;
    case INPAR::SOLVER::PreconditionerType::multigrid_muelu_beamsolid:
      beloslist.set("Preconditioner Type", "BeamSolid");
      break;
    case INPAR::SOLVER::PreconditionerType::multigrid_muelu_fsi:
      beloslist.set("Preconditioner Type", "FSI");
      break;
    case INPAR::SOLVER::PreconditionerType::multigrid_nxn:
      beloslist.set("Preconditioner Type", "AMGnxn");
      break;
    case INPAR::SOLVER::PreconditionerType::block_gauss_seidel_2x2:
      beloslist.set("Preconditioner Type", "ML");
      break;
    case INPAR::SOLVER::PreconditionerType::cheap_simple:
      beloslist.set("Preconditioner Type", "CheapSIMPLE");
      break;
    default:
      dserror("Unknown preconditioner for Belos");
      break;
  }

  // set scaling of linear problem
  switch (Teuchos::getIntegralValue<INPAR::SOLVER::ScalingStrategy>(inparams, "AZSCAL"))
  {
    case INPAR::SOLVER::ScalingStrategy::none:
      beloslist.set("scaling", "none");
      break;
    case INPAR::SOLVER::ScalingStrategy::symmetric:
      beloslist.set("scaling", "symmetric");
      break;
    case INPAR::SOLVER::ScalingStrategy::infnorm:
      beloslist.set("scaling", "infnorm");
      break;
    default:
      dserror("No valid scaling method selected. Choose between \"none\", \"sym\" or \"infnorm\".");
  }

  // set permutation of linear problem
  switch (Teuchos::getIntegralValue<INPAR::SOLVER::PermutationStrategy>(inparams, "PERMUTE_SYSTEM"))
  {
    case INPAR::SOLVER::PermutationStrategy::algebraic:
      beloslist.set("permutation strategy", "Algebraic");
      break;
    case INPAR::SOLVER::PermutationStrategy::local:
      beloslist.set("permutation strategy", "Local");
      break;
    case INPAR::SOLVER::PermutationStrategy::none:
    default:
      beloslist.set("permutation strategy", "none");
      break;
  }
  double nonDiagDominance = inparams.get<double>("NON_DIAGDOMINANCE_RATIO");
  beloslist.set("diagonal dominance ratio", nonDiagDominance);

  // set parameters for Ifpack if used
  if (azprectyp == INPAR::SOLVER::PreconditionerType::ilu ||
      azprectyp == INPAR::SOLVER::PreconditionerType::icc)
  {
    Teuchos::ParameterList& ifpacklist = outparams.sublist("IFPACK Parameters");
    ifpacklist = CORE::LINALG::Solver::TranslateBACIToIfpack(inparams);
  }

  // set parameters for CheapSIMPLE if used
  if (azprectyp == INPAR::SOLVER::PreconditionerType::cheap_simple)
  {
    Teuchos::ParameterList& simplelist = outparams.sublist("CheapSIMPLE Parameters");
    simplelist.set("Prec Type", "CheapSIMPLE");  // not used
    Teuchos::ParameterList& predictList = simplelist.sublist("Inverse1");
    predictList = TranslateSolverParameters(
        ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
    std::cout << predictList << std::endl;
    Teuchos::ParameterList& schurList = simplelist.sublist("Inverse2");
    schurList = TranslateSolverParameters(
        ::DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
  }

  // set parameters for ML if used
  if (azprectyp == INPAR::SOLVER::PreconditionerType::multigrid_ml ||
      azprectyp == INPAR::SOLVER::PreconditionerType::multigrid_ml_fluid ||
      azprectyp == INPAR::SOLVER::PreconditionerType::multigrid_ml_fluid2)
  {
    Teuchos::ParameterList& mllist = outparams.sublist("ML Parameters");
    mllist = CORE::LINALG::Solver::TranslateBACIToML(inparams, &beloslist);
  }
  if (azprectyp == INPAR::SOLVER::PreconditionerType::multigrid_muelu)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu Parameters");
    muelulist = CORE::LINALG::Solver::TranslateBACIToMuelu(inparams, &beloslist);
  }
  if (azprectyp == INPAR::SOLVER::PreconditionerType::multigrid_muelu_fluid)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu (Fluid) Parameters");
    muelulist = CORE::LINALG::Solver::TranslateBACIToMuelu(inparams, &beloslist);
  }
  if (azprectyp == INPAR::SOLVER::PreconditionerType::multigrid_muelu_tsi)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu (TSI) Parameters");
    muelulist = CORE::LINALG::Solver::TranslateBACIToMuelu(inparams, &beloslist);
  }
  if (azprectyp == INPAR::SOLVER::PreconditionerType::multigrid_muelu_contactsp)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu (Contact) Parameters");
    muelulist = CORE::LINALG::Solver::TranslateBACIToMuelu(inparams, &beloslist);
  }
  if (azprectyp == INPAR::SOLVER::PreconditionerType::multigrid_muelu_beamsolid)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu (BeamSolid) Parameters");
    muelulist = CORE::LINALG::Solver::TranslateBACIToMuelu(inparams, &beloslist);
  }
  if (azprectyp == INPAR::SOLVER::PreconditionerType::multigrid_muelu_fsi)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu (FSI) Parameters");
    muelulist = CORE::LINALG::Solver::TranslateBACIToMuelu(inparams, &beloslist);
  }
  if (azprectyp == INPAR::SOLVER::PreconditionerType::block_gauss_seidel_2x2)
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
  if (azprectyp == INPAR::SOLVER::PreconditionerType::multigrid_nxn)
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
const Teuchos::ParameterList CORE::LINALG::Solver::TranslateSolverParameters(
    const Teuchos::ParameterList& inparams)
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::Solver:  0)   TranslateSolverParameters");

  Teuchos::ParameterList outparams;
  outparams.set<std::string>("name", inparams.get<std::string>("NAME"));

  switch (Teuchos::getIntegralValue<INPAR::SOLVER::SolverType>(inparams, "SOLVER"))
  {
    case INPAR::SOLVER::SolverType::undefined:
      std::cout << "undefined solver! Set " << inparams.name() << "  in your dat file!"
                << std::endl;
      dserror("fix your dat file");
      break;
    case INPAR::SOLVER::SolverType::umfpack:
      outparams.set("solver", "umfpack");
      break;
    case INPAR::SOLVER::SolverType::superlu:
      outparams.set("solver", "superlu");
      break;
    case INPAR::SOLVER::SolverType::belos:
      outparams = CORE::LINALG::Solver::TranslateBACIToBelos(inparams);
      break;
    default:
      dserror("Unsupported type of solver");
      break;
  }

  return outparams;
}
