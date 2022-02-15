/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation

\level 0

*-----------------------------------------------------------------------*/

#undef WRITEOUTSTATISTICS
#ifdef WRITEOUTSTATISTICS
#include "Teuchos_Time.hpp"
#endif

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <Epetra_MpiComm.h>
#include <Epetra_LinearProblem.h>

#include "../drt_inpar/inpar_solver.H"
#include "../drt_lib/drt_globalproblem.H"  // access global problem. can we avoid this?
#include "linalg_solver.H"
#include "linalg_sparsematrix.H"

// some more Trilinos headers
#include <BelosTypes.hpp>                 // for Belos verbosity codes
#include <ml_MultiLevelPreconditioner.h>  // includes for ML parameter list validation
#include <az_aztec_defs.h>                // for translation of parameters

// solver interfaces
#include "../solver/solver_directsolver.H"  // Amesos
#include "../solver/solver_aztecsolver.H"   // AztecOO
#include "../solver/solver_belossolver.H"   // Belos

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::Solver::Solver(
    Teuchos::RCP<Teuchos::ParameterList> params, const Epetra_Comm& comm, FILE* outfile)
    : comm_(comm), params_(params), outfile_(outfile)
{
  Setup();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::Solver::Solver(const Epetra_Comm& comm, FILE* outfile)
    : comm_(comm), params_(Teuchos::rcp(new Teuchos::ParameterList())), outfile_(outfile)
{
  Params().set("solver", "klu");
  Params().set("symmetric", false);
  Setup();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::Solver::Solver(
    const Teuchos::ParameterList& inparams, const Epetra_Comm& comm, FILE* outfile)
    : comm_(comm), params_(Teuchos::rcp(new Teuchos::ParameterList())), outfile_(outfile)
{
  *params_ = TranslateSolverParameters(inparams);
  Setup();
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Solver::Setup() { solver_ = Teuchos::null; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
LINALG::Solver::~Solver() { Reset(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Solver::Reset() { solver_ = Teuchos::null; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
std::ostream& operator<<(std::ostream& os, const LINALG::Solver& solver)
{
  solver.Print(os);
  return os;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Solver::Print(std::ostream& os) const
{
  if (Comm().MyPID() == 0)
  {
    os << "============================LINALG::Solver Parameter List\n";
    os << *params_;
    os << "========================end LINALG::Solver Parameter List\n";
  }
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::Solver::getNumIters() const { return solver_->getNumIters(); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Solver::AdaptTolerance(
    const double desirednlnres, const double currentnlnres, const double better)
{
  if (!Params().isSublist("Aztec Parameters")) return;
  const int myrank = Comm().MyPID();
  Teuchos::ParameterList& azlist = Params().sublist("Aztec Parameters");
  int output = azlist.get<int>("AZ_output", 1);
  int convtest = azlist.get<int>("AZ_conv", AZ_noscaled);
  if (convtest != AZ_r0) dserror("Using convergence adaptivity: Use AZ_r0 in input file");
  bool havesavedvalue = azlist.isParameter("AZ_tol save");
  if (!havesavedvalue)
  {
    if (!azlist.isParameter("AZ_tol"))
    {
      std::cout << azlist;
      dserror("No Aztec tolerance in ParameterList");
    }
    azlist.set<double>("AZ_tol save", azlist.get<double>("AZ_tol", 1.e-8));
  }
  double tol = azlist.get<double>("AZ_tol save", 1.e-8);
  if (!myrank && output)
    printf("                --- Aztec input   relative tolerance %10.3E\n", tol);
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
      printf("                *** Aztec adapted relative tolerance %10.3E\n", tolnew);
    azlist.set<double>("AZ_tol", tolnew);
  }

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Solver::ResetTolerance()
{
  if (!Params().isSublist("Aztec Parameters")) return;
  Teuchos::ParameterList& azlist = Params().sublist("Aztec Parameters");
  bool havesavedvalue = azlist.isParameter("AZ_tol save");
  if (!havesavedvalue) return;
  azlist.set<double>("AZ_tol", azlist.get<double>("AZ_tol save", 1.e-8));
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Solver::Setup(Teuchos::RCP<Epetra_Operator> matrix, Teuchos::RCP<Epetra_MultiVector> x,
    Teuchos::RCP<Epetra_MultiVector> b, bool refactor, bool reset,
    Teuchos::RCP<LINALG::KrylovProjector> projector)
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::Solver:  1)   Setup");

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

    if ("aztec" == solvertype)
    {
      solver_ = Teuchos::rcp(new LINALG::SOLVER::AztecSolver(comm_, Params(), outfile_));
    }
    else if ("belos" == solvertype)
    {
      solver_ = Teuchos::rcp(new LINALG::SOLVER::BelosSolver(comm_, Params(), outfile_));
    }
    else if ("klu" == solvertype or "umfpack" == solvertype or "superlu" == solvertype or
             "lapack" == solvertype)
    {
      solver_ = Teuchos::rcp(new LINALG::SOLVER::DirectSolver(solvertype));
    }
    else
      dserror("Unknown type of solver");
  }

  solver_->Setup(matrix, x, b, refactor, reset, projector);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::Solver::Solve()
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::Solver:  2)   Solve");
  return solver_->Solve();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::Solver::ApplyInverse(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  return solver_->ApplyInverse(X, Y);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::Solver::Solve(Teuchos::RCP<Epetra_Operator> matrix, Teuchos::RCP<Epetra_MultiVector> x,
    Teuchos::RCP<Epetra_MultiVector> b, bool refactor, bool reset,
    Teuchos::RCP<LINALG::KrylovProjector> projector)
{
  Setup(matrix, x, b, refactor, reset, projector);
  return Solve();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int LINALG::Solver::NoxSolve(Epetra_LinearProblem& linProblem, bool refactor, bool reset,
    Teuchos::RCP<LINALG::KrylovProjector> projector)
{
  Teuchos::RCP<Epetra_Operator> matrix = Teuchos::rcp(linProblem.GetOperator(), false);
  Teuchos::RCP<Epetra_MultiVector> x = Teuchos::rcp(linProblem.GetLHS(), false);
  Teuchos::RCP<Epetra_MultiVector> b = Teuchos::rcp(linProblem.GetRHS(), false);

  return this->Solve(matrix, x, b, refactor, reset, projector);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void LINALG::Solver::FixMLNullspace(std::string field, const Epetra_Map& oldmap,
    const Epetra_Map& newmap, Teuchos::ParameterList& solveparams)
{
  // there is no ML list, do nothing
  if (!solveparams.isSublist("ML Parameters") && !solveparams.isSublist("MueLu Parameters")) return;

  // find the ML or MueLu list
  Teuchos::RCP<Teuchos::ParameterList> params_ptr = Teuchos::null;
  if (solveparams.isSublist("ML Parameters"))
    params_ptr = Teuchos::rcp(&(solveparams.sublist("ML Parameters")), false);
  else if (solveparams.isSublist("MueLu Parameters"))
    params_ptr = Teuchos::rcp(&(solveparams.sublist("MueLu Parameters")), false);
  else
    return;
  Teuchos::ParameterList& params = *params_ptr;

  const int ndim = params.get("null space: dimension", -1);
  if (ndim == -1) dserror("List does not contain nullspace dimension");

  Teuchos::RCP<std::vector<double>> ns =
      params.get<Teuchos::RCP<std::vector<double>>>("nullspace", Teuchos::null);
  if (ns == Teuchos::null) dserror("List does not contain nullspace");
  double* ons = &((*ns)[0]);

  const int olength = (int)ns->size() / ndim;

  const int nlength = newmap.NumMyElements();

  if (olength == nlength) return;  // everything should be ok, do nothing

  if (olength != oldmap.NumMyElements())
    dserror("Nullspace does not match old map length, %d %d", olength, oldmap.NumMyElements());

  if (nlength > olength)
    dserror("New problem size larger than old - full rebuild of nullspace neccessary");

  // Allocate a new nullspace and fill it
  Teuchos::RCP<std::vector<double>> nsnew =
      Teuchos::rcp(new std::vector<double>(nlength * ndim, 0.0));
  double* nns = &((*nsnew)[0]);

  for (int i = 0; i < nlength; ++i)
  {
    int gid = newmap.GID(i);
    int olid = oldmap.LID(gid);
    if (olid == -1) continue;

    // transfer entries for this dof to new nullspace vector
    for (int j = 0; j < ndim; ++j) nns[j * nlength + i] = ons[j * olength + olid];
  }

  // put new nullspace in parameter list
  // this print message can go away at some point
  if (!oldmap.Comm().MyPID()) printf("Fixing %s ML Nullspace\n", field.c_str());
  params.set<Teuchos::RCP<std::vector<double>>>("nullspace", nsnew);
  params.set<double*>("null space: vectors", nns);

  return;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList LINALG::Solver::TranslateBACIToIfpack(
    const Teuchos::ParameterList& inparams)
{
  Teuchos::ParameterList ifpacklist;

  ifpacklist.set("relaxation: damping factor", inparams.get<double>("AZOMEGA"));
  ifpacklist.set("fact: drop tolerance", inparams.get<double>("AZDROP"));
  ifpacklist.set("fact: level-of-fill", inparams.get<int>("IFPACKGFILL"));
  ifpacklist.set("fact: ilut level-of-fill", inparams.get<double>("IFPACKFILL"));
  ifpacklist.set("partitioner: overlap", inparams.get<int>("IFPACKOVERLAP"));
  ifpacklist.set("schwarz: combine mode",
      inparams.get<std::string>("IFPACKCOMBINE"));    // can be "Zero", "Add", "Insert"
  ifpacklist.set("schwarz: reordering type", "rcm");  // "rcm" or "metis" or "amd"
  ifpacklist.set("amesos: solver type",
      "Amesos_Klu");  // can be "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu"

  //------------------------------------- set type of preconditioner
  const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(inparams, "AZPREC");
  switch (prectyp)
  {
    case INPAR::SOLVER::azprec_Jacobi:
    {
      ifpacklist.set("relaxation: type", "Jacobi");
      ifpacklist.set("relaxation: sweeps", inparams.get<int>("IFPACKGFILL"));
    }
    break;
    case INPAR::SOLVER::azprec_SymmGaussSeidel:
    {
      ifpacklist.set("relaxation: type", "symmetric Gauss-Seidel");
      ifpacklist.set(
          "relaxation: sweeps", inparams.get<int>("IFPACKGFILL"));  // misuse IFPACKGFILL parameter
    }
    break;
    case INPAR::SOLVER::azprec_GaussSeidel:
    {
      ifpacklist.set("relaxation: type", "Gauss-Seidel");
      ifpacklist.set(
          "relaxation: sweeps", inparams.get<int>("IFPACKGFILL"));  // misuse IFPACKGFILL parameter
    }
    break;
    case INPAR::SOLVER::azprec_DownwindGaussSeidel:
    {
      // in case of downwinding prevent ifpack from again reordering
      ifpacklist.set("schwarz: reordering type", "none");
      ifpacklist.set("relaxation: type", "Gauss-Seidel");
      ifpacklist.set(
          "relaxation: sweeps", inparams.get<int>("IFPACKGFILL"));  // misuse IFPACKGFILL parameter
    }
    break;
    case INPAR::SOLVER::azprec_Chebyshev:
    {
      ifpacklist.set("chebyshev: degree", inparams.get<int>("IFPACKGFILL"));
    }
    break;
  }

  return ifpacklist;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList LINALG::Solver::TranslateBACIToML(
    const Teuchos::ParameterList& inparams, Teuchos::ParameterList* azlist)
{
  Teuchos::ParameterList mllist;

  ML_Epetra::SetDefaults("SA", mllist);
  const int prectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(inparams, "AZPREC");

  switch (prectyp)
  {
    case INPAR::SOLVER::azprec_ML:  // do nothing, this is standard
      break;
    case INPAR::SOLVER::azprec_MLfluid:  // unsymmetric, unsmoothed restriction
      mllist.set("aggregation: use tentative restriction", true);
      break;
    case INPAR::SOLVER::azprec_MLfluid2:  // full Pretrov-Galerkin unsymmetric smoothed
      mllist.set("energy minimization: enable", true);
      mllist.set("energy minimization: type", 3);  // 1,2,3 cheap -> expensive
      mllist.set("aggregation: block scaling", false);
      break;
    default:
      dserror("Unknown type of ml preconditioner");
      break;
  }

  // set repartitioning parameters
  // En-/Disable ML repartitioning. Note: ML requires parameter to be set as integer.
  bool doRepart = DRT::INPUT::IntegralValue<bool>(inparams, "ML_REBALANCE");
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
  switch (DRT::INPUT::IntegralValue<int>(inparams, "ML_COARSEN"))
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
      type = DRT::INPUT::IntegralValue<int>(inparams, "ML_SMOOTHERFINE");
      damp = inparams.get<double>("ML_DAMPFINE");
    }
    else if (i < mlmaxlevel - 1)
    {
      type = DRT::INPUT::IntegralValue<int>(inparams, "ML_SMOOTHERMED");
      damp = inparams.get<double>("ML_DAMPMED");
    }
    /*else
    {
      type = DRT::INPUT::IntegralValue<int>(inparams,"ML_SMOOTHERCOARSE");
      damp = inparams.get<double>("ML_DAMPCOARSE");
    }*/
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
      case 8:  // DGS
        smolevelsublist.set("smoother: type", "Gauss-Seidel");
        smolevelsublist.set("smoother: sweeps", mlsmotimessteps[i]);
        smolevelsublist.set("smoother: damping factor", damp);
        if (azlist != NULL)
        {
          azlist->set<bool>("downwinding", true);
          azlist->set<double>("downwinding tau", inparams.get<double>("DWINDTAU"));
        }
        else
        {
          std::cout << "WARNING: cannot set parameters for Downwinding Gauss Seidel" << std::endl;
        }
        {
          Teuchos::ParameterList& ifpacklist = mllist.sublist("smoother: ifpack list");
          ifpacklist.set("schwarz: reordering type", "true");
        }
        break;
      case 1:  // Jacobi
        smolevelsublist.set("smoother: type", "Jacobi");
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
            DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
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
            DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
        Teuchos::ParameterList& SchurCompList = smolevelsublist.sublist("smoother: SchurComp list");
        SchurCompList = TranslateSolverParameters(
            DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
      }
      break;
      case 13:  // IBD: indefinite block diagonal preconditioner
      {
        smolevelsublist.set("smoother: type", "IBD");
        smolevelsublist.set("smoother: sweeps", mlsmotimessteps[i]);
        smolevelsublist.set("smoother: damping factor", damp);
        Teuchos::ParameterList& predictList = smolevelsublist.sublist("smoother: Predictor list");
        predictList = TranslateSolverParameters(
            DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
        Teuchos::ParameterList& SchurCompList = smolevelsublist.sublist("smoother: SchurComp list");
        SchurCompList = TranslateSolverParameters(
            DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
      }
      break;
      case 14:  // Uzawa: inexact Uzawa smoother
      {
        smolevelsublist.set("smoother: type", "Uzawa");
        smolevelsublist.set("smoother: sweeps", mlsmotimessteps[i]);
        smolevelsublist.set("smoother: damping factor", damp);
        Teuchos::ParameterList& predictList = smolevelsublist.sublist("smoother: Predictor list");
        predictList = TranslateSolverParameters(
            DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
        Teuchos::ParameterList& SchurCompList = smolevelsublist.sublist("smoother: SchurComp list");
        SchurCompList = TranslateSolverParameters(
            DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
      }
      break;
      default:
        dserror("Unknown type of smoother for ML: tuple %d", type);
        break;
    }  // switch (type)
  }    // for (int i=0; i<azvar->mlmaxlevel-1; ++i)

  // set coarse grid solver
  const int coarse = mlmaxlevel - 1;
  switch (DRT::INPUT::IntegralValue<int>(inparams, "ML_SMOOTHERCOARSE"))
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
    case 8:
      mllist.set("coarse: type", "Gauss-Seidel");
      mllist.set("coarse: sweeps", mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor", inparams.get<double>("ML_DAMPCOARSE"));
      if (azlist != NULL)
      {
        azlist->set<bool>("downwinding", true);
        azlist->set<double>("downwinding tau", inparams.get<double>("DWINDTAU"));
      }
      else
      {
        std::cout << "WARNING: cannot set parameters for Downwinding Gauss Seidel" << std::endl;
      }
      {
        Teuchos::ParameterList& ifpacklist = mllist.sublist("smoother: ifpack list");
        ifpacklist.set("schwarz: reordering type", "true");
      }
      break;
    case 1:
      mllist.set("coarse: type", "Jacobi");
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
          DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
    }
    break;
    case 11:  // SIMPLE smoother  (only for MueLu with BlockedOperators)
    case 12:  // SIMPLEC smoother (only for MueLu with BlockedOperators)
    {
      int type = DRT::INPUT::IntegralValue<int>(inparams, "ML_SMOOTHERCOARSE");
      if (type == 11)
        mllist.set("coarse: type", "SIMPLE");
      else if (type == 12)
        mllist.set("coarse: type", "SIMPLEC");
      mllist.set("coarse: sweeps", mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor", inparams.get<double>("ML_DAMPCOARSE"));
      Teuchos::ParameterList& predictList = mllist.sublist("coarse: Predictor list");
      predictList = TranslateSolverParameters(
          DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
      Teuchos::ParameterList& SchurCompList = mllist.sublist("coarse: SchurComp list");
      SchurCompList = TranslateSolverParameters(
          DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
    }
    break;
    case 13:  // IBD: indefinite block diagonal preconditioner
    {
      mllist.set("coarse: type", "IBD");
      mllist.set("coarse: sweeps", mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor", inparams.get<double>("ML_DAMPCOARSE"));
      Teuchos::ParameterList& predictList = mllist.sublist("coarse: Predictor list");
      predictList = TranslateSolverParameters(
          DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
      Teuchos::ParameterList& SchurCompList = mllist.sublist("coarse: SchurComp list");
      SchurCompList = TranslateSolverParameters(
          DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
    }
    break;
    case 14:  // Uzawa: inexact Uzawa smoother
    {
      mllist.set("coarse: type", "Uzawa");
      mllist.set("coarse: sweeps", mlsmotimessteps[coarse]);
      mllist.set("coarse: damping factor", inparams.get<double>("ML_DAMPCOARSE"));
      Teuchos::ParameterList& predictList = mllist.sublist("coarse: Predictor list");
      predictList = TranslateSolverParameters(
          DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
      Teuchos::ParameterList& SchurCompList = mllist.sublist("coarse: SchurComp list");
      SchurCompList = TranslateSolverParameters(
          DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
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
const Teuchos::ParameterList LINALG::Solver::TranslateBACIToMuelu(
    const Teuchos::ParameterList& inparams, Teuchos::ParameterList* azlist)
{
  Teuchos::ParameterList muelulist;

  std::string xmlfile = inparams.get<std::string>("MUELU_XML_FILE");
  if (xmlfile != "none") muelulist.set("MUELU_XML_FILE", xmlfile);

  muelulist.set<bool>(
      "MUELU_XML_ENFORCE", DRT::INPUT::IntegralValue<bool>(inparams, "MUELU_XML_ENFORCE"));
  muelulist.set<bool>("LINALG::MueLu_Preconditioner", true);

  muelulist.set("aggregation: threshold", inparams.get<double>("ML_PROLONG_THRES"));

  int doRepart = DRT::INPUT::IntegralValue<int>(inparams, "MueLu_REBALANCE");
  if (doRepart > 2)
  {
    muelulist.set("repartition: enable", 1);
  }
  else
  {
    muelulist.set("repartition: enable", 0);
  }
  muelulist.set("repartition: partitioner", "ParMETIS");
  muelulist.set(
      "repartition: max min ratio", inparams.get<double>("MueLu_REBALANCE_NONZEROIMBALANCE"));
  muelulist.set("repartition: min per proc", inparams.get<int>("MueLu_REBALANCE_MINROWS"));
  muelulist.set("aggregation: min nodes per aggregate", inparams.get<int>("MueLu_MIN_AGG_SIZE"));

  return muelulist;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList LINALG::Solver::TranslateBACIToBelos(
    const Teuchos::ParameterList& inparams)
{
  Teuchos::ParameterList outparams;
  outparams.set("solver", "belos");
  outparams.set("symmetric", false);
  Teuchos::ParameterList& beloslist = outparams.sublist("Belos Parameters");

  //--------------------------------------------- set type of solver
  switch (DRT::INPUT::IntegralValue<INPAR::SOLVER::AzSolverType>(inparams, "AZSOLVE"))
  {
    case INPAR::SOLVER::azsolv_CG:
      beloslist.set("Solver Type", "CG");
      break;
    case INPAR::SOLVER::azsolv_GMRES:
      beloslist.set("Solver Type", "GMRES");
      break;
    case INPAR::SOLVER::azsolv_BiCGSTAB:
      beloslist.set("Solver Type", "BiCGSTAB");
      break;
    default:
    {
      dserror("Flag '%s'! \nUnknown solver for Belos.",
          DRT::INPUT::IntegralValue<INPAR::SOLVER::AzSolverType>(inparams, "AZSOLVE"));
      break;
    }
  }

  //------------------------------------- set type of preconditioner
  const int azprectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(inparams, "AZPREC");
  switch (azprectyp)
  {
    case INPAR::SOLVER::azprec_none:
      beloslist.set("Preconditioner Type", "none");
      break;
    case INPAR::SOLVER::azprec_ILUT:
      // using ifpack
      beloslist.set("Preconditioner Type", "ILUT");
      break;
    case INPAR::SOLVER::azprec_ILU:
      // using ifpack
      beloslist.set("Preconditioner Type", "ILU");
      break;
    case INPAR::SOLVER::azprec_Jacobi:
    case INPAR::SOLVER::azprec_SymmGaussSeidel:
    case INPAR::SOLVER::azprec_GaussSeidel:
      // using ifpack
      beloslist.set("Preconditioner Type", "point relaxation");
      break;
    case INPAR::SOLVER::azprec_DownwindGaussSeidel:
      // using ifpack
      beloslist.set("Preconditioner Type", "point relaxation");
      beloslist.set<bool>("downwinding", true);
      beloslist.set<double>("downwinding tau", inparams.get<double>("DWINDTAU"));
      std::cout << "warning: downwinding not supported?" << std::endl;
      break;
    case INPAR::SOLVER::azprec_LU:
      // using ifpack
      beloslist.set("Preconditioner Type", "Amesos");
      break;
    case INPAR::SOLVER::azprec_ML:
    case INPAR::SOLVER::azprec_MLfluid:
    case INPAR::SOLVER::azprec_MLfluid2:
    case INPAR::SOLVER::azprec_MueLuAMG:
      beloslist.set("Preconditioner Type", "ML");
      break;
    case INPAR::SOLVER::azprec_MueLuAMG_fluid:
      beloslist.set("Preconditioner Type", "Fluid");
      break;
    case INPAR::SOLVER::azprec_MueLuAMG_tsi:
      beloslist.set("Preconditioner Type", "TSI");
      break;
    case INPAR::SOLVER::azprec_MueLuAMG_contactSP:
      beloslist.set("Preconditioner Type", "ContactSP");
      break;
    case INPAR::SOLVER::azprec_BGS2x2:
      beloslist.set("Preconditioner Type", "ML");
      break;
      break;
    case INPAR::SOLVER::azprec_CheapSIMPLE:
      beloslist.set("Preconditioner Type", "CheapSIMPLE");
      break;
    default:
      dserror("Unknown preconditioner for Belos");
      break;
  }

  //------------------------------------- set other belos parameters
  beloslist.set("Num Blocks", inparams.get<int>("AZSUB"));
  // beloslist.set("Block Size", 1); // TODO blocksize
  beloslist.set("Orthogonalization", /*"DGKS"*/ /* ICGS,*/ "IMGS");
  beloslist.set("Maximum Iterations", inparams.get<int>("AZITER"));
  // beloslist.set("Adaptive Block Size",true);
  int outputfrequency = inparams.get<int>("AZOUTPUT");
  beloslist.set("Output Frequency", outputfrequency);
  int verbosity = inparams.get<int>("VERBOSITY");
  if (verbosity > 9)
    beloslist.set("Verbosity",
        Belos::Errors + Belos::Warnings + Belos::StatusTestDetails + Belos::TimingDetails);
  else if (verbosity > 4)
    beloslist.set("Verbosity", Belos::Errors + Belos::Warnings + Belos::StatusTestDetails);
  else if (verbosity > 2)
    beloslist.set("Verbosity", Belos::Errors + Belos::Warnings);
  else if (verbosity > 0)
    beloslist.set("Verbosity", Belos::Errors);
  // beloslist.set("allow permutation", DRT::INPUT::IntegralValue<int>(inparams,"PERMUTE_SYSTEM"));
  const int PermutationStrategy =
      DRT::INPUT::IntegralValue<INPAR::SOLVER::PermutationStrategy>(inparams, "PERMUTE_SYSTEM");

  switch (PermutationStrategy)
  {
    case INPAR::SOLVER::Permutation_algebraic:
      beloslist.set("permutation strategy", "Algebraic");
      break;
    case INPAR::SOLVER::Permutation_local:
      beloslist.set("permutation strategy", "Local");
      break;
    case INPAR::SOLVER::Permutation_none:
    default:
      beloslist.set("permutation strategy", "none");
      break;
  }
  // beloslist.set("allow permutation", bAllowPermutation);
  double nonDiagDominance = inparams.get<double>("NON_DIAGDOMINANCE_RATIO");
  beloslist.set("diagonal dominance ratio", nonDiagDominance);
  beloslist.set("Output Style", Belos::Brief);
  beloslist.set("Convergence Tolerance", inparams.get<double>("AZTOL"));

  //-------------------------------- set parameters for Ifpack if used
  if (azprectyp == INPAR::SOLVER::azprec_ILU || azprectyp == INPAR::SOLVER::azprec_ILUT ||
      azprectyp == INPAR::SOLVER::azprec_ICC || azprectyp == INPAR::SOLVER::azprec_LU ||
      azprectyp == INPAR::SOLVER::azprec_SymmGaussSeidel ||
      azprectyp == INPAR::SOLVER::azprec_GaussSeidel ||
      azprectyp == INPAR::SOLVER::azprec_DownwindGaussSeidel ||
      azprectyp == INPAR::SOLVER::azprec_Jacobi)
  {
    Teuchos::ParameterList& ifpacklist = outparams.sublist("IFPACK Parameters");
    ifpacklist = LINALG::Solver::TranslateBACIToIfpack(inparams);
  }

  //------------------------------------- set parameters for CheapSIMPLE if used
  if (azprectyp == INPAR::SOLVER::azprec_CheapSIMPLE)
  {
    Teuchos::ParameterList& simplelist = outparams.sublist("CheapSIMPLE Parameters");
    simplelist.set("Prec Type", "CheapSIMPLE");  // not used
    Teuchos::ParameterList& predictList = simplelist.sublist("Inverse1");
    predictList = TranslateSolverParameters(
        DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
    std::cout << predictList << std::endl;
    Teuchos::ParameterList& schurList = simplelist.sublist("Inverse2");
    schurList = TranslateSolverParameters(
        DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
  }

  //------------------------------------- set parameters for ML if used
  if (azprectyp == INPAR::SOLVER::azprec_ML || azprectyp == INPAR::SOLVER::azprec_MLfluid ||
      azprectyp == INPAR::SOLVER::azprec_MLfluid2)
  {
    Teuchos::ParameterList& mllist = outparams.sublist("ML Parameters");
    mllist = LINALG::Solver::TranslateBACIToML(inparams, &beloslist);
  }
  if (azprectyp == INPAR::SOLVER::azprec_MueLuAMG)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu Parameters");
    muelulist = LINALG::Solver::TranslateBACIToMuelu(inparams, &beloslist);
  }
  if (azprectyp == INPAR::SOLVER::azprec_MueLuAMG_fluid)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu (Fluid) Parameters");
    muelulist = LINALG::Solver::TranslateBACIToMuelu(inparams, &beloslist);
  }
  if (azprectyp == INPAR::SOLVER::azprec_MueLuAMG_tsi)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu (TSI) Parameters");
    muelulist = LINALG::Solver::TranslateBACIToMuelu(inparams, &beloslist);
  }
  if (azprectyp == INPAR::SOLVER::azprec_MueLuAMG_contactSP)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu (Contact) Parameters");
    muelulist = LINALG::Solver::TranslateBACIToMuelu(inparams, &beloslist);
  }
  if (azprectyp == INPAR::SOLVER::azprec_BGS2x2)
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

  return outparams;
}

/*------------------------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------------------------*/
const Teuchos::ParameterList LINALG::Solver::TranslateBACIToAztec(
    const Teuchos::ParameterList& inparams)
{
  Teuchos::ParameterList outparams;
  outparams.set("solver", "aztec");
  outparams.set("symmetric", false);
  Teuchos::ParameterList& azlist = outparams.sublist("Aztec Parameters");

  //--------------------------------- set scaling of linear problem
  const int azscal = DRT::INPUT::IntegralValue<int>(inparams, "AZSCAL");
  if (azscal == 1)
    azlist.set("scaling", "symmetric");
  else if (azscal == 2)
    azlist.set("scaling", "infnorm");
  else
    azlist.set("scaling", "none");

  //--------------------------------------------- set type of solver
  switch (DRT::INPUT::IntegralValue<INPAR::SOLVER::AzSolverType>(inparams, "AZSOLVE"))
  {
    case INPAR::SOLVER::azsolv_CG:
      azlist.set("AZ_solver", AZ_cg);
      break;
    case INPAR::SOLVER::azsolv_GMRES:
      azlist.set("AZ_solver", AZ_gmres);
      break;
    case INPAR::SOLVER::azsolv_GMRES_CONDEST:
      azlist.set("AZ_solver", AZ_gmres_condnum);
      break;
    case INPAR::SOLVER::azsolv_BiCGSTAB:
      azlist.set("AZ_solver", AZ_bicgstab);
      break;
    default:
    {
      dserror("Flag '%s'! Unknown solver for AztecOO",
          DRT::INPUT::IntegralValue<INPAR::SOLVER::AzSolverType>(inparams, "AZSOLVE"));
      break;
    }
  }

  //------------------------------------- set type of preconditioner
  const int azprectyp = DRT::INPUT::IntegralValue<INPAR::SOLVER::AzPrecType>(inparams, "AZPREC");
  switch (azprectyp)
  {
    case INPAR::SOLVER::azprec_none:
      azlist.set("AZ_precond", AZ_none);
      azlist.set("AZ_subdomain_solve", AZ_none);
      azlist.set("AZ_precond", AZ_none);
      azlist.set("Preconditioner Type", "none");
      break;
    case INPAR::SOLVER::azprec_ILUT:
      // using ifpack
      azlist.set("AZ_precond", AZ_user_precond);
      azlist.set("Preconditioner Type", "ILUT");
      break;
    case INPAR::SOLVER::azprec_ILU:
      // using ifpack
      azlist.set("AZ_precond", AZ_user_precond);
      azlist.set("Preconditioner Type", "ILU");
      break;
    case INPAR::SOLVER::azprec_Neumann:
      azlist.set("AZ_precond", AZ_Neumann);
      break;
    case INPAR::SOLVER::azprec_Least_Squares:
      azlist.set("AZ_precond", AZ_ls);
      break;
    case INPAR::SOLVER::azprec_Jacobi:
      // using ifpack
      azlist.set("AZ_precond", AZ_user_precond);
      azlist.set("Preconditioner Type", "point relaxation");
      break;
    case INPAR::SOLVER::azprec_Chebyshev:
      // using ifpack
      azlist.set("AZ_precond", AZ_user_precond);
      azlist.set("Preconditioner Type", "Chebyshev");
      break;
    case INPAR::SOLVER::azprec_SymmGaussSeidel:
      // using ifpack
      azlist.set("AZ_precond", AZ_user_precond);
      azlist.set("Preconditioner Type", "point relaxation");
      break;
    case INPAR::SOLVER::azprec_GaussSeidel:
      // using ifpack
      azlist.set("AZ_precond", AZ_user_precond);
      azlist.set("Preconditioner Type", "point relaxation");
      break;
    case INPAR::SOLVER::azprec_DownwindGaussSeidel:
      // using ifpack
      azlist.set("AZ_precond", AZ_user_precond);
      azlist.set("Preconditioner Type", "point relaxation");
      azlist.set<bool>("downwinding", true);
      azlist.set<double>("downwinding tau", inparams.get<double>("DWINDTAU"));
      break;
    case INPAR::SOLVER::azprec_LU:
      // using ifpack
      azlist.set("AZ_precond", AZ_user_precond);
      azlist.set("Preconditioner Type", "Amesos");
      break;
    case INPAR::SOLVER::azprec_RILU:
      azlist.set("AZ_precond", AZ_dom_decomp);
      azlist.set("AZ_subdomain_solve", AZ_rilu);
      azlist.set("AZ_graph_fill", inparams.get<int>("IFPACKGFILL"));
      break;
    case INPAR::SOLVER::azprec_ICC:
      // using ifpack
      azlist.set("AZ_precond", AZ_user_precond);
      azlist.set("Preconditioner Type", "IC");
      break;
    case INPAR::SOLVER::azprec_ML:
    case INPAR::SOLVER::azprec_MLfluid:
    case INPAR::SOLVER::azprec_MLfluid2:
    case INPAR::SOLVER::azprec_BGS2x2:
    case INPAR::SOLVER::azprec_MueLuAMG:
    case INPAR::SOLVER::azprec_AMGnxn:
      azlist.set("AZ_precond", AZ_user_precond);
      break;
    case INPAR::SOLVER::azprec_MueLuAMG_fluid:
      azlist.set("AZ_precond", AZ_user_precond);
      azlist.set("Preconditioner Type", "Fluid");
      break;
    case INPAR::SOLVER::azprec_MueLuAMG_tsi:
      azlist.set("AZ_precond", AZ_user_precond);
      azlist.set("Preconditioner Type", "TSI");
      break;
    case INPAR::SOLVER::azprec_MueLuAMG_contactSP:
      azlist.set("AZ_precond", AZ_user_precond);
      azlist.set("Preconditioner Type", "ContactSP");
      break;
    case INPAR::SOLVER::azprec_CheapSIMPLE:
      azlist.set("AZ_precond", AZ_user_precond);
      azlist.set("Preconditioner Type", "CheapSIMPLE");
      break;
    default:
      dserror("Unknown preconditioner for AztecOO");
      break;
  }

  //------------------------------------- set other aztec parameters
  azlist.set("AZ_kspace", inparams.get<int>("AZSUB"));
  azlist.set("AZ_max_iter", inparams.get<int>("AZITER"));
  azlist.set("AZ_overlap", inparams.get<int>("IFPACKOVERLAP"));
  azlist.set("AZ_type_overlap", AZ_symmetric);
  azlist.set("AZ_poly_ord", inparams.get<int>("AZPOLY"));
  const int azoutput = inparams.get<int>("AZOUTPUT");
  if (!azoutput)
    azlist.set("AZ_output", AZ_none);  // AZ_none AZ_all AZ_warnings AZ_last 10
  else
    azlist.set("AZ_output", azoutput);
  azlist.set("AZ_diagnostics", inparams.get<int>("AZBDIAG"));  // AZ_none AZ_all
  azlist.set("AZ_conv", DRT::INPUT::IntegralValue<int>(inparams, "AZCONV"));
  azlist.set("AZ_tol", inparams.get<double>("AZTOL"));
  azlist.set("AZ_drop", inparams.get<double>("AZDROP"));
  azlist.set("AZ_scaling", AZ_none);
  azlist.set("AZ_keep_info", 0);
  // set reuse parameters
  azlist.set("ncall", 0);                             // counting number of solver calls
  azlist.set("reuse", inparams.get<int>("AZREUSE"));  // reuse info for n solver calls
  // bool bAllowPermutation = DRT::INPUT::IntegralValue<bool>(inparams,"PERMUTE_SYSTEM");
  // azlist.set("allow permutation", bAllowPermutation);
  const int PermutationStrategy =
      DRT::INPUT::IntegralValue<INPAR::SOLVER::PermutationStrategy>(inparams, "PERMUTE_SYSTEM");

  switch (PermutationStrategy)
  {
    case INPAR::SOLVER::Permutation_algebraic:
      azlist.set("permutation strategy", "Algebraic");
      break;
    case INPAR::SOLVER::Permutation_local:
      azlist.set("permutation strategy", "Local");
      break;
    case INPAR::SOLVER::Permutation_none:
    default:
      azlist.set("permutation strategy", "none");
      break;
  }
  double nonDiagDominance = inparams.get<double>("NON_DIAGDOMINANCE_RATIO");
  azlist.set("diagonal dominance ratio", nonDiagDominance);
  azlist.set("verbosity", inparams.get<int>("VERBOSITY"));  // this is not an official Aztec flag

  //-------------------------------- set parameters for Ifpack if used
  if (azprectyp == INPAR::SOLVER::azprec_ILU || azprectyp == INPAR::SOLVER::azprec_ILUT ||
      azprectyp == INPAR::SOLVER::azprec_ICC || azprectyp == INPAR::SOLVER::azprec_LU ||
      azprectyp == INPAR::SOLVER::azprec_SymmGaussSeidel ||
      azprectyp == INPAR::SOLVER::azprec_GaussSeidel ||
      azprectyp == INPAR::SOLVER::azprec_DownwindGaussSeidel ||
      azprectyp == INPAR::SOLVER::azprec_Jacobi || azprectyp == INPAR::SOLVER::azprec_Chebyshev)
  {
    Teuchos::ParameterList& ifpacklist = outparams.sublist("IFPACK Parameters");
    ifpacklist = LINALG::Solver::TranslateBACIToIfpack(inparams);
  }

  //------------------------------------- set parameters for CheapSIMPLE if used
  if (azprectyp == INPAR::SOLVER::azprec_CheapSIMPLE)
  {
    Teuchos::ParameterList& simplelist = outparams.sublist("CheapSIMPLE Parameters");
    simplelist.set("Prec Type", "CheapSIMPLE");  // not used
    Teuchos::ParameterList& predictList = simplelist.sublist("Inverse1");
    predictList = TranslateSolverParameters(
        DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER1")));
    Teuchos::ParameterList& schurList = simplelist.sublist("Inverse2");
    schurList = TranslateSolverParameters(
        DRT::Problem::Instance()->SolverParams(inparams.get<int>("SUB_SOLVER2")));
  }

  //------------------------------------- set parameters for ML if used
  if (azprectyp == INPAR::SOLVER::azprec_ML || azprectyp == INPAR::SOLVER::azprec_MLfluid ||
      azprectyp == INPAR::SOLVER::azprec_MLfluid2)
  {
    Teuchos::ParameterList& mllist = outparams.sublist("ML Parameters");
    mllist = LINALG::Solver::TranslateBACIToML(inparams, &azlist);
  }
  if (azprectyp == INPAR::SOLVER::azprec_MueLuAMG)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu Parameters");
    muelulist = LINALG::Solver::TranslateBACIToMuelu(inparams, &azlist);
  }
  if (azprectyp == INPAR::SOLVER::azprec_MueLuAMG_fluid)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu (Fluid) Parameters");
    muelulist = LINALG::Solver::TranslateBACIToMuelu(inparams, &azlist);
  }
  if (azprectyp == INPAR::SOLVER::azprec_MueLuAMG_tsi)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu (TSI) Parameters");
    muelulist = LINALG::Solver::TranslateBACIToMuelu(inparams, &azlist);
  }
  if (azprectyp == INPAR::SOLVER::azprec_MueLuAMG_contactSP)
  {
    Teuchos::ParameterList& muelulist = outparams.sublist("MueLu (Contact) Parameters");
    muelulist = LINALG::Solver::TranslateBACIToMuelu(inparams, &azlist);
  }
  if (azprectyp == INPAR::SOLVER::azprec_BGS2x2)
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
  if (azprectyp == INPAR::SOLVER::azprec_AMGnxn)
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
const Teuchos::ParameterList LINALG::Solver::TranslateSolverParameters(
    const Teuchos::ParameterList& inparams)
{
  TEUCHOS_FUNC_TIME_MONITOR("LINALG::Solver:  0)   TranslateSolverParameters");

  Teuchos::ParameterList outparams;
  outparams.set<std::string>("name", inparams.get<std::string>("NAME"));

  switch (DRT::INPUT::IntegralValue<INPAR::SOLVER::SolverType>(inparams, "SOLVER"))
  {
    case INPAR::SOLVER::undefined:
      std::cout << "undefined solver! Set " << inparams.name() << "  in your dat file!"
                << std::endl;
      dserror("fix your dat file");
      break;
    case INPAR::SOLVER::superlu:
      outparams.set("solver", "superlu");
      outparams.set("symmetric", false);
      break;
    case INPAR::SOLVER::amesos_klu_sym:
      outparams.set("solver", "klu");
      outparams.set("symmetric", true);
      break;
    case INPAR::SOLVER::amesos_klu_nonsym:
      outparams.set("solver", "klu");
      outparams.set("symmetric", false);
      break;
    case INPAR::SOLVER::umfpack:
      outparams.set("solver", "umfpack");
      outparams.set("symmetric", false);
      break;
    case INPAR::SOLVER::lapack_sym:
      outparams.set("solver", "lapack");
      outparams.set("symmetric", true);
      break;
    case INPAR::SOLVER::lapack_nonsym:
      outparams.set("solver", "lapack");
      outparams.set("symmetric", false);
      break;
    case INPAR::SOLVER::belos:
      outparams = LINALG::Solver::TranslateBACIToBelos(inparams);
      break;
    case INPAR::SOLVER::aztec_msr:
      outparams = LINALG::Solver::TranslateBACIToAztec(inparams);
      break;
    default:
      dserror("Unsupported type of solver");
      break;
  }

  return outparams;
}
