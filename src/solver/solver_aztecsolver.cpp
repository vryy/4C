/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation to AztecOO solver package

\level 0

\maintainer Martin Kronbichler
*/
/*---------------------------------------------------------------------*/

#ifdef HAVE_MueLu

#include <MueLu_ConfigDefs.hpp>

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>

#include <MueLu.hpp>
#include <MueLu_FactoryBase.hpp>
#include <MueLu_PermutationFactory.hpp>
#include <MueLu_SmootherPrototype.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_DirectSolver.hpp>  // remove me
#include <MueLu_HierarchyHelpers.hpp>
#include <MueLu_VerboseObject.hpp>

#endif  // HAVE_MueLu

// Aztec headers
#include "AztecOO.h"
#include "AztecOO_StatusTestResNorm.h"
#include "AztecOO_StatusTestCombo.h"
#include "AztecOO_StatusTestMaxIters.h"

// BACI headers
#include "../drt_lib/drt_dserror.H"
#include "solver_aztecsolver.H"
#include "solver_aztecsolver_projectedresidual.H"
#include "../linalg/linalg_krylov_projector.H"

// Read a parameter value from a parameter list and copy it into a new parameter list (with another
// parameter name)
#define LINALG_COPY_PARAM(paramList, paramStr, varType, defaultValue, outParamList, outParamStr) \
  if (paramList.isParameter(paramStr))                                                           \
    outParamList.set<varType>(outParamStr, paramList.get<varType>(paramStr));                    \
  else                                                                                           \
    outParamList.set<varType>(outParamStr, defaultValue);

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::AztecSolver::AztecSolver(
    const Epetra_Comm& comm, Teuchos::ParameterList& params, FILE* outfile)
    : KrylovSolver(comm, params, outfile), numiters_(-1)
{
  ncall_ = 0;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
LINALG::SOLVER::AztecSolver::~AztecSolver()
{
  preconditioner_ = Teuchos::null;
  A_ = Teuchos::null;
  x_ = Teuchos::null;
  b_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
void LINALG::SOLVER::AztecSolver::Setup(Teuchos::RCP<Epetra_Operator> matrix,
    Teuchos::RCP<Epetra_MultiVector> x, Teuchos::RCP<Epetra_MultiVector> b, const bool refactor,
    const bool reset, Teuchos::RCP<LINALG::KrylovProjector> projector)
{
  if (!Params().isSublist("Aztec Parameters")) dserror("Do not have aztec parameter list");
  Teuchos::ParameterList& azlist = Params().sublist("Aztec Parameters");
  // int azoutput = azlist.get<int>("AZ_output",0);

  // see whether operator is a Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(matrix);

#ifdef HAVE_MueLu
  permutationStrategy_ = azlist.get<std::string>("permutation strategy", "none");
  diagDominanceRatio_ = azlist.get<double>("diagonal dominance ratio", 1.0);
  if (permutationStrategy_ == "none")
    bAllowPermutation_ = false;
  else
    bAllowPermutation_ = true;
#endif

  // decide whether we recreate preconditioners
  // after this call, the solver can access the preconditioner object using the
  // Preconditioner() function.
  int reuse = azlist.get("reuse", 0);
  const bool create = AllowReusePreconditioner(reuse, reset) == false;
  if (create)
  {
    ncall_ = 0;
    projector_ = projector;
    CreatePreconditioner(azlist, A != Teuchos::null, projector_);
  }

  // feed preconditioner with more information about linear system using
  // the "Linear System properties" sublist in the preconditioner's
  // paramter list
  if (Preconditioner() != NULL)
  {
    const std::string precondParamListName = Preconditioner()->getParameterListName();
    if (Params().isSublist(precondParamListName))
    {
      Teuchos::ParameterList& precondParams = Params().sublist(precondParamListName);
      Teuchos::ParameterList& linSystemProps = precondParams.sublist("Linear System properties");

      LINALG_COPY_PARAM(Params().sublist("Aztec Parameters").sublist("Linear System properties"),
          "contact slaveDofMap", Teuchos::RCP<Epetra_Map>, Teuchos::null, linSystemProps,
          "contact slaveDofMap");
      LINALG_COPY_PARAM(Params().sublist("Aztec Parameters").sublist("Linear System properties"),
          "contact masterDofMap", Teuchos::RCP<Epetra_Map>, Teuchos::null, linSystemProps,
          "contact masterDofMap");
      LINALG_COPY_PARAM(Params().sublist("Aztec Parameters").sublist("Linear System properties"),
          "contact innerDofMap", Teuchos::RCP<Epetra_Map>, Teuchos::null, linSystemProps,
          "contact innerDofMap");
      LINALG_COPY_PARAM(Params().sublist("Aztec Parameters").sublist("Linear System properties"),
          "contact activeDofMap", Teuchos::RCP<Epetra_Map>, Teuchos::null, linSystemProps,
          "contact activeDofMap");
      LINALG_COPY_PARAM(Params().sublist("Aztec Parameters").sublist("Linear System properties"),
          "ProblemType", std::string, "contact", linSystemProps, "ProblemType");
      LINALG_COPY_PARAM(Params().sublist("Aztec Parameters").sublist("Linear System properties"),
          "time step", int, -1, linSystemProps, "time step");
      LINALG_COPY_PARAM(Params().sublist("Aztec Parameters").sublist("Linear System properties"),
          "iter", int, -1, linSystemProps, "iter");
    }
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ////////////////////////////////////// permutation stuff
#ifdef HAVE_MueLu
  if (bAllowPermutation_)
  {
    // extract (user-given) additional information about linear system from
    // "Aztec Parameters" -> "Linear System properties"
    Teuchos::RCP<Epetra_Map> epSlaveDofMap =
        ExtractPermutationMap("Aztec Parameters", "contact slaveDofMap");

    // build permutation operators
    // permP, permQT and A = permQ^T A permP
    // all variables and information is stored in data_
    // note: we only allow permutations for rows in epSlaveDofMap
    //       the idea is not to disturb the matrix in regions which are
    //       known to work perfectly (no contact)
    BuildPermutationOperator(A, epSlaveDofMap);

    // TODO decide whether to permute linear system or not using the information of
    //      the permuted system matrix A

    // decide whether to permute linear system or not.
    // set all information corresponding to the decision.
    bPermuteLinearSystem_ = DecideAboutPermutation(A);
  }

  if (bAllowPermutation_ && bPermuteLinearSystem_)
  {
    // set
    // b_ = permP * b;
    // A_ = permQ^T * A * permP
    PermuteLinearSystem(A, b);

    // calculate (permQT)^T * b_f where b_f is the fine level null space (multi)vector
    // PermuteNullSpace(A);  // TODO think about this
    // do not permute null space to preserve pattern of null space for transfer operators
    // important e.g. for one pt aggregates?
  }
  else
  {
#endif  // HAVE_MueLu
    b_ = b;
    A_ = matrix;  // we cannot use A, since it could be Teuchos::null (for blocked operators)
#ifdef HAVE_MueLu
  }
#endif
  x_ = x;
  ////

#ifdef WRITEOUTSTATISTICS
  tttcreate.ResetStartTime();
#endif

  preconditioner_->Setup(create, &*A_, &*x_, &*b_);

#ifdef WRITEOUTSTATISTICS
  dtimeprecondsetup = tttcreate.ElapsedTime();
#endif
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
int LINALG::SOLVER::AztecSolver::Solve()
{
#ifdef WRITEOUTSTATISTICS
  Epetra_Time ttt(Comm());  // time measurement for whole routine
  ttt.ResetStartTime();
#endif

  Teuchos::ParameterList& azlist = Params().sublist("Aztec Parameters");

  // Allocate an aztec solver with default parameters
  // We do this every time because reusing the solver object
  // does lead to crashes that are not understood

  // create an aztec solver
  AztecOO aztec;
  aztec.SetAztecDefaults();

  // tell aztec to which stream to write
  aztec.SetOutputStream(std::cout);
  aztec.SetErrorStream(std::cerr);

  // Don't want linear problem to alter our aztec parameters (idiot feature!)
  // this is why we set our list here AFTER the linear problem has been set
  aztec.SetProblem(preconditioner_->LinearProblem());

  {
    // We don't want to use Aztec's scaling capabilities as we prefer to do
    // the scaling ourselves (so we precisely know what happens)
    // Therefore set scaling parameter to none and reset it after aztec has made
    // its internal copy of the parameter list
    std::string scaling = azlist.get("scaling", "none");
    azlist.set("scaling", "none");
    aztec.SetParameters(azlist, false);
    azlist.set("scaling", scaling);
  }

  aztec.SetPrecOperator(preconditioner_->PrecOperator());

  // iterate on the solution
  int iter = azlist.get("AZ_max_iter", 500);
  double tol = azlist.get("AZ_tol", 1.0e-6);

  // bool to return error code
  bool we_have_a_problem = false;

  // This hurts! It supresses error messages. This needs to be fixed.
  if (projector_ != Teuchos::null)
  {
    Epetra_Operator* op = aztec.GetProblem()->GetOperator();
    Epetra_Vector* rhs = static_cast<Epetra_Vector*>(aztec.GetProblem()->GetRHS());
    Epetra_Vector* lhs = static_cast<Epetra_Vector*>(aztec.GetProblem()->GetLHS());
    // max iterations
    aztest_maxiter_ = Teuchos::rcp(new AztecOO_StatusTestMaxIters(iter));
    // L2 norm of projected residual
    aztest_norm2_ =
        Teuchos::rcp(new AztecOO_StatusTestProjResNorm(*op, *lhs, *rhs, projector_, tol));
    aztest_norm2_->DefineResForm(
        AztecOO_StatusTestResNorm::Explicit, AztecOO_StatusTestResNorm::TwoNorm);
    aztest_norm2_->DefineScaleForm(
        AztecOO_StatusTestResNorm::NormOfInitRes, AztecOO_StatusTestResNorm::TwoNorm);
    // maxiters OR L2 norm of projected residual
    aztest_combo_ = Teuchos::rcp(new AztecOO_StatusTestCombo(AztecOO_StatusTestCombo::OR));
    aztest_combo_->AddStatusTest(*aztest_maxiter_);
    aztest_combo_->AddStatusTest(*aztest_norm2_);
    // set status test
    aztec.SetStatusTest(aztest_combo_.get());
  }
#if 0
  // create an aztec convergence test as combination of
  // L2-norm and Inf-Norm to be both satisfied where we demand
  // L2 < tol and Linf < 10*tol
  {
    Epetra_Operator* op  = aztec.GetProblem()->GetOperator();
    Epetra_Vector*   rhs = static_cast<Epetra_Vector*>(aztec.GetProblem()->GetRHS());
    Epetra_Vector*   lhs = static_cast<Epetra_Vector*>(aztec.GetProblem()->GetLHS());
    // max iterations
    aztest_maxiter_ = Teuchos::rcp(new AztecOO_StatusTestMaxIters(iter));
    // L2 norm
    aztest_norm2_ = Teuchos::rcp(new AztecOO_StatusTestResNorm(*op,*lhs,*rhs,tol));
    aztest_norm2_->DefineResForm(AztecOO_StatusTestResNorm::Implicit,
                                 AztecOO_StatusTestResNorm::TwoNorm);
    aztest_norm2_->DefineScaleForm(AztecOO_StatusTestResNorm::NormOfInitRes,
                                   AztecOO_StatusTestResNorm::TwoNorm);
    // Linf norm (demanded to be 1.0 times L2-norm now, to become an input parameter?)
    aztest_norminf_ = Teuchos::rcp(new AztecOO_StatusTestResNorm(*op,*lhs,*rhs,1.0*tol));
    aztest_norminf_->DefineResForm(AztecOO_StatusTestResNorm::Implicit,
                                   AztecOO_StatusTestResNorm::InfNorm);
    aztest_norminf_->DefineScaleForm(AztecOO_StatusTestResNorm::NormOfInitRes,
                                     AztecOO_StatusTestResNorm::InfNorm);
    // L2 AND Linf
    aztest_combo1_ = Teuchos::rcp(new AztecOO_StatusTestCombo(AztecOO_StatusTestCombo::SEQ));
    // maxiters OR (L2 AND Linf)
    aztest_combo2_ = Teuchos::rcp(new AztecOO_StatusTestCombo(AztecOO_StatusTestCombo::OR));
    aztest_combo1_->AddStatusTest(*aztest_norm2_);
    aztest_combo1_->AddStatusTest(*aztest_norminf_);
    aztest_combo2_->AddStatusTest(*aztest_maxiter_);
    aztest_combo2_->AddStatusTest(*aztest_combo1_);
    // set status test
    aztec.SetStatusTest(aztest_combo2_.get());
  }
#endif

  // if you want to get some information on eigenvalues of the Hessenberg matrix/the
  // estimated condition number of the preconditioned system, uncomment the following
  // line and set AZOUTPUT>0 in your .dat-file
  //  aztec.SetAztecOption(AZ_solver, AZ_gmres_condnum);

  //------------------------------- just do it----------------------------------------
  aztec.Iterate(iter, tol);
  //----------------------------------------------------------------------------------

  // store number of iterations
  numiters_ = aztec.NumIters();

  preconditioner_->Finish(&*A_, &*x_, &*b_);

  // check status of solution process
  const double* status = aztec.GetAztecStatus();
#if 0
  AztecOO_StatusType stat = aztest_combo2_->GetStatus();
  if (stat!=Converged)
  {
    bool resolve = false;
    if (stat==Unconverged)
    {
      if (comm_.MyPID()==0) printf("Max iterations reached in AztecOO\n");
    }
    else if (stat==Failed || stat==NaN || stat==PartialFailed)
    {
      if (comm_.MyPID()==0) printf("Numerical breakdown in AztecOO\n");
    }
    else dserror("Aztec returned unknown nonzero status %d",(int)stat);
  }
#else
  if (status[AZ_why] != AZ_normal)
  {
    we_have_a_problem = true;
    if (status[AZ_why] == AZ_breakdown)
    {
      if (comm_.MyPID() == 0) printf("Numerical breakdown in AztecOO\n");
    }
    else if (status[AZ_why] == AZ_ill_cond)
    {
      if (comm_.MyPID() == 0) printf("Problem is near singular in AztecOO\n");
    }
    else if (status[AZ_why] == AZ_loss)
    {
      if (comm_.MyPID() == 0) printf("Numerical loss of precision occurred in AztecOO\n");
    }
    else if (status[AZ_why] == AZ_maxits)
    {
      if (comm_.MyPID() == 0) printf("Max iterations reached in AztecOO\n");
    }
  }  // if (status[AZ_why] != AZ_normal)
#endif

#ifdef WRITEOUTSTATISTICS
  if (outfile_)
  {
    fprintf(outfile_,
        "LinIter %i\tNumGlobalElements %i\tAZ_solve_time %f\tAztecSolveTime %f\tAztecPrecondSetup "
        "%f\t\n",
        (int)status[AZ_its], A_->OperatorRangeMap().NumGlobalElements(), status[AZ_solve_time],
        dtimeprecondsetup_ + ttt.ElapsedTime(), dtimeprecondsetup_);
    fflush(outfile_);
  }
#endif

#ifdef HAVE_MueLu
  GlobalOrdinal rowperm = 0;
  GlobalOrdinal colperm = 0;
  GlobalOrdinal lrowperm = 0;
  GlobalOrdinal lcolperm = 0;
  int nonDiagDomRows = 0;
  int nonPermutedZeros = 0;
  int PermutedZeros = 0;
  int PermutedNearZeros = 0;
  int NonPermutedNearZeros = 0;

  if (bAllowPermutation_ && bPermuteLinearSystem_)
  {
    // repermutate solution vector
    this->ReTransformSolution();
    rowperm = data_->Get<GlobalOrdinal>("#RowPermutations", PermFact_.get());
    colperm = data_->Get<GlobalOrdinal>("#ColPermutations", PermFact_.get());
    lrowperm = data_->Get<GlobalOrdinal>("#WideRangeRowPermutations", PermFact_.get());
    lcolperm = data_->Get<GlobalOrdinal>("#WideRangeColPermutations", PermFact_.get());
  }
  if (data_->IsAvailable("nonDiagDomRows")) nonDiagDomRows = data_->Get<int>("nonDiagDomRows");
  if (data_->IsAvailable("NonPermutedZerosOnDiagonal"))
    nonPermutedZeros = data_->Get<int>("NonPermutedZerosOnDiagonal");
  if (data_->IsAvailable("PermutedZerosOnDiagonal"))
    PermutedZeros = data_->Get<int>("PermutedZerosOnDiagonal");
  if (data_->IsAvailable("PermutedNearZeros"))
    PermutedNearZeros = data_->Get<int>("PermutedNearZeros");
  if (data_->IsAvailable("NonPermutedNearZeros"))
    NonPermutedNearZeros = data_->Get<int>("NonPermutedNearZeros");

  // print some output if desired
  if (comm_.MyPID() == 0 && outfile_)
  {
    fprintf(outfile_,
        "AztecOO: "
        "unknowns/iterations/time/rowpermutations/colpermutations/lrowperm/lcolperm/nonDiagDomRows "
        "%d  %d  %f %d %d %d %d %d NonPermutedZeros/PermutedZeros %d %d bPermuted %d "
        "nonPermNearZeros/PermNearZeros %d %d\n",
        A_->OperatorRangeMap().NumGlobalElements(), (int)status[AZ_its], status[AZ_solve_time],
        rowperm, colperm, lrowperm, lcolperm, nonDiagDomRows, nonPermutedZeros, PermutedZeros,
        bPermuteLinearSystem_ ? 1 : 0, NonPermutedNearZeros, PermutedNearZeros);
    fflush(outfile_);
  }
#endif  // HAVE_MueLu

  ncall_ += 1;
  if (we_have_a_problem)
    return 1;
  else  // everything is fine
    return 0;
}
