/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation to AztecOO solver package

\level 0

*/
/*---------------------------------------------------------------------*/

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
#include <MueLu_HierarchyUtils.hpp>
#include <MueLu_VerboseObject.hpp>

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

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
LINALG::SOLVER::AztecSolver<MatrixType, VectorType>::AztecSolver(
    const Epetra_Comm& comm, Teuchos::ParameterList& params, FILE* outfile)
    : KrylovSolver<MatrixType, VectorType>(comm, params, outfile), numiters_(-1)
{
  this->ncall_ = 0;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
void LINALG::SOLVER::AztecSolver<MatrixType, VectorType>::Setup(Teuchos::RCP<MatrixType> matrix,
    Teuchos::RCP<VectorType> x, Teuchos::RCP<VectorType> b, const bool refactor, const bool reset,
    Teuchos::RCP<LINALG::KrylovProjector> projector)
{
  if (!this->Params().isSublist("Aztec Parameters")) dserror("Do not have aztec parameter list");
  Teuchos::ParameterList& azlist = this->Params().sublist("Aztec Parameters");

  // see whether operator is a Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(matrix);

  this->permutationStrategy_ = azlist.get<std::string>("permutation strategy", "none");
  this->diagDominanceRatio_ = azlist.get<double>("diagonal dominance ratio", 1.0);
  if (this->permutationStrategy_ == "none")
    this->bAllowPermutation_ = false;
  else
    this->bAllowPermutation_ = true;

  // decide whether we recreate preconditioners
  // after this call, the solver can access the preconditioner object using the
  // Preconditioner() function.
  int reuse = azlist.get("reuse", 0);
  const bool create = this->AllowReusePreconditioner(reuse, reset) == false;
  if (create)
  {
    this->ncall_ = 0;
    projector_ = projector;
    this->CreatePreconditioner(azlist, A != Teuchos::null, projector_);
  }

  // feed preconditioner with more information about linear system using
  // the "Linear System properties" sublist in the preconditioner's
  // paramter list
  {
    const std::string precondParamListName = this->Preconditioner().getParameterListName();
    if (this->Params().isSublist(precondParamListName))
    {
      Teuchos::ParameterList& precondParams = this->Params().sublist(precondParamListName);
      Teuchos::ParameterList& linSystemProps = precondParams.sublist("Linear System properties");
      this->template copyParams<Teuchos::RCP<Epetra_Map>>(
          this->Params().sublist("Aztec Parameters").sublist("Linear System properties"),
          "contact slaveDofMap", Teuchos::null, linSystemProps, "contact slaveDofMap");
      this->template copyParams<Teuchos::RCP<Epetra_Map>>(
          this->Params().sublist("Aztec Parameters").sublist("Linear System properties"),
          "contact masterDofMap", Teuchos::null, linSystemProps, "contact masterDofMap");
      this->template copyParams<Teuchos::RCP<Epetra_Map>>(
          this->Params().sublist("Aztec Parameters").sublist("Linear System properties"),
          "contact innerDofMap", Teuchos::null, linSystemProps, "contact innerDofMap");
      this->template copyParams<Teuchos::RCP<Epetra_Map>>(
          this->Params().sublist("Aztec Parameters").sublist("Linear System properties"),
          "contact activeDofMap", Teuchos::null, linSystemProps, "contact activeDofMap");
      this->template copyParams<std::string>(
          this->Params().sublist("Aztec Parameters").sublist("Linear System properties"),
          "ProblemType", "contact", linSystemProps, "ProblemType");
      this->template copyParams<int>(
          this->Params().sublist("Aztec Parameters").sublist("Linear System properties"),
          "time step", -1, linSystemProps, "time step");
      this->template copyParams<int>(
          this->Params().sublist("Aztec Parameters").sublist("Linear System properties"), "iter",
          -1, linSystemProps, "iter");
    }
  }

  ///////////////////////////////////////////////////////////////////////
  ///////////////////////////////////////////////////////////////////////
  ////////////////////////////////////// permutation stuff
  if (this->bAllowPermutation_)
  {
    // extract (user-given) additional information about linear system from
    // "Aztec Parameters" -> "Linear System properties"
    Teuchos::RCP<Epetra_Map> epSlaveDofMap =
        this->ExtractPermutationMap("Aztec Parameters", "contact slaveDofMap");

    // build permutation operators
    // permP, permQT and A = permQ^T A permP
    // all variables and information is stored in data_
    // note: we only allow permutations for rows in epSlaveDofMap
    //       the idea is not to disturb the matrix in regions which are
    //       known to work perfectly (no contact)
    this->BuildPermutationOperator(A, epSlaveDofMap);

    // TODO decide whether to permute linear system or not using the information of
    //      the permuted system matrix A

    // decide whether to permute linear system or not.
    // set all information corresponding to the decision.
    this->bPermuteLinearSystem_ = this->DecideAboutPermutation(A);
  }

  if (this->bAllowPermutation_ && this->bPermuteLinearSystem_)
  {
    // set
    // b_ = permP * b;
    // A_ = permQ^T * A * permP
    this->PermuteLinearSystem(A, b);

    // calculate (permQT)^T * b_f where b_f is the fine level null space (multi)vector
    // PermuteNullSpace(A);  // TODO think about this
    // do not permute null space to preserve pattern of null space for transfer operators
    // important e.g. for one pt aggregates?
  }
  else
  {
    this->b_ = b;
    this->A_ = matrix;  // we cannot use A, since it could be Teuchos::null (for blocked operators)
  }
  this->x_ = x;
  ////

#ifdef WRITEOUTSTATISTICS
  tttcreate.ResetStartTime();
#endif

  this->preconditioner_->Setup(create, this->A_.get(), this->x_.get(), this->b_.get());

#ifdef WRITEOUTSTATISTICS
  dtimeprecondsetup = tttcreate.ElapsedTime();
#endif
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
int LINALG::SOLVER::AztecSolver<MatrixType, VectorType>::Solve()
{
#ifdef WRITEOUTSTATISTICS
  Epetra_Time ttt(Comm());  // time measurement for whole routine
  ttt.ResetStartTime();
#endif

  Teuchos::ParameterList& azlist = this->Params().sublist("Aztec Parameters");

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
  aztec.SetProblem(this->preconditioner_->LinearProblem());

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

  aztec.SetPrecOperator(this->preconditioner_->PrecOperator());

  // iterate on the solution
  int iter = azlist.get("AZ_max_iter", 500);
  double tol = azlist.get("AZ_tol", 1.0e-6);

  // bool to return error code
  bool we_have_a_problem = false;

  // This hurts! It supresses error messages. This needs to be fixed.
  if (projector_ != Teuchos::null)
  {
    MatrixType* op = aztec.GetProblem()->GetOperator();
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

  //------------------------------- just do it----------------------------------------
  aztec.Iterate(iter, tol);
  //----------------------------------------------------------------------------------

  // store number of iterations
  numiters_ = aztec.NumIters();

  this->preconditioner_->Finish(this->A_.get(), this->x_.get(), this->b_.get());

  // check status of solution process
  const double* status = aztec.GetAztecStatus();
  if (status[AZ_why] != AZ_normal)
  {
    we_have_a_problem = true;
    if (status[AZ_why] == AZ_breakdown)
    {
      if (this->comm_.MyPID() == 0) printf("Numerical breakdown in AztecOO\n");
    }
    else if (status[AZ_why] == AZ_ill_cond)
    {
      if (this->comm_.MyPID() == 0) printf("Problem is near singular in AztecOO\n");
    }
    else if (status[AZ_why] == AZ_loss)
    {
      if (this->comm_.MyPID() == 0) printf("Numerical loss of precision occurred in AztecOO\n");
    }
    else if (status[AZ_why] == AZ_maxits)
    {
      if (this->comm_.MyPID() == 0) printf("Max iterations reached in AztecOO\n");
    }
  }

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

  GlobalOrdinal rowperm = 0;
  GlobalOrdinal colperm = 0;
  GlobalOrdinal lrowperm = 0;
  GlobalOrdinal lcolperm = 0;
  int nonDiagDomRows = 0;
  int nonPermutedZeros = 0;
  int PermutedZeros = 0;
  int PermutedNearZeros = 0;
  int NonPermutedNearZeros = 0;

  if (this->bAllowPermutation_ && this->bPermuteLinearSystem_)
  {
    // repermutate solution vector
    this->ReTransformSolution();
    rowperm = this->data_->template Get<GlobalOrdinal>("#RowPermutations", this->PermFact_.get());
    colperm = this->data_->template Get<GlobalOrdinal>("#ColPermutations", this->PermFact_.get());
    lrowperm = this->data_->template Get<GlobalOrdinal>(
        "#WideRangeRowPermutations", this->PermFact_.get());
    lcolperm = this->data_->template Get<GlobalOrdinal>(
        "#WideRangeColPermutations", this->PermFact_.get());
  }
  if (this->data_->IsAvailable("nonDiagDomRows"))
    nonDiagDomRows = this->data_->template Get<int>("nonDiagDomRows");
  if (this->data_->IsAvailable("NonPermutedZerosOnDiagonal"))
    nonPermutedZeros = this->data_->template Get<int>("NonPermutedZerosOnDiagonal");
  if (this->data_->IsAvailable("PermutedZerosOnDiagonal"))
    PermutedZeros = this->data_->template Get<int>("PermutedZerosOnDiagonal");
  if (this->data_->IsAvailable("PermutedNearZeros"))
    PermutedNearZeros = this->data_->template Get<int>("PermutedNearZeros");
  if (this->data_->IsAvailable("NonPermutedNearZeros"))
    NonPermutedNearZeros = this->data_->template Get<int>("NonPermutedNearZeros");

  // print some output if desired
  if (this->comm_.MyPID() == 0 && this->outfile_)
  {
    fprintf(this->outfile_,
        "AztecOO: "
        "unknowns/iterations/time/rowpermutations/colpermutations/lrowperm/lcolperm/nonDiagDomRows "
        "%d  %d  %f %d %d %d %d %d NonPermutedZeros/PermutedZeros %d %d bPermuted %d "
        "nonPermNearZeros/PermNearZeros %d %d\n",
        this->A_->OperatorRangeMap().NumGlobalElements(), (int)status[AZ_its],
        status[AZ_solve_time], rowperm, colperm, lrowperm, lcolperm, nonDiagDomRows,
        nonPermutedZeros, PermutedZeros, this->bPermuteLinearSystem_ ? 1 : 0, NonPermutedNearZeros,
        PermutedNearZeros);
    fflush(this->outfile_);
  }

  this->ncall_ += 1;
  if (we_have_a_problem)
    return 1;
  else  // everything is fine
    return 0;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// explicit initialization
template class LINALG::SOLVER::AztecSolver<Epetra_Operator, Epetra_MultiVector>;