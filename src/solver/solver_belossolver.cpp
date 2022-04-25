/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of interface to Belos solver package

\level 1

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
#include <MueLu_DirectSolver.hpp>
#include <Trilinos_version.h>
#include <MueLu_HierarchyUtils.hpp>
#include <MueLu_VerboseObject.hpp>

// Belos headers
#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockCGSolMgr.hpp"
#include "BelosBlockGmresSolMgr.hpp"
#include "BelosBiCGStabSolMgr.hpp"

// BACI headers
#include "solver_belossolver.H"
#include "solver_pointpreconditioner.H"
#include "solver_blockpreconditioners.H"
#include "solver_krylovprojectionpreconditioner.H"
#include "solver_ifpackpreconditioner.H"
#include "solver_mlpreconditioner.H"

#include "../linalg/linalg_solver.H"

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
LINALG::SOLVER::BelosSolver<MatrixType, VectorType>::BelosSolver(
    const Epetra_Comm& comm, Teuchos::ParameterList& params, FILE* outfile)
    : KrylovSolver<MatrixType, VectorType>(comm, params, outfile)
{
  this->ncall_ = 0;
  this->preconditioner_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
void LINALG::SOLVER::BelosSolver<MatrixType, VectorType>::Setup(Teuchos::RCP<MatrixType> matrix,
    Teuchos::RCP<VectorType> x, Teuchos::RCP<VectorType> b, const bool refactor, const bool reset,
    Teuchos::RCP<LINALG::KrylovProjector> projector)
{
  // see whether operator is a Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(matrix);

  if (!this->Params().isSublist("Belos Parameters")) dserror("Do not have belos parameter list");
  Teuchos::ParameterList& belist = this->Params().sublist("Belos Parameters");

  this->permutationStrategy_ = belist.get<std::string>("permutation strategy", "none");
  this->diagDominanceRatio_ = belist.get<double>("diagonal dominance ratio", 1.0);
  if (this->permutationStrategy_ == "none")
    this->bAllowPermutation_ = false;
  else
    this->bAllowPermutation_ = true;

  int reuse = belist.get("reuse", 0);
  bool create = this->AllowReusePreconditioner(reuse, reset) == false;
  if (create)
  {
    this->ncall_ = 0;
    this->CreatePreconditioner(belist, A != Teuchos::null, projector);
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
          this->Params().sublist("Belos Parameters").sublist("Linear System properties"),
          "contact slaveDofMap", Teuchos::null, linSystemProps, "contact slaveDofMap");
      this->template copyParams<Teuchos::RCP<Epetra_Map>>(
          this->Params().sublist("Belos Parameters").sublist("Linear System properties"),
          "contact masterDofMap", Teuchos::null, linSystemProps, "contact masterDofMap");
      this->template copyParams<Teuchos::RCP<Epetra_Map>>(
          this->Params().sublist("Belos Parameters").sublist("Linear System properties"),
          "contact innerDofMap", Teuchos::null, linSystemProps, "contact innerDofMap");
      this->template copyParams<Teuchos::RCP<Epetra_Map>>(
          this->Params().sublist("Belos Parameters").sublist("Linear System properties"),
          "contact activeDofMap", Teuchos::null, linSystemProps, "contact activeDofMap");
      this->template copyParams<std::string>(
          this->Params().sublist("Belos Parameters").sublist("Linear System properties"),
          "ProblemType", "contact", linSystemProps, "ProblemType");
      this->template copyParams<int>(
          this->Params().sublist("Belos Parameters").sublist("Linear System properties"),
          "time step", -1, linSystemProps, "time step");
      this->template copyParams<int>(
          this->Params().sublist("Belos Parameters").sublist("Linear System properties"), "iter",
          -1, linSystemProps, "iter");
    }
  }

  ////////////////////////////////////// permutation stuff
  if (this->bAllowPermutation_)
  {
    // extract (user-given) additional information about linear system from
    // "Aztec Parameters" -> "Linear System properties"
    Teuchos::RCP<Epetra_Map> epSlaveDofMap =
        this->ExtractPermutationMap("Belos Parameters", "contact slaveDofMap");

    // build permutation operators
    // permP, permQT and A = permQ^T A permP
    // all variables and information is stored in data_
    // note: we only allow permutations for rows in epSlaveDofMap
    //       the idea is not to disturb the matrix in regions which are
    //       known to work perfectly (no contact)
    this->BuildPermutationOperator(A, epSlaveDofMap);

    // decide whether to permute linear system or not.
    // set all information corresponding to the decision.
    this->bPermuteLinearSystem_ = this->DecideAboutPermutation(A);
  }

  // set linear system
  if (this->bAllowPermutation_ && this->bPermuteLinearSystem_)
  {
    // set
    // b_ = permP * b;
    // A_ = permQ^T * A * permP
    this->PermuteLinearSystem(A, b);
  }
  else
  {
    this->b_ = b;
    this->A_ =
        matrix;  // we cannot use A here, since it could be Teuchos::null (for blocked operators);
  }
  this->x_ = x;


  // call setup of preconditioner
  this->preconditioner_->Setup(create, this->A_.get(), this->x_.get(), this->b_.get());
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
int LINALG::SOLVER::BelosSolver<MatrixType, VectorType>::Solve()
{
  Teuchos::ParameterList& belist = this->Params().sublist("Belos Parameters");

  // build Belos linear problem
  Teuchos::RCP<Belos::LinearProblem<double, VectorType, MatrixType>> problem = Teuchos::rcp(
      new Belos::LinearProblem<double, VectorType, MatrixType>(this->A_, this->x_, this->b_));
  bool we_have_a_problem = false;
  // TODO support for left preconditioner?
  if (this->preconditioner_ != Teuchos::null)
  {
    // prepare preconditioner in preconditioner_->PrecOperator() for Belos
    Teuchos::RCP<Belos::EpetraPrecOp> belosPrec = Teuchos::rcp(
        new Belos::EpetraPrecOp(Teuchos::rcp(this->preconditioner_->PrecOperator(), false)));
    problem->setRightPrec(belosPrec);
  }
  bool set = problem->setProblem();
  if (set == false)
  {
    std::cout << std::endl
              << "ERROR: Belos::LinearProblem failed to set up correctly!" << std::endl;
  }

  // create iterative solver manager
  Teuchos::RCP<Belos::SolverManager<double, VectorType, MatrixType>> newSolver;
  std::string solverType = belist.get<std::string>("Solver Type");
  if (solverType == "GMRES")
    newSolver = Teuchos::rcp(new Belos::BlockGmresSolMgr<double, VectorType, MatrixType>(
        problem, Teuchos::rcp(&belist, false)));
  else if (solverType == "CG")
    newSolver = Teuchos::rcp(new Belos::BlockCGSolMgr<double, VectorType, MatrixType>(
        problem, Teuchos::rcp(&belist, false)));
  else if (solverType == "BiCGSTAB")
    newSolver = Teuchos::rcp(new Belos::BiCGStabSolMgr<double, VectorType, MatrixType>(
        problem, Teuchos::rcp(&belist, false)));
  else
    dserror("unknown solver type for Belos");

  //
  // Perform solve
  //
  Belos::ReturnType ret = newSolver->solve();

  // TODO: check me -> access solution x from linear problem???
  if (this->preconditioner_ != Teuchos::null)
    this->preconditioner_->Finish(this->A_.get(), this->x_.get(), this->b_.get());

  if (ret != Belos::Converged)
  {
    std::cout << std::endl << "WARNING: Belos did not converge!" << std::endl;
    we_have_a_problem = true;
  }

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
        "Belos: "
        "unknowns/iterations/time/rowpermutations/colpermutations/lrowperm/lcolperm/nonDiagDomRows "
        "%d  %d  %f %d %d %d %d %d NonPermutedZeros/PermutedZeros %d %d bPermuted %d "
        "nonPermNearZeros/PermNearZeros %d %d\n",
        this->A_->OperatorRangeMap().NumGlobalElements(), (int)newSolver->getNumIters(), -1.0,
        rowperm, colperm, lrowperm, lcolperm, nonDiagDomRows, nonPermutedZeros, PermutedZeros,
        this->bPermuteLinearSystem_ ? 1 : 0, NonPermutedNearZeros, PermutedNearZeros);
    fflush(this->outfile_);
  }

  this->ncall_ += 1;  // increment counter of solver calls
  if (we_have_a_problem)
    return 1;
  else  // everything is fine
    return 0;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// explicit initialization
template class LINALG::SOLVER::BelosSolver<Epetra_Operator, Epetra_MultiVector>;
