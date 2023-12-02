/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of interface to Belos solver package

\level 1

*/
/*---------------------------------------------------------------------*/

#include "baci_linear_solver_method_belos.H"

#include "baci_linear_solver_method_linalg.H"
#include "baci_linear_solver_preconditioner_block.H"
#include "baci_linear_solver_preconditioner_ifpack.H"
#include "baci_linear_solver_preconditioner_krylovprojection.H"
#include "baci_linear_solver_preconditioner_ml.H"
#include "baci_linear_solver_preconditioner_point.H"

#include <BelosBiCGStabSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosConfigDefs.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosLinearProblem.hpp>
#include <MueLu.hpp>
#include <MueLu_ConfigDefs.hpp>
#include <MueLu_DirectSolver.hpp>
#include <MueLu_FactoryBase.hpp>
#include <MueLu_HierarchyUtils.hpp>
#include <MueLu_PermutationFactory.hpp>
#include <MueLu_SmootherFactory.hpp>
#include <MueLu_SmootherPrototype.hpp>
#include <MueLu_VerboseObject.hpp>
#include <Trilinos_version.h>
#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
CORE::LINEAR_SOLVER::BelosSolver<MatrixType, VectorType>::BelosSolver(
    const Epetra_Comm& comm, Teuchos::ParameterList& params)
    : KrylovSolver<MatrixType, VectorType>(comm, params), numiters_(-1)
{
  this->ncall_ = 0;
  this->preconditioner_ = Teuchos::null;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
template <class MatrixType, class VectorType>
void CORE::LINEAR_SOLVER::BelosSolver<MatrixType, VectorType>::Setup(
    Teuchos::RCP<MatrixType> matrix, Teuchos::RCP<VectorType> x, Teuchos::RCP<VectorType> b,
    const bool refactor, const bool reset, Teuchos::RCP<CORE::LINALG::KrylovProjector> projector)
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
    // extract (user-given) additional information about linear system
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
int CORE::LINEAR_SOLVER::BelosSolver<MatrixType, VectorType>::Solve()
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

  // store number of iterations
  numiters_ = newSolver->getNumIters();

  // TODO: check me -> access solution x from linear problem???
  if (this->preconditioner_ != Teuchos::null)
    this->preconditioner_->Finish(this->A_.get(), this->x_.get(), this->b_.get());

  // communicate non convergence to proc 0 and print warning
  int my_error = 0;
  if (ret != Belos::Converged)
  {
    my_error = 1;
    we_have_a_problem = true;
  }

  int glob_error = 0;
  this->comm_.SumAll(&my_error, &glob_error, 1);

  if (glob_error > 0 and this->comm_.MyPID() == 0)
    std::cout << std::endl << "WARNING: Belos did not converge!" << std::endl;

  if (this->bAllowPermutation_ && this->bPermuteLinearSystem_)
  {
    // repermutate solution vector
    this->ReTransformSolution();
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
template class CORE::LINEAR_SOLVER::BelosSolver<Epetra_Operator, Epetra_MultiVector>;
