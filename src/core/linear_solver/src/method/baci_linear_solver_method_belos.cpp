/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of interface to Belos solver package

\level 1

*/
/*---------------------------------------------------------------------*/

#include "baci_linear_solver_method_belos.hpp"

#include <BelosBiCGStabSolMgr.hpp>
#include <BelosBlockCGSolMgr.hpp>
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosConfigDefs.hpp>
#include <BelosEpetraAdapter.hpp>
#include <BelosLinearProblem.hpp>

BACI_NAMESPACE_OPEN

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
          "GLOBAL::ProblemType", "contact", linSystemProps, "GLOBAL::ProblemType");
      this->template copyParams<int>(
          this->Params().sublist("Belos Parameters").sublist("Linear System properties"),
          "time step", -1, linSystemProps, "time step");
      this->template copyParams<int>(
          this->Params().sublist("Belos Parameters").sublist("Linear System properties"), "iter",
          -1, linSystemProps, "iter");
    }
  }

  this->b_ = b;
  this->A_ =
      matrix;  // we cannot use A here, since it could be Teuchos::null (for blocked operators);
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

  Teuchos::RCP<Belos::LinearProblem<double, VectorType, MatrixType>> problem = Teuchos::rcp(
      new Belos::LinearProblem<double, VectorType, MatrixType>(this->A_, this->x_, this->b_));
  bool we_have_a_problem = false;

  if (this->preconditioner_ != Teuchos::null)
  {
    Teuchos::RCP<Belos::EpetraPrecOp> belosPrec =
        Teuchos::rcp(new Belos::EpetraPrecOp(this->preconditioner_->PrecOperator()));
    problem->setRightPrec(belosPrec);
  }

  bool set = problem->setProblem();
  if (set == false)
    dserror("CORE::LINEAR_SOLVER::BelosSolver: Iterative solver failed to set up correctly.");

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
    dserror("CORE::LINEAR_SOLVER::BelosSolver: Unknown iterative solver solver type chosen.");

  Belos::ReturnType ret = newSolver->solve();

  numiters_ = newSolver->getNumIters();

  if (this->preconditioner_ != Teuchos::null)
    this->preconditioner_->Finish(this->A_.get(), this->x_.get(), this->b_.get());

  int my_error = 0;
  if (ret != Belos::Converged)
  {
    my_error = 1;
    we_have_a_problem = true;
  }

  int glob_error = 0;
  this->comm_.SumAll(&my_error, &glob_error, 1);

  if (glob_error > 0 and this->comm_.MyPID() == 0)
    std::cout << std::endl
              << "CORE::LINEAR_SOLVER::BelosSolver: WARNING: Iterative solver did not converge!"
              << std::endl;

  this->ncall_ += 1;

  if (we_have_a_problem)
    return 1;
  else
    return 0;
}

//----------------------------------------------------------------------------------
//----------------------------------------------------------------------------------
// explicit initialization
template class CORE::LINEAR_SOLVER::BelosSolver<Epetra_Operator, Epetra_MultiVector>;

BACI_NAMESPACE_CLOSE
