/*----------------------------------------------------------------------*/
/*! \file

\brief FSI linear system interface to the nonlinear solver NOX

\level 3
*/

/*----------------------------------------------------------------------*/

#include <vector>

#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include "baci_linalg_serialdensematrix.H"
#include "baci_linalg_serialdensevector.H"
#include <Epetra_VbrMatrix.h>
#include <Epetra_Vector.h>

#include "baci_fsi_nox_linearsystem.H"
#include "baci_lib_globalproblem.H"
#include "baci_linear_solver_method_linalg.H"
#include "baci_linalg_blocksparsematrix.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<Epetra_Operator>& J, const NOX::Epetra::Vector& cloneVector,
    Teuchos::RCP<CORE::LINALG::Solver> solver, const Teuchos::RCP<NOX::Epetra::Scaling> s)
    : utils_(printParams),
      jacInterfacePtr_(iJac),
      jacType_(EpetraOperator),
      precType_(EpetraOperator),
      jacPtr_(J),
      scaling_(s),
      conditionNumberEstimate_(0.0),
      callcount_(0),
      solver_(solver),
      timer_("", true),
      timeApplyJacbianInverse_(0.0)
{
  tmpVectorPtr_ = Teuchos::rcp(new NOX::Epetra::Vector(cloneVector));

  jacType_ = getOperatorType(*jacPtr_);

  reset(linearSolverParams);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::LinearSystem::OperatorType NOX::FSI::LinearSystem::getOperatorType(
    const Epetra_Operator& Op)
{
  // check via dynamic cast, which type of Jacobian was broadcast
  const Epetra_Operator* testOperator = nullptr;

  testOperator = dynamic_cast<
      const CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>*>(&Op);
  if (testOperator != 0) return BlockSparseMatrix;

  testOperator = dynamic_cast<const CORE::LINALG::SparseMatrix*>(&Op);
  if (testOperator != 0) return SparseMatrix;

  testOperator = dynamic_cast<const Epetra_CrsMatrix*>(&Op);
  if (testOperator != 0) return EpetraCrsMatrix;

  testOperator = dynamic_cast<const Epetra_VbrMatrix*>(&Op);
  if (testOperator != 0) return EpetraVbrMatrix;

  testOperator = dynamic_cast<const Epetra_RowMatrix*>(&Op);
  if (testOperator != 0) return EpetraRowMatrix;

  return EpetraOperator;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearSystem::reset(Teuchos::ParameterList& linearSolverParams)
{
  zeroInitialGuess_ = linearSolverParams.get("Zero Initial Guess", false);
  manualScaling_ = linearSolverParams.get("Compute Scaling Manually", true);
  outputSolveDetails_ = linearSolverParams.get("Output Solver Details", true);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearSystem::applyJacobian(
    const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const
{
  jacPtr_->SetUseTranspose(false);
  int status = jacPtr_->Apply(input.getEpetraVector(), result.getEpetraVector());

  return status == 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearSystem::applyJacobianTranspose(
    const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const
{
  jacPtr_->SetUseTranspose(true);
  int status = jacPtr_->Apply(input.getEpetraVector(), result.getEpetraVector());
  jacPtr_->SetUseTranspose(false);

  return status == 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearSystem::applyJacobianInverse(
    Teuchos::ParameterList& p, const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result)
{
  // Zero out the delta X of the linear problem if requested by user.
  if (zeroInitialGuess_) result.init(0.0);

  const int maxit = p.get("Max Iterations", 30);
  const double tol = p.get("Tolerance", 1.0e-10);

  Teuchos::RCP<Epetra_Vector> fres = Teuchos::rcp(new Epetra_Vector(input.getEpetraVector()));
  Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(&(result.getEpetraVector()), false);

  // get the hopefully adaptive linear solver convergence tolerance
  solver_->Params()
      .sublist("Belos Parameters")
      .set("Convergence Tolerance", p.get("Tolerance", 1.0e-10));
  solver_->Solve(jacPtr_, disi, fres, true, callcount_ == 0);

  callcount_ += 1;

  // Set the output parameters in the "Output" sublist
  if (outputSolveDetails_)
  {
    Teuchos::ParameterList& outputList = p.sublist("Output");
    const int prevLinIters = outputList.get("Total Number of Linear Iterations", 0);
    const int curLinIters = maxit;
    const double achievedTol = tol;

    outputList.set("Number of Linear Iterations", curLinIters);
    outputList.set("Total Number of Linear Iterations", (prevLinIters + curLinIters));
    outputList.set("Achieved Tolerance", achievedTol);
  }

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearSystem::applyRightPreconditioning(bool useTranspose,
    Teuchos::ParameterList& params, const NOX::Epetra::Vector& input,
    NOX::Epetra::Vector& result) const
{
  if (&result != &input) result = input;

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::Scaling> NOX::FSI::LinearSystem::getScaling() { return scaling_; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearSystem::resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& scalingObject)
{
  scaling_ = scalingObject;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearSystem::computeJacobian(const NOX::Epetra::Vector& x)
{
  bool success = jacInterfacePtr_->computeJacobian(x.getEpetraVector(), *jacPtr_);
  return success;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearSystem::createPreconditioner(
    const NOX::Epetra::Vector& x, Teuchos::ParameterList& p, bool recomputeGraph) const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearSystem::destroyPreconditioner() const { return false; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearSystem::recomputePreconditioner(
    const NOX::Epetra::Vector& x, Teuchos::ParameterList& linearSolverParams) const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::LinearSystem::PreconditionerReusePolicyType
NOX::FSI::LinearSystem::getPreconditionerPolicy(bool advanceReuseCounter)
{
  return PRPT_REBUILD;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearSystem::isPreconditionerConstructed() const { return false; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearSystem::hasPreconditioner() const { return false; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Operator> NOX::FSI::LinearSystem::getJacobianOperator() const
{
  return jacPtr_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> NOX::FSI::LinearSystem::getJacobianOperator() { return jacPtr_; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Operator> NOX::FSI::LinearSystem::getGeneratedPrecOperator() const
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> NOX::FSI::LinearSystem::getGeneratedPrecOperator()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearSystem::setJacobianOperatorForSolve(
    const Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
  jacPtr_ = Teuchos::rcp_const_cast<Epetra_Operator>(solveJacOp);
  jacType_ = getOperatorType(*solveJacOp);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearSystem::setPrecOperatorForSolve(
    const Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
  throwError("setPrecOperatorForSolve", "no preconditioner supported");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearSystem::throwError(
    const std::string& functionName, const std::string& errorMsg) const
{
  if (utils_.isPrintType(NOX::Utils::Error))
    utils_.out() << "NOX::FSI::LinearSystem::" << functionName << " - " << errorMsg << std::endl;

  throw "NOX Error";
}
