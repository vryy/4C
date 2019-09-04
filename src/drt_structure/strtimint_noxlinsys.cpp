/*----------------------------------------------------------------------*/
/*! \file
\brief Use #NOX as non-linear solution technique for implicit
       structureal time integration

\maintainer Matthias Mayr

\level 3
*/

/*----------------------------------------------------------------------*/
/* macros */

/*----------------------------------------------------------------------*/
/* headers */
#include <vector>

#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_VbrMatrix.h>
#include <Epetra_Vector.h>

#include "strtimint_noxlinsys.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_blocksparsematrix.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::STR::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<Epetra_Operator>& J, const NOX::Epetra::Vector& cloneVector,
    Teuchos::RCP<LINALG::Solver> structure_solver, const Teuchos::RCP<NOX::Epetra::Scaling> s)
    : utils_(printParams),
      jacInterfacePtr_(iJac),
      jacType_(EpetraOperator),
      precType_(EpetraOperator),
      jacPtr_(J),
      scaling_(s),
      conditionNumberEstimate_(0.0),
      callcount_(0),
      structureSolver_(structure_solver),
      timer_(cloneVector.getEpetraVector().Comm()),
      timeApplyJacbianInverse_(0.0)
{
  tmpVectorPtr_ = Teuchos::rcp(new NOX::Epetra::Vector(cloneVector));

  // std::cout << "STRUCTURE SOLVER: " << *structureSolver_ << " " << structureSolver_ << std::endl;

  // Jacobian operator is supplied.
  // get type of it
  jacType_ = getOperatorType(*jacPtr_);

  reset(linearSolverParams);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::STR::LinearSystem::~LinearSystem() {}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::STR::LinearSystem::OperatorType NOX::STR::LinearSystem::getOperatorType(
    const Epetra_Operator& Op)
{
  // check per dynamik cast, which type of Jacobian was broadcast

  const Epetra_Operator* testOperator = 0;

  testOperator =
      dynamic_cast<const LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>*>(&Op);
  if (testOperator != 0) return BlockSparseMatrix;

  testOperator = dynamic_cast<const LINALG::SparseMatrix*>(&Op);
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
void NOX::STR::LinearSystem::reset(Teuchos::ParameterList& linearSolverParams)
{
  zeroInitialGuess_ = linearSolverParams.get("Zero Initial Guess", false);
  manualScaling_ = linearSolverParams.get("Compute Scaling Manually", true);
  outputSolveDetails_ = linearSolverParams.get("Output Solver Details", true);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::STR::LinearSystem::applyJacobian(
    const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const
{
  jacPtr_->SetUseTranspose(false);
  int status = jacPtr_->Apply(input.getEpetraVector(), result.getEpetraVector());
  return (status == 0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::STR::LinearSystem::applyJacobianTranspose(
    const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const
{
  jacPtr_->SetUseTranspose(true);
  int status = jacPtr_->Apply(input.getEpetraVector(), result.getEpetraVector());
  jacPtr_->SetUseTranspose(false);
  return (status == 0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::STR::LinearSystem::applyJacobianInverse(
    Teuchos::ParameterList& p, const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result)
{
  double startTime = timer_.WallTime();

  // Zero out the delta X of the linear problem if requested by user.
  if (zeroInitialGuess_) result.init(0.0);

  int maxit = p.get("Max Iterations", 30);
  double tol = p.get("Tolerance", 1.0e-10);

  // Structure
  if (jacType_ == SparseMatrix)
  {
    Teuchos::RCP<Epetra_Vector> fres = Teuchos::rcp(new Epetra_Vector(input.getEpetraVector()));
    Teuchos::RCP<Epetra_Vector> disi = Teuchos::rcp(&(result.getEpetraVector()), false);
    LINALG::SparseMatrix* J = dynamic_cast<LINALG::SparseMatrix*>(jacPtr_.get());
    structureSolver_->Solve(J->EpetraMatrix(), disi, fres, true, callcount_ == 0);
    callcount_ += 1;
  }
  else
  {
    dserror("Cannot deal with Epetra_Operator of type %d", jacType_);
  }

  // Set the output parameters in the "Output" sublist
  if (outputSolveDetails_)
  {
    Teuchos::ParameterList& outputList = p.sublist("Output");
    int prevLinIters = outputList.get("Total Number of Linear Iterations", 0);
    int curLinIters = maxit;
    double achievedTol = tol;

    outputList.set("Number of Linear Iterations", curLinIters);
    outputList.set("Total Number of Linear Iterations", (prevLinIters + curLinIters));
    outputList.set("Achieved Tolerance", achievedTol);
  }

  double endTime = timer_.WallTime();
  timeApplyJacbianInverse_ += (endTime - startTime);

  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::STR::LinearSystem::applyRightPreconditioning(bool useTranspose,
    Teuchos::ParameterList& params, const NOX::Epetra::Vector& input,
    NOX::Epetra::Vector& result) const
{
  if (&result != &input) result = input;
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::Scaling> NOX::STR::LinearSystem::getScaling() { return scaling_; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::STR::LinearSystem::resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& scalingObject)
{
  scaling_ = scalingObject;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::STR::LinearSystem::computeJacobian(const NOX::Epetra::Vector& x)
{
  bool success = jacInterfacePtr_->computeJacobian(x.getEpetraVector(), *jacPtr_);
  return success;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::STR::LinearSystem::createPreconditioner(
    const NOX::Epetra::Vector& x, Teuchos::ParameterList& p, bool recomputeGraph) const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::STR::LinearSystem::destroyPreconditioner() const { return false; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::STR::LinearSystem::recomputePreconditioner(
    const NOX::Epetra::Vector& x, Teuchos::ParameterList& linearSolverParams) const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::STR::LinearSystem::PreconditionerReusePolicyType
NOX::STR::LinearSystem::getPreconditionerPolicy(bool advanceReuseCounter)
{
  return PRPT_REBUILD;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::STR::LinearSystem::isPreconditionerConstructed() const { return false; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::STR::LinearSystem::hasPreconditioner() const { return false; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Operator> NOX::STR::LinearSystem::getJacobianOperator() const
{
  return jacPtr_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> NOX::STR::LinearSystem::getJacobianOperator() { return jacPtr_; }


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Operator> NOX::STR::LinearSystem::getGeneratedPrecOperator() const
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> NOX::STR::LinearSystem::getGeneratedPrecOperator()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::STR::LinearSystem::setJacobianOperatorForSolve(
    const Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
  jacPtr_ = Teuchos::rcp_const_cast<Epetra_Operator>(solveJacOp);
  jacType_ = getOperatorType(*solveJacOp);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::STR::LinearSystem::setPrecOperatorForSolve(
    const Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
  throwError("setPrecOperatorForSolve", "no preconditioner supported");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::STR::LinearSystem::throwError(
    const std::string& functionName, const std::string& errorMsg) const
{
  if (utils_.isPrintType(NOX::Utils::Error))

  {
    utils_.out() << "NOX::STR::LinearSystem::" << functionName << " - " << errorMsg << std::endl;
  }
  throw "NOX Error";
}

/*----------------------------------------------------------------------*/
