#ifdef CCADISCRET

#include <vector>

#include <Epetra_Vector.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_VbrMatrix.h>
#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>

#include "fsi_nox_linearsystem_gcr.H"

NOX::FSI::LinearSystemGCR::LinearSystemGCR(
  Teuchos::ParameterList& printParams,
  Teuchos::ParameterList& linearSolverParams,
  const Teuchos::RefCountPtr< NOX::Epetra::Interface::Required>& iReq,
  const Teuchos::RefCountPtr< NOX::Epetra::Interface::Jacobian>& iJac,
  const Teuchos::RefCountPtr<Epetra_Operator>& jacobian,
  const NOX::Epetra::Vector& cloneVector,
  const Teuchos::RefCountPtr< NOX::Epetra::Scaling> s)
  : utils(printParams),
    jacInterfacePtr(iJac),
    jacType(EpetraOperator),
    jacPtr(jacobian),
    scaling(s),
    conditionNumberEstimate(0.0),
    timer(cloneVector.getEpetraVector().Comm()),
    timeApplyJacbianInverse(0.0)
{
  // Allocate solver
  tmpVectorPtr = Teuchos::rcp(new NOX::Epetra::Vector(cloneVector));

  // Jacobian operator is supplied
  jacType = getOperatorType(*jacPtr);

  reset(linearSolverParams);
}


void NOX::FSI::LinearSystemGCR::reset(Teuchos::ParameterList& linearSolverParams)
{
  zeroInitialGuess =
    linearSolverParams.get("Zero Initial Guess", false);

  manualScaling =
    linearSolverParams.get("Compute Scaling Manually", true);

  // Place linear solver details in the "Output" sublist of the
  // "Linear Solver" parameter list
  outputSolveDetails =
    linearSolverParams.get("Output Solver Details", true);
}


bool NOX::FSI::LinearSystemGCR::applyJacobian(const NOX::Epetra::Vector& input,
                                              NOX::Epetra::Vector& result) const
{
  jacPtr->SetUseTranspose(false);
  int status = jacPtr->Apply(input.getEpetraVector(), result.getEpetraVector());
  return status == 0;
}


bool NOX::FSI::LinearSystemGCR::applyJacobianTranspose(const NOX::Epetra::Vector& input,
                                                       NOX::Epetra::Vector& result) const
{
  jacPtr->SetUseTranspose(true);
  int status = jacPtr->Apply(input.getEpetraVector(), result.getEpetraVector());
  jacPtr->SetUseTranspose(false);

  return status == 0;
}


bool NOX::FSI::LinearSystemGCR::applyJacobianInverse(Teuchos::ParameterList &p,
                                                     const NOX::Epetra::Vector &input,
                                                     NOX::Epetra::Vector &result)
{
  double startTime = timer.WallTime();

  // Need non-const version of the input vector
  // Epetra_LinearProblem requires non-const versions so we can perform
  // scaling of the linear problem.
  NOX::Epetra::Vector& nonConstInput = const_cast< NOX::Epetra::Vector&>(input);

  // Zero out the delta X of the linear problem if requested by user.
  if (zeroInitialGuess)
    result.init(0.0);

  // Create Epetra linear problem object for the linear solve
  Epetra_LinearProblem Problem(jacPtr.get(),
  			       &(result.getEpetraVector()),
			       &(nonConstInput.getEpetraVector()));

  // ************* Begin linear system scaling *******************
  if ( !Teuchos::is_null(scaling) ) {

    if ( !manualScaling )
      scaling->computeScaling(Problem);

    scaling->scaleLinearSystem(Problem);

    if (utils.isPrintType(NOX::Utils::Details)) {
      utils.out() << *scaling << endl;
    }
  }
  // ************* End linear system scaling *******************

  // Get linear solver convergence parameters
  int maxit = p.get("Max Iterations", 400);
  double tol = p.get("Tolerance", 1.0e-6);

  int status = -1;

  // solve using GCR
  status = SolveGCR(input, result, maxit, tol);

  // Unscale the linear system
  if ( !Teuchos::is_null(scaling) )
    scaling->unscaleLinearSystem(Problem);

  // Set the output parameters in the "Output" sublist
  if (outputSolveDetails)
  {
    Teuchos::ParameterList& outputList = p.sublist("Output");
    int prevLinIters =
      outputList.get("Total Number of Linear Iterations", 0);
    int curLinIters = 0;
    double achievedTol = -1.0;
    //curLinIters = aztecSolverPtr->NumIters();
    //achievedTol = aztecSolverPtr->ScaledResidual();

    outputList.set("Number of Linear Iterations", curLinIters);
    outputList.set("Total Number of Linear Iterations", (prevLinIters + curLinIters));
    outputList.set("Achieved Tolerance", achievedTol);
  }

  double endTime = timer.WallTime();
  timeApplyJacbianInverse += (endTime - startTime);

  if (status != 0)
    return false;

  return true;
}


int NOX::FSI::LinearSystemGCR::SolveGCR(const NOX::Epetra::Vector &b,
                                        NOX::Epetra::Vector &x,
                                        int maxit, double tol)
{
  NOX::Epetra::Vector r(x, NOX::ShapeCopy);
  NOX::Epetra::Vector tmp(x, NOX::ShapeCopy);
  if (not zeroInitialGuess)
  {
    // calculate initial residual
    if (not applyJacobian(x, r))
      throwError("SolveGCR", "applyJacobian failed");
    r.update(1., b, -1.);
  }
  else
  {
    r = b;
  }

  double error = r.norm();
  if (error < tol)
  {
    utils.out() << "initial error=" << error << " tol=" << tol << endl;
    return 0;
  }

  std::vector<Teuchos::RefCountPtr< NOX::Epetra::Vector > > u;
  std::vector<Teuchos::RefCountPtr< NOX::Epetra::Vector > > c;
  int k=0;

  while (error>=tol)
  {
    // this is GCR, not GMRESR
    u.push_back(Teuchos::rcp(new NOX::Epetra::Vector(r)));
    if (not applyJacobian(r, tmp))
      throwError("SolveGCR", "applyJacobian failed");
    c.push_back(Teuchos::rcp(new NOX::Epetra::Vector(tmp)));

    for (int i=0; i<k; ++i)
    {
      double beta = c.back()->innerProduct(*c[i]);

      c.back()->update(-beta, *c[i], 1.);
      u.back()->update(-beta, *u[i], 1.);
    }

    double nc = c.back()->norm();
    c.back()->scale(1./nc);
    u.back()->scale(1./nc);

    double alpha = c.back()->innerProduct(r);
    x.update( alpha, *u.back(), 1.);
    r.update(-alpha, *c.back(), 1.);
    error = r.norm();

    k += 1;

    utils.out() << "gmresr |r|=" << error
                << " |dx|=" << u.back()->norm()*alpha
                << " tol=" << tol << endl;
  }

  return 0;
}


bool NOX::FSI::LinearSystemGCR::applyRightPreconditioning(bool useTranspose,
                                                          Teuchos::ParameterList& params,
                                                          const NOX::Epetra::Vector& input,
                                                          NOX::Epetra::Vector& result) const
{
  if (&result != &input)
    result = input;
  return true;
}


Teuchos::RefCountPtr< NOX::Epetra::Scaling> NOX::FSI::LinearSystemGCR::getScaling()
{
  return scaling;
}


void NOX::FSI::LinearSystemGCR::resetScaling(const Teuchos::RefCountPtr< NOX::Epetra::Scaling>& scalingObject)
{
  scaling = scalingObject;
}


bool NOX::FSI::LinearSystemGCR::computeJacobian(const NOX::Epetra::Vector& x)
{
  bool success = jacInterfacePtr->computeJacobian(x.getEpetraVector(),
						  *jacPtr);
  return success;
}


bool NOX::FSI::LinearSystemGCR::createPreconditioner(const NOX::Epetra::Vector& x,
                                      Teuchos::ParameterList& p,
                                      bool recomputeGraph) const
{
  return false;
}


bool NOX::FSI::LinearSystemGCR::destroyPreconditioner() const
{
  return false;
}


bool NOX::FSI::LinearSystemGCR::recomputePreconditioner(const NOX::Epetra::Vector& x,
                                                        Teuchos::ParameterList& linearSolverParams) const
{
  return false;
}


NOX::FSI::LinearSystemGCR::PreconditionerReusePolicyType NOX::FSI::LinearSystemGCR::getPreconditionerPolicy(bool advanceReuseCounter)
{
  return PRPT_REBUILD;
}


bool NOX::FSI::LinearSystemGCR::isPreconditionerConstructed() const
{
  return false;
}


bool NOX::FSI::LinearSystemGCR::hasPreconditioner() const
{
  return false;
}


Teuchos::RefCountPtr<const Epetra_Operator> NOX::FSI::LinearSystemGCR::getJacobianOperator() const
{
  return jacPtr;
}


Teuchos::RefCountPtr<Epetra_Operator> NOX::FSI::LinearSystemGCR::getJacobianOperator()
{
  return jacPtr;
}


Teuchos::RefCountPtr<const Epetra_Operator> NOX::FSI::LinearSystemGCR::getGeneratedPrecOperator() const
{
  return Teuchos::null;
}


Teuchos::RefCountPtr<Epetra_Operator> NOX::FSI::LinearSystemGCR::getGeneratedPrecOperator()
{
  return Teuchos::null;
}


void NOX::FSI::LinearSystemGCR::setJacobianOperatorForSolve(const Teuchos::RefCountPtr<const Epetra_Operator>& solveJacOp)
{
  jacPtr = Teuchos::rcp_const_cast<Epetra_Operator>(solveJacOp);
  jacType = getOperatorType(*solveJacOp);
}


void NOX::FSI::LinearSystemGCR::setPrecOperatorForSolve(const Teuchos::RefCountPtr<const Epetra_Operator>& solvePrecOp)
{
  throwError("setPrecOperatorForSolve", "no preconditioner supported");
}


void NOX::FSI::LinearSystemGCR::throwError(const string& functionName, const string& errorMsg) const
{
  if (utils.isPrintType(NOX::Utils::Error))
  {
    utils.out() << "NOX::FSI::LinearSystemGCR::" << functionName
                << " - " << errorMsg << endl;
  }
  throw "NOX Error";
}


NOX::FSI::LinearSystemGCR::OperatorType NOX::FSI::LinearSystemGCR::getOperatorType(const Epetra_Operator& Op)
{
  //***************
  //*** NOTE: The order in which the following tests occur is important!
  //***************

  const Epetra_Operator* testOperator = 0;

  // Is it an Epetra_CrsMatrix ?
  testOperator = dynamic_cast<const Epetra_CrsMatrix*>(&Op);
  if (testOperator != 0)
    return EpetraCrsMatrix;

  // Is it an Epetra_VbrMatrix ?
  testOperator = dynamic_cast<const Epetra_VbrMatrix*>(&Op);
  if (testOperator != 0)
    return EpetraVbrMatrix;

  // Is it an Epetra_RowMatrix ?
  testOperator = dynamic_cast<const Epetra_RowMatrix*>(&Op);
  if (testOperator != 0)
    return EpetraRowMatrix;

  // Otherwise it must be an Epetra_Operator!
  return EpetraOperator;
}


#endif
