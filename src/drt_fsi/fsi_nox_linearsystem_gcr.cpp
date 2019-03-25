/*----------------------------------------------------------------------*/
/*!
\file fsi_nox_linearsystem_gcr.cpp

\brief Generalized conjugate residual linear system solver for FSI

\level 1

\maintainer Matthias Mayr
*/
/*----------------------------------------------------------------------*/

#include <vector>

#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_VbrMatrix.h>
#include <Epetra_Vector.h>

#include "../linalg/linalg_serialdensevector.H"
#include "../linalg/linalg_serialdensematrix.H"
#include "fsi_nox_linearsystem_gcr.H"


NOX::FSI::LinearSystemGCR::LinearSystemGCR(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<Epetra_Operator>& jacobian, const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> s)
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
  zeroInitialGuess = linearSolverParams.get("Zero Initial Guess", false);

  manualScaling = linearSolverParams.get("Compute Scaling Manually", true);

  // Place linear solver details in the "Output" sublist of the
  // "Linear Solver" parameter list
  outputSolveDetails = linearSolverParams.get("Output Solver Details", true);

  // so we have a new time step and start anew
  u_.clear();
  c_.clear();
}


bool NOX::FSI::LinearSystemGCR::applyJacobian(
    const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const
{
  jacPtr->SetUseTranspose(false);
  int status = jacPtr->Apply(input.getEpetraVector(), result.getEpetraVector());
  return status == 0;
}


bool NOX::FSI::LinearSystemGCR::applyJacobianTranspose(
    const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const
{
  jacPtr->SetUseTranspose(true);
  int status = jacPtr->Apply(input.getEpetraVector(), result.getEpetraVector());
  jacPtr->SetUseTranspose(false);

  return status == 0;
}


bool NOX::FSI::LinearSystemGCR::applyJacobianInverse(
    Teuchos::ParameterList& p, const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result)
{
  double startTime = timer.WallTime();

  // Need non-const version of the input vector
  // Epetra_LinearProblem requires non-const versions so we can perform
  // scaling of the linear problem.
  NOX::Epetra::Vector& nonConstInput = const_cast<NOX::Epetra::Vector&>(input);

  // Zero out the delta X of the linear problem if requested by user.
  if (zeroInitialGuess) result.init(0.0);

  // Create Epetra linear problem object for the linear solve
  Epetra_LinearProblem Problem(
      jacPtr.get(), &(result.getEpetraVector()), &(nonConstInput.getEpetraVector()));

  // ************* Begin linear system scaling *******************
  if (!Teuchos::is_null(scaling))
  {
    if (!manualScaling) scaling->computeScaling(Problem);

    scaling->scaleLinearSystem(Problem);

    if (utils.isPrintType(NOX::Utils::Details))
    {
      utils.out() << *scaling << std::endl;
    }
  }
  // ************* End linear system scaling *******************

  // Get linear solver convergence parameters
  int maxit = p.get("Max Iterations", 400);
  double tol = p.get("Tolerance", 1.0e-6);

  int status = -1;

  // solve using GCR
  std::string linearSolver = p.get("Solver", "GMRES");
  if (linearSolver == "GMRES")
    status = SolveGMRES(input, result, maxit, tol, p.get("Size of Krylov Subspace", 300));
  else if (linearSolver == "GCR")
    status = SolveGCR(input, result, maxit, tol);
  else
  {
    utils.out() << "ERROR: NOX::FSI::LinearSystemGCR::applyJacobianInverse" << std::endl
                << "\"Solver\" parameter \"" << linearSolver << "\" is invalid!" << std::endl;
    throw "NOX Error";
  }

  // Unscale the linear system
  if (!Teuchos::is_null(scaling)) scaling->unscaleLinearSystem(Problem);

  // Set the output parameters in the "Output" sublist
  if (outputSolveDetails)
  {
    Teuchos::ParameterList& outputList = p.sublist("Output");
    int prevLinIters = outputList.get("Total Number of Linear Iterations", 0);
    int curLinIters = maxit;
    double achievedTol = tol;

    outputList.set("Number of Linear Iterations", curLinIters);
    outputList.set("Total Number of Linear Iterations", (prevLinIters + curLinIters));
    outputList.set("Achieved Tolerance", achievedTol);
  }

  double endTime = timer.WallTime();
  timeApplyJacbianInverse += (endTime - startTime);

  if (status != 0) return false;

  return true;
}


int NOX::FSI::LinearSystemGCR::SolveGCR(
    const NOX::Epetra::Vector& b, NOX::Epetra::Vector& x, int& maxit, double& tol)
{
  NOX::Epetra::Vector r(x, NOX::ShapeCopy);
  NOX::Epetra::Vector tmp(x, NOX::ShapeCopy);
  if (not zeroInitialGuess)
  {
    // calculate initial residual
    if (not applyJacobian(x, r)) throwError("SolveGCR", "applyJacobian failed");
    r.update(1., b, -1.);
  }
  else
  {
    r = b;
  }

  double normb = b.norm();
  double error0 = r.norm() / normb;

  std::vector<Teuchos::RCP<NOX::Epetra::Vector>>& u = u_;
  std::vector<Teuchos::RCP<NOX::Epetra::Vector>>& c = c_;

  // reset krylov space
  u.clear();
  c.clear();

  int k = static_cast<int>(u.size());

  // use the available vectors
  for (int i = 0; i < k; ++i)
  {
    double alpha = c[i]->innerProduct(r);
    x.update(alpha, *u[i], 1.);
    r.update(-alpha, *c[i], 1.);
  }

  double error = 1. * normb;
  // while (error>=tol*error0)
  while (error / normb >= tol)
  {
    // this is GCR, not GMRESR
    u.push_back(Teuchos::rcp(new NOX::Epetra::Vector(r)));
    if (not applyJacobian(r, tmp)) throwError("SolveGCR", "applyJacobian failed");
    c.push_back(Teuchos::rcp(new NOX::Epetra::Vector(tmp)));

    for (int i = 0; i < k; ++i)
    {
      double beta = c.back()->innerProduct(*c[i]);

      c.back()->update(-beta, *c[i], 1.);
      u.back()->update(-beta, *u[i], 1.);
    }

    double nc = c.back()->norm();
    c.back()->scale(1. / nc);
    u.back()->scale(1. / nc);

    double alpha = c.back()->innerProduct(r);
    x.update(alpha, *u.back(), 1.);
    r.update(-alpha, *c.back(), 1.);
    error = r.norm();

    k += 1;

    utils.out() << "gcr |r|=" << error << " |r0|=" << error0 << " |dx|=" << u.back()->norm() * alpha
                << " |b|=" << normb << " tol=" << tol << std::endl;
  }

  maxit = k;
  tol = error;

  return 0;
}


int NOX::FSI::LinearSystemGCR::SolveGMRES(
    const NOX::Epetra::Vector& b, NOX::Epetra::Vector& x, int& max_iter, double& tol, int m)
{
  double resid = 0;
  LINALG::SerialDenseVector s(m + 1, true);
  LINALG::SerialDenseVector cs(m + 1, true);
  LINALG::SerialDenseVector sn(m + 1, true);
  LINALG::SerialDenseMatrix H(m + 1, m, true);

  NOX::Epetra::Vector r(x, NOX::ShapeCopy);
  NOX::Epetra::Vector w(x, NOX::ShapeCopy);
  if (not zeroInitialGuess)
  {
    // calculate initial residual
    if (not applyJacobian(x, r)) throwError("SolveGMRES", "applyJacobian failed");
    r.update(1., b, -1.);
  }
  else
  {
    r = b;
  }

  double normb = b.norm();
  double beta = r.norm();

  if (normb == 0.0) normb = 1;

  if ((resid = beta / normb) <= tol)
  {
    tol = resid;
    max_iter = 0;
    return 0;
  }

  std::vector<Teuchos::RCP<NOX::Epetra::Vector>> v;
  v.reserve(m + 1);

  int j = 1;
  while (j <= max_iter)
  {
    v.clear();
    v.push_back(Teuchos::rcp(new NOX::Epetra::Vector(r, NOX::ShapeCopy)));
    v[0]->update(1. / beta, r, 0.);
    s.Scale(0.0);
    s(0) = beta;

    for (int i = 0; i < m and j <= max_iter; i++, j++)
    {
      Epetra_Time t(x.getEpetraVector().Comm());
      // w = M.solve(A * v[i]);
      if (not applyJacobian(*v[i], w)) throwError("SolveGMRES", "applyJacobian failed");
      for (int k = 0; k <= i; k++)
      {
        H(k, i) = w.innerProduct(*v[k]);
        w.update(H(k, i), *v[k], -1.);
      }
      H(i + 1, i) = w.norm();
      v.push_back(Teuchos::rcp(new NOX::Epetra::Vector(w)));
      v.back()->scale(1.0 / H(i + 1, i));

      for (int k = 0; k < i; k++) ApplyPlaneRotation(H(k, i), H(k + 1, i), cs(k), sn(k));

      GeneratePlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
      ApplyPlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
      ApplyPlaneRotation(s(i), s(i + 1), cs(i), sn(i));

      utils.out() << "gmres |r|=" << std::scientific << fabs(s(i + 1))
                  << "   |b|=" << std::scientific << normb << "   tol=" << std::scientific << tol
                  << "   time=" << std::scientific << t.ElapsedTime() << std::endl;

      if ((resid = fabs(s(i + 1)) / normb) < tol)
      {
        // Update(x, i, H, s, v);
        LINALG::SerialDenseVector y(s);

        // Backsolve:
        for (int l = i; l >= 0; l--)
        {
          y(l) /= H(l, l);
          for (int k = l - 1; k >= 0; k--) y(k) -= H(k, l) * y(l);
        }

        for (int k = 0; k <= i; k++) x.update(y(k), *v[k], 1.);

        tol = resid;
        max_iter = j;
        return 0;
      }
    }

    // Update(x, m - 1, H, s, v);
    LINALG::SerialDenseVector y(s);

    // Backsolve:
    for (int i = m - 1; i >= 0; i--)
    {
      y(i) /= H(i, i);
      for (int k = i - 1; k >= 0; k--) y(k) -= H(k, i) * y(i);
    }

    for (int k = 0; k <= m - 1; k++) x.update(y(k), *v[k], 1.);

    // Isn't there a cheaper way to calculate that?
    if (not applyJacobian(x, r)) throwError("SolveGMRES", "applyJacobian failed");
    r.update(1., b, -1.);
    beta = r.norm();
    if ((resid = beta / normb) < tol)
    {
      tol = resid;
      max_iter = j;
      return 0;
    }
  }

  tol = resid;
  return 1;
}


void NOX::FSI::LinearSystemGCR::GeneratePlaneRotation(
    double& dx, double& dy, double& cs, double& sn)
{
  if (dy == 0.0)
  {
    cs = 1.0;
    sn = 0.0;
  }
  else if (fabs(dy) > fabs(dx))
  {
    double temp = dx / dy;
    sn = 1.0 / sqrt(1.0 + temp * temp);
    cs = temp * sn;
  }
  else
  {
    double temp = dy / dx;
    cs = 1.0 / sqrt(1.0 + temp * temp);
    sn = temp * cs;
  }
}


void NOX::FSI::LinearSystemGCR::ApplyPlaneRotation(double& dx, double& dy, double& cs, double& sn)
{
  double temp = cs * dx + sn * dy;
  dy = -sn * dx + cs * dy;
  dx = temp;
}


bool NOX::FSI::LinearSystemGCR::applyRightPreconditioning(bool useTranspose,
    Teuchos::ParameterList& params, const NOX::Epetra::Vector& input,
    NOX::Epetra::Vector& result) const
{
  if (&result != &input) result = input;
  return true;
}


Teuchos::RCP<NOX::Epetra::Scaling> NOX::FSI::LinearSystemGCR::getScaling() { return scaling; }


void NOX::FSI::LinearSystemGCR::resetScaling(
    const Teuchos::RCP<NOX::Epetra::Scaling>& scalingObject)
{
  scaling = scalingObject;
}


bool NOX::FSI::LinearSystemGCR::computeJacobian(const NOX::Epetra::Vector& x)
{
  bool success = jacInterfacePtr->computeJacobian(x.getEpetraVector(), *jacPtr);
  return success;
}


bool NOX::FSI::LinearSystemGCR::createPreconditioner(
    const NOX::Epetra::Vector& x, Teuchos::ParameterList& p, bool recomputeGraph) const
{
  return false;
}


bool NOX::FSI::LinearSystemGCR::destroyPreconditioner() const { return false; }


bool NOX::FSI::LinearSystemGCR::recomputePreconditioner(
    const NOX::Epetra::Vector& x, Teuchos::ParameterList& linearSolverParams) const
{
  return false;
}


NOX::FSI::LinearSystemGCR::PreconditionerReusePolicyType
NOX::FSI::LinearSystemGCR::getPreconditionerPolicy(bool advanceReuseCounter)
{
  return PRPT_REBUILD;
}


bool NOX::FSI::LinearSystemGCR::isPreconditionerConstructed() const { return false; }


bool NOX::FSI::LinearSystemGCR::hasPreconditioner() const { return false; }


Teuchos::RCP<const Epetra_Operator> NOX::FSI::LinearSystemGCR::getJacobianOperator() const
{
  return jacPtr;
}


Teuchos::RCP<Epetra_Operator> NOX::FSI::LinearSystemGCR::getJacobianOperator() { return jacPtr; }


Teuchos::RCP<const Epetra_Operator> NOX::FSI::LinearSystemGCR::getGeneratedPrecOperator() const
{
  return Teuchos::null;
}


Teuchos::RCP<Epetra_Operator> NOX::FSI::LinearSystemGCR::getGeneratedPrecOperator()
{
  return Teuchos::null;
}


void NOX::FSI::LinearSystemGCR::setJacobianOperatorForSolve(
    const Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
  jacPtr = Teuchos::rcp_const_cast<Epetra_Operator>(solveJacOp);
  jacType = getOperatorType(*solveJacOp);
}


void NOX::FSI::LinearSystemGCR::setPrecOperatorForSolve(
    const Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
  throwError("setPrecOperatorForSolve", "no preconditioner supported");
}


void NOX::FSI::LinearSystemGCR::throwError(
    const std::string& functionName, const std::string& errorMsg) const
{
  if (utils.isPrintType(NOX::Utils::Error))
  {
    utils.out() << "NOX::FSI::LinearSystemGCR::" << functionName << " - " << errorMsg << std::endl;
  }
  throw "NOX Error";
}


NOX::FSI::LinearSystemGCR::OperatorType NOX::FSI::LinearSystemGCR::getOperatorType(
    const Epetra_Operator& Op)
{
  //***************
  //*** NOTE: The order in which the following tests occur is important!
  //***************

  const Epetra_Operator* testOperator = 0;

  // Is it an Epetra_CrsMatrix ?
  testOperator = dynamic_cast<const Epetra_CrsMatrix*>(&Op);
  if (testOperator != 0) return EpetraCrsMatrix;

  // Is it an Epetra_VbrMatrix ?
  testOperator = dynamic_cast<const Epetra_VbrMatrix*>(&Op);
  if (testOperator != 0) return EpetraVbrMatrix;

  // Is it an Epetra_RowMatrix ?
  testOperator = dynamic_cast<const Epetra_RowMatrix*>(&Op);
  if (testOperator != 0) return EpetraRowMatrix;

  // Otherwise it must be an Epetra_Operator!
  return EpetraOperator;
}
