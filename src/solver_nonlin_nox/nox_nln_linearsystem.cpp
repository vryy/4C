/*-----------------------------------------------------------*/
/*!
\file nox_nln_linearsystem.cpp

\brief %NOX::NLN extension of the %NOX::Epetra::LinearSystem.

\maintainer Michael Hiermeier

\date Jul 2, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#include "nox_nln_linearsystem.H"
#include "nox_nln_interface_jacobian.H"
#include "nox_nln_interface_required.H"
#include "nox_nln_linearsystem_prepostoperator.H"
#include "nox_nln_solver_ptc.H"
#include "nox_nln_aux.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_serialdensevector.H"

#include "../drt_lib/epetra_utils.H"

#include <unordered_map>

#include <Epetra_Vector.h>
#include <Epetra_LinearProblem.h>

#include <Teuchos_ParameterList.hpp>

#include <NOX_Epetra_Scaling.H>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include "../drt_structure_new/str_nln_linearsystem_scaling.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& jacobian,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const Teuchos::RCP<LINALG::SparseOperator>& preconditioner,
    const NOX::Epetra::Vector& cloneVector, const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject)
    : utils_(printParams),
      solvers_(solvers),
      reqInterfacePtr_(iReq),
      jacInterfacePtr_(iJac),
      jacType_(NOX::NLN::LinSystem::LinalgSparseOperator),
      precInterfacePtr_(iPrec),
      precType_(NOX::NLN::LinSystem::LinalgSparseOperator),
      precPtr_(preconditioner),
      precMatrixSource_(SeparateMatrix),
      scaling_(scalingObject),
      conditionNumberEstimate_(0.0),
      timer_(cloneVector.getEpetraVector().Comm()),
      timeCreatePreconditioner_(0.0),
      timeApplyJacbianInverse_(0.0),
      resNorm2_(0.0),
      prePostOperatorPtr_(Teuchos::null),
      jacPtr_(jacobian)
{
  // Jacobian operator is supplied
  jacType_ = NOX::NLN::AUX::GetOperatorType(Jacobian());

  reset(linearSolverParams);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& jacobian,
    const Teuchos::RCP<NOX::Epetra::Interface::Preconditioner>& iPrec,
    const Teuchos::RCP<LINALG::SparseOperator>& preconditioner,
    const NOX::Epetra::Vector& cloneVector)
    : utils_(printParams),
      solvers_(solvers),
      reqInterfacePtr_(iReq),
      jacInterfacePtr_(iJac),
      jacType_(NOX::NLN::LinSystem::LinalgSparseOperator),
      precInterfacePtr_(iPrec),
      precType_(NOX::NLN::LinSystem::LinalgSparseOperator),
      precPtr_(preconditioner),
      precMatrixSource_(SeparateMatrix),
      scaling_(Teuchos::null),
      conditionNumberEstimate_(0.0),
      timer_(cloneVector.getEpetraVector().Comm()),
      timeCreatePreconditioner_(0.0),
      timeApplyJacbianInverse_(0.0),
      resNorm2_(0.0),
      prePostOperatorPtr_(Teuchos::null),
      jacPtr_(jacobian)
{
  // Jacobian operator is supplied
  jacType_ = NOX::NLN::AUX::GetOperatorType(Jacobian());

  reset(linearSolverParams);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& jacobian, const NOX::Epetra::Vector& cloneVector,
    const Teuchos::RCP<NOX::Epetra::Scaling> scalingObject)
    : utils_(printParams),
      solvers_(solvers),
      reqInterfacePtr_(iReq),
      jacInterfacePtr_(iJac),
      jacType_(NOX::NLN::LinSystem::LinalgSparseOperator),
      precInterfacePtr_(Teuchos::null),
      precType_(NOX::NLN::LinSystem::LinalgSparseOperator),
      precPtr_(Teuchos::null),
      precMatrixSource_(SeparateMatrix),
      scaling_(scalingObject),
      conditionNumberEstimate_(0.0),
      timer_(cloneVector.getEpetraVector().Comm()),
      timeCreatePreconditioner_(0.0),
      timeApplyJacbianInverse_(0.0),
      resNorm2_(0.0),
      prePostOperatorPtr_(Teuchos::null),
      jacPtr_(jacobian)
{
  // Jacobian operator is supplied
  jacType_ = NOX::NLN::AUX::GetOperatorType(Jacobian());

  reset(linearSolverParams);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::NLN::LinearSystem::LinearSystem(Teuchos::ParameterList& printParams,
    Teuchos::ParameterList& linearSolverParams,
    const std::map<NOX::NLN::SolutionType, Teuchos::RCP<LINALG::Solver>>& solvers,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& iReq,
    const Teuchos::RCP<NOX::Epetra::Interface::Jacobian>& iJac,
    const Teuchos::RCP<LINALG::SparseOperator>& jacobian, const NOX::Epetra::Vector& cloneVector)
    : utils_(printParams),
      solvers_(solvers),
      reqInterfacePtr_(iReq),
      jacInterfacePtr_(iJac),
      jacType_(NOX::NLN::LinSystem::LinalgSparseOperator),
      precInterfacePtr_(Teuchos::null),
      precType_(NOX::NLN::LinSystem::LinalgSparseOperator),
      precPtr_(Teuchos::null),
      precMatrixSource_(SeparateMatrix),
      scaling_(Teuchos::null),
      conditionNumberEstimate_(0.0),
      timer_(cloneVector.getEpetraVector().Comm()),
      timeCreatePreconditioner_(0.0),
      timeApplyJacbianInverse_(0.0),
      resNorm2_(0.0),
      prePostOperatorPtr_(Teuchos::null),
      jacPtr_(jacobian)
{
  // Jacobian operator is supplied
  jacType_ = NOX::NLN::AUX::GetOperatorType(Jacobian());

  reset(linearSolverParams);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::reset(Teuchos::ParameterList& p)
{
  zeroInitialGuess_ = p.get<bool>("Zero Initial Guess", false);

  manualScaling_ = p.get<bool>("Compute Scaling Manually", true);

  // Place linear solver details in the "Output" sublist of the
  // "Linear Solver" parameter list
  outputSolveDetails_ = p.get<bool>("Output Solver Details", true);

  // set the pre/post-operator
  resetPrePostOperator(p);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::resetPrePostOperator(Teuchos::ParameterList& p)
{
  if (prePostOperatorPtr_.is_null())
    prePostOperatorPtr_ = Teuchos::rcp(new NOX::NLN::LinSystem::PrePostOperator(p));
  else
    prePostOperatorPtr_->reset(p);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::applyJacobianBlock(const NOX::Epetra::Vector& input,
    Teuchos::RCP<NOX::Epetra::Vector>& result, unsigned rbid, unsigned cbid) const
{
  const LINALG::SparseMatrix& blockc = getJacobianBlock(rbid, cbid);
  LINALG::SparseMatrix& block = const_cast<LINALG::SparseMatrix&>(blockc);
  const Epetra_Map& domainmap = block.DomainMap();
  const Epetra_Map& rangemap = block.RangeMap();

  const Epetra_Vector& input_epetra = input.getEpetraVector();
  Teuchos::RCP<const Epetra_Vector> input_apply = Teuchos::null;

  if (not input_epetra.Map().SameAs(domainmap))
  {
    input_apply = LINALG::ExtractMyVector(input_epetra, domainmap);
  }
  else
  {
    input_apply = Teuchos::rcpFromRef(input_epetra);
  }

  Teuchos::RCP<Epetra_Vector> result_apply = Teuchos::rcp(new Epetra_Vector(rangemap, true));

  block.SetUseTranspose(false);
  int status = block.Apply(*input_apply, *result_apply);

  result = Teuchos::rcp(new NOX::Epetra::Vector(result_apply, NOX::Epetra::Vector::CreateView));

  return (status == 0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::applyJacobian(
    const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const
{
  Jacobian().SetUseTranspose(false);
  int status = Jacobian().Apply(input.getEpetraVector(), result.getEpetraVector());
  return (status == 0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::applyJacobianTranspose(
    const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result) const
{
  // Apply the Jacobian
  Jacobian().SetUseTranspose(true);
  int status = Jacobian().Apply(input.getEpetraVector(), result.getEpetraVector());
  Jacobian().SetUseTranspose(false);

  return (status == 0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::SetLinearProblemForSolve(Epetra_LinearProblem& linear_problem,
    LINALG::SparseOperator& jac, Epetra_Vector& lhs, Epetra_Vector& rhs) const
{
  linear_problem.SetOperator(jac.EpetraOperator().get());
  linear_problem.SetLHS(&lhs);
  linear_problem.SetRHS(&rhs);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::CompleteSolutionAfterSolve(
    const Epetra_LinearProblem& linProblem, Epetra_Vector& lhs) const
{
  /* nothing to do in the default case */
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::applyJacobianInverse(Teuchos::ParameterList& linearSolverParams,
    const NOX::Epetra::Vector& input, NOX::Epetra::Vector& result)
{
  /* Need non-const version of the input vector
   * Epetra_LinearProblem requires non-const versions so we can perform
   * scaling of the linear problem.
   * Same is valid for the prePostOperator. We want to have the
   * possibility to change the linear system. */
  NOX::Epetra::Vector& nonConstInput = const_cast<NOX::Epetra::Vector&>(input);

  prePostOperatorPtr_->runPreApplyJacobianInverse(nonConstInput, Jacobian(), *this);

  double startTime = timer_.WallTime();

  // calculate the residual norm
  resNorm2_ = nonConstInput.norm(NOX::Abstract::Vector::TwoNorm);

  // Zero out the delta X of the linear problem if requested by user.
  if (zeroInitialGuess_) result.init(0.0);

  // Create Epetra linear problem object for the linear solve
  /* Note: We switch from LINALG_objects to pure Epetra_objects.
   * This is necessary for the linear solver.
   *     LINALG::SparseMatrix ---> Epetra_CrsMatrix */
  Epetra_LinearProblem linProblem;
  SetLinearProblemForSolve(
      linProblem, Jacobian(), result.getEpetraVector(), nonConstInput.getEpetraVector());

  // ************* Begin linear system scaling *****************
  if (!Teuchos::is_null(scaling_))
  {
    if (!manualScaling_) scaling_->computeScaling(linProblem);

    scaling_->scaleLinearSystem(linProblem);

    if (utils_.isPrintType(NOX::Utils::Details)) utils_.out() << *scaling_ << std::endl;
  }
  // ************* End linear system scaling *******************

  // get current linear solver from the std_map
  Teuchos::RCP<LINALG::Solver> currSolver;
  NOX::NLN::SolutionType solType = GetActiveLinSolver(solvers_, currSolver);

  // set solver options if necessary
  SetSolverOptions(linearSolverParams, currSolver, solType);

  // solve
  int iter = linearSolverParams.get<int>("Number of Nonlinear Iterations", -10);
  if (iter == -10)
    throwError("applyJacobianInverse", "\"Number of Nonlinear Iterations\" was not specified");

  const int linsol_status = currSolver->NoxSolve(linProblem, true, iter == 0);
  if (linsol_status)
  {
    if (utils_.isPrintType(NOX::Utils::Warning))
      utils_.out() << "NOX::NLN::LinearSystem::applyJacobianInverse -- "
                      "linear solve failed (err = "
                   << linsol_status << ")\n";
  }

  // ************* Begin linear system unscaling *************
  if (!Teuchos::is_null(scaling_)) scaling_->unscaleLinearSystem(linProblem);
  // ************* End linear system unscaling ***************

  CompleteSolutionAfterSolve(linProblem, result.getEpetraVector());

  double endTime = timer_.WallTime();
  timeApplyJacbianInverse_ += (endTime - startTime);

  prePostOperatorPtr_->runPostApplyJacobianInverse(result, nonConstInput, Jacobian(), *this);

  return (linsol_status == 0);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::applyRightPreconditioning(bool useTranspose,
    Teuchos::ParameterList& params, const NOX::Epetra::Vector& input,
    NOX::Epetra::Vector& result) const
{
  if (&result != &input) result = input;
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::computeJacobian(const NOX::Epetra::Vector& x)
{
  prePostOperatorPtr_->runPreComputeJacobian(Jacobian(), x.getEpetraVector(), *this);

  bool success = jacInterfacePtr_->computeJacobian(x.getEpetraVector(), Jacobian());

  prePostOperatorPtr_->runPostComputeJacobian(Jacobian(), x.getEpetraVector(), *this);
  return success;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::computeFandJacobian(
    const NOX::Epetra::Vector& x, NOX::Epetra::Vector& rhs)
{
  prePostOperatorPtr_->runPreComputeFandJacobian(
      rhs.getEpetraVector(), Jacobian(), x.getEpetraVector(), *this);

  const bool success =
      Teuchos::rcp_dynamic_cast<NOX::NLN::Interface::Jacobian>(jacInterfacePtr_, true)
          ->computeFandJacobian(x.getEpetraVector(), rhs.getEpetraVector(), Jacobian());

  prePostOperatorPtr_->runPostComputeFandJacobian(
      rhs.getEpetraVector(), Jacobian(), x.getEpetraVector(), *this);
  return success;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::computeCorrectionSystem(const enum CorrectionType type,
    const NOX::Abstract::Group& grp, const NOX::Epetra::Vector& x, NOX::Epetra::Vector& rhs)
{
  prePostOperatorPtr_->runPreComputeFandJacobian(
      rhs.getEpetraVector(), Jacobian(), x.getEpetraVector(), *this);

  const bool success =
      Teuchos::rcp_dynamic_cast<NOX::NLN::Interface::Jacobian>(jacInterfacePtr_, true)
          ->computeCorrectionSystem(
              type, grp, x.getEpetraVector(), rhs.getEpetraVector(), Jacobian());

  prePostOperatorPtr_->runPostComputeFandJacobian(
      rhs.getEpetraVector(), Jacobian(), x.getEpetraVector(), *this);

  return success;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::Scaling> NOX::NLN::LinearSystem::getScaling() { return scaling_; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& scalingObject)
{
  scaling_ = scalingObject;
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::adjustPseudoTimeStep(double& delta, const double& stepSize,
    const NOX::Epetra::Vector& dir, const NOX::Epetra::Vector& rhs,
    const NOX::NLN::Solver::PseudoTransient& ptcsolver)
{
  const Epetra_Vector& scalingDiagOp = ptcsolver.getScalingDiagOperator();
  // ---------------------------------------------------------------------
  // first undo the modification of the jacobian
  // ---------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> v = Teuchos::rcp(new Epetra_Vector(scalingDiagOp));
  v->Scale(ptcsolver.getInversePseudoTimeStep());
  Teuchos::RCP<LINALG::SparseMatrix> jac =
      Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(JacobianPtr());
  if (jac.is_null()) throwError("adjustPseudoTimeStep()", "Cast to LINALG::SparseMatrix failed!");
  // get the diagonal terms of the jacobian
  Teuchos::RCP<Epetra_Vector> diag = LINALG::CreateVector(jac->RowMap(), false);
  jac->ExtractDiagonalCopy(*diag);
  diag->Update(-1.0, *v, 1.0);
  // Finally undo the changes
  jac->ReplaceDiagonalValues(*diag);

  // ---------------------------------------------------------------------
  // calculate the least squares approximated corrected pseudo time step
  // ---------------------------------------------------------------------
  /* evaluate the first vector:
   *    eta^{-1} F_{n-1} + (\nabla_{x} F_{n-1})^{T} d_{n-1}             */
  double stepSizeInv = 1.0 / stepSize;
  Teuchos::RCP<Epetra_Vector> vec_1 = LINALG::CreateVector(jac->RowMap(), true);
  Teuchos::RCP<Epetra_Vector> vec_2 = Teuchos::rcp(new Epetra_Vector(rhs.getEpetraVector()));
  jac->Multiply(false, dir.getEpetraVector(), *vec_1);
  vec_2->Scale(stepSizeInv);
  vec_1->Update(1.0, *vec_2, 1.0);
  /* evaluate the second vector:              d^{T} V                   */
  vec_2->Multiply(1.0, scalingDiagOp, dir.getEpetraVector(), 0.0);

  // finally evaluate the scalar product
  double numerator = 0.0;
  double denominator = 0.0;
  vec_2->Dot(*vec_1, &numerator);
  vec_1->Dot(*vec_1, &denominator);

  // ---------------------------------------------------------------------
  // show the error (L2-norm)
  // ---------------------------------------------------------------------
  Teuchos::RCP<Epetra_Vector> vec_err = LINALG::CreateVector(jac->RowMap(), true);
  vec_err->Update(delta, *vec_1, 1.0, *vec_2, 0.0);
  double error_start = 0.0;
  vec_err->Norm2(&error_start);

  delta = -numerator / denominator;

  // ---------------------------------------------------------------------
  // show the actual remaining error (L2-norm)
  // ---------------------------------------------------------------------
  vec_err->Update(delta, *vec_1, 1.0, *vec_2, 0.0);
  double error_end = 0.0;
  vec_err->Norm2(&error_end);
  if (utils_.isPrintType(NOX::Utils::Details))
  {
    utils_.out() << "| Error: " << std::setw(5) << std::setprecision(3) << std::scientific
                 << error_start << " -> " << std::setw(5) << std::setprecision(3) << std::scientific
                 << error_end << " |" << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::throwError(
    const std::string& functionName, const std::string& errorMsg) const
{
  if (utils_.isPrintType(NOX::Utils::Error))
  {
    utils_.out() << "NOX::NLN::LinearSystem::" << functionName << " - " << errorMsg << std::endl;
  }
  throw "NOX Error";
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const NOX::Epetra::Interface::Required> NOX::NLN::LinearSystem::getRequiredInterface()
    const
{
  return reqInterfacePtr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const NOX::Epetra::Interface::Jacobian> NOX::NLN::LinearSystem::getJacobianInterface()
    const
{
  return jacInterfacePtr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const NOX::Epetra::Interface::Preconditioner>
NOX::NLN::LinearSystem::getPrecInterface() const
{
  return precInterfacePtr_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Operator> NOX::NLN::LinearSystem::getJacobianOperator() const
{
  return JacobianPtr();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> NOX::NLN::LinearSystem::getJacobianOperator()
{
  return JacobianPtr();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const enum NOX::NLN::LinSystem::OperatorType& NOX::NLN::LinearSystem::getJacobianOperatorType()
    const
{
  return jacType_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::setJacobianOperatorForSolve(
    const Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
  const Teuchos::RCP<const LINALG::SparseOperator>& linalgSprOp =
      Teuchos::rcp_dynamic_cast<const LINALG::SparseOperator>(solveJacOp);
  if (linalgSprOp.is_null())
    throwError("setJacobianOperatorForSolve", "dynamic_cast to LINALG_SparseOperator failed!");

  SetJacobianOperatorForSolve(linalgSprOp);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::SetJacobianOperatorForSolve(
    const Teuchos::RCP<const LINALG::SparseOperator>& solveJacOp)
{
  if (jacType_ != NOX::NLN::AUX::GetOperatorType(*solveJacOp))
    throwError("SetJacobianOperatorForSolve", "wrong operator type!");

  jacPtr_ = Teuchos::rcp_const_cast<LINALG::SparseOperator>(solveJacOp);
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::createPreconditioner(
    const NOX::Epetra::Vector& x, Teuchos::ParameterList& p, bool recomputeGraph) const
{
  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::destroyPreconditioner() const { return false; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::recomputePreconditioner(
    const NOX::Epetra::Vector& x, Teuchos::ParameterList& linearSolverParams) const
{
  return false;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::Epetra::LinearSystem::PreconditionerReusePolicyType
NOX::NLN::LinearSystem::getPreconditionerPolicy(bool advanceReuseCounter)
{
  return PRPT_REBUILD;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::isPreconditionerConstructed() const { return false; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::hasPreconditioner() const { return false; }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Operator> NOX::NLN::LinearSystem::getGeneratedPrecOperator() const
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> NOX::NLN::LinearSystem::getGeneratedPrecOperator()
{
  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::setPrecOperatorForSolve(
    const Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
  throwError("setPrecOperatorForSolve", "no preconditioner supported");
  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::NLN::LinearSystem::DestroyJacobian()
{
  jacPtr_ = Teuchos::null;
  return true;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const LINALG::SparseMatrix& NOX::NLN::LinearSystem::getJacobianBlock(
    unsigned rbid, unsigned cbid) const
{
  switch (jacType_)
  {
    case LinSystem::LinalgBlockSparseMatrix:
    {
      typedef LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy> linalg_bsm;
      const linalg_bsm& jac_block = dynamic_cast<const linalg_bsm&>(Jacobian());

      if (rbid >= static_cast<unsigned>(jac_block.Rows()) or
          cbid >= static_cast<unsigned>(jac_block.Cols()))
        dserror(
            "The given row/column block ids exceed the block dimension of "
            "the jacobian matrix.");

      return jac_block.Matrix(rbid, cbid);
    }
    case LinSystem::LinalgSparseMatrix:
    {
      if (rbid != 0 or cbid != 0) dserror("There is only one block for a LINALG::SparseMatrix!");

      const LINALG::SparseMatrix& jac_sm = dynamic_cast<const LINALG::SparseMatrix&>(Jacobian());

      return jac_sm;
    }
    default:
    {
      dserror("Unsupported LinSystem::OperatorType: %d | %s", jacType_,
          NOX::NLN::LinSystem::OperatorType2String(jacType_).c_str());
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const Epetra_Map& NOX::NLN::LinearSystem::getJacobianRangeMap(unsigned rbid, unsigned cbid) const
{
  return getJacobianBlock(rbid, cbid).RangeMap();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> NOX::NLN::LinearSystem::getDiagonalOfJacobian(unsigned diag_bid) const
{
  const LINALG::SparseMatrix& diag_block = getJacobianBlock(diag_bid, diag_bid);
  const Epetra_Map& rmap = diag_block.RangeMap();

  Teuchos::RCP<Epetra_Vector> diag_copy = Teuchos::rcp(new Epetra_Vector(rmap, true));

  diag_block.ExtractDiagonalCopy(*diag_copy);

  return diag_copy;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::replaceDiagonalOfJacobian(
    const Epetra_Vector& new_diag, unsigned diag_bid)
{
  const LINALG::SparseMatrix& diag_block = getJacobianBlock(diag_bid, diag_bid);
  LINALG::SparseMatrix& mod_diag_block = const_cast<LINALG::SparseMatrix&>(diag_block);

  CATCH_EPETRA_ERROR(mod_diag_block.ReplaceDiagonalValues(new_diag));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double NOX::NLN::LinearSystem::computeSerialConditionNumberOfJacobian(
    const LinSystem::ConditionNumber condnum_type) const
{
  if (Jacobian().Comm().NumProc() > 1) dserror("Currently only one processor is supported!");

  LINALG::SerialDenseMatrix dense_jac;
  convertJacobianToDenseMatrix(dense_jac);

  Epetra_LAPACK lapack;

  const int N = dense_jac.N();
  const int M = dense_jac.M();
  double* A = dense_jac.A();
  const int LDA = dense_jac.LDA();
  double anorm = 0.0;
  char ctype = char();
  switch (condnum_type)
  {
    case LinSystem::ConditionNumber::one_norm:
      anorm = dense_jac.NormOne();
      ctype = '1';
      break;
    case LinSystem::ConditionNumber::inf_norm:
      anorm = dense_jac.InfNorm();
      ctype = 'I';
      break;
    default:
      dserror("Unsupported");
      exit(EXIT_FAILURE);
  }
  double rcond = 0.0;

  std::vector<double> work(4 * N);
  std::vector<int> iwork(N);
  int info = 0;

  // pre-step: compute the factorization of the matrix A (LU-decomposition)
  std::vector<int> ipiv(std::min(M, N));
  lapack.GETRF(M, N, A, LDA, ipiv.data(), &info);

  if (info) dserror("Error detected in LAPACK::GETRF");

  // compute the condition number
  lapack.GECON(ctype, N, A, LDA, anorm, &rcond, work.data(), iwork.data(), &info);

  if (info) dserror("Error detected in LAPACK::GECON");

  if (rcond < std::numeric_limits<double>::min())
    dserror(
        "The reciprocal condition number is zero. The jacobian"
        " seems to be singular.");

  return 1.0 / rcond;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::computeSerialEigenvaluesOfJacobian(
    LINALG::SerialDenseVector& reigenvalues, LINALG::SerialDenseVector& ieigenvalues) const
{
  if (Jacobian().Comm().NumProc() > 1) dserror("Currently only one processor is supported!");

  LINALG::SerialDenseMatrix dense_jac;
  convertJacobianToDenseMatrix(dense_jac);

  throwIfZeroRow(dense_jac);
  solveNonSymmEigenValueProblem(dense_jac, reigenvalues, ieigenvalues);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::prepareBlockDenseMatrix(
    const LINALG::BlockSparseMatrixBase& block_sparse, LINALG::SerialDenseMatrix& block_dense) const
{
  const int grows = block_sparse.FullRangeMap().NumGlobalElements();
  const int gcols = block_sparse.FullDomainMap().NumGlobalElements();

  block_dense.Reshape(grows, gcols);
  if (block_dense.N() != block_dense.M())
    dserror("The complete block dense matrix is not quadratic!");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::throwIfZeroRow(const LINALG::SerialDenseMatrix& block_dense) const
{
  for (int r = 0; r < block_dense.N(); ++r)
  {
    double csum = 0.0;
    for (int c = 0; c < block_dense.M(); ++c) csum += std::abs(block_dense(r, c));

    if (std::abs(csum) < std::numeric_limits<double>::epsilon()) dserror("Zero row detected!");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::convertJacobianToDenseMatrix(
    LINALG::SerialDenseMatrix& dense_jac) const
{
  switch (getJacobianOperatorType())
  {
    case NOX::NLN::LinSystem::LinalgBlockSparseMatrix:
    {
      typedef LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>
          linalg_blocksparsematrix;
      const linalg_blocksparsematrix& block_sparse =
          dynamic_cast<const linalg_blocksparsematrix&>(Jacobian());

      const Epetra_Map& full_rangemap = block_sparse.FullRangeMap();
      const Epetra_Map& full_domainmap = block_sparse.FullDomainMap();

      prepareBlockDenseMatrix(block_sparse, dense_jac);

      for (int r = 0; r < block_sparse.Rows(); ++r)
      {
        for (int c = 0; c < block_sparse.Cols(); ++c)
        {
          const LINALG::SparseMatrix& block = block_sparse.Matrix(r, c);
          convertSparseToDenseMatrix(block, dense_jac, full_rangemap, full_domainmap);
        }
      }

      break;
    }
    case NOX::NLN::LinSystem::LinalgSparseMatrix:
    {
      const LINALG::SparseMatrix& jac = dynamic_cast<const LINALG::SparseMatrix&>(Jacobian());

      convertSparseToDenseMatrix(jac, dense_jac, jac.RangeMap(), jac.DomainMap());

      break;
    }
    default:
    {
      dserror("Unsupported jacobian operator type!");
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::convertSparseToDenseMatrix(const LINALG::SparseMatrix& sparse,
    LINALG::SerialDenseMatrix& dense, const Epetra_Map& full_rangemap,
    const Epetra_Map& full_domainmap) const
{
  if (not sparse.Filled()) dserror("The sparse matrix must be filled!");
  Teuchos::RCP<const Epetra_CrsMatrix> crs_mat = sparse.EpetraMatrix();

  if (dense.N() == 0 or dense.M() == 0)
  {
    const int numgrows = sparse.RowMap().NumGlobalElements();
    const int numgcols = sparse.ColMap().NumGlobalElements();
    dense.Reshape(numgrows, numgcols);
    dense.Zero();
  }

  const int num_myrows = sparse.RowMap().NumMyElements();
  const int* rgids = sparse.RowMap().MyGlobalElements();

  for (int rlid = 0; rlid < num_myrows; ++rlid)
  {
    // number of entries in this row
    int numentries = 0;
    // values of the entries in this row
    double* rvals = NULL;
    // local indices
    int* indices = NULL;

    crs_mat->ExtractMyRowView(rlid, numentries, rvals, indices);

    const int rgid = rgids[rlid];
    const int full_rlid = full_rangemap.LID(rgid);
    if (full_rlid == -1) dserror("Row/Range: Couldn't find the corresponding LID to GID %d", rgid);

    double csum = 0.0;
    for (int i = 0; i < numentries; ++i)
    {
      const int cgid = sparse.ColMap().GID(indices[i]);
      if (cgid == -1)
        dserror("Column/Domain: Couldn't find the corresponding GID to LID %d", indices[i]);

      const int full_clid = full_domainmap.LID(cgid);
      if (full_clid == -1)
        dserror("Column/Domain: Couldn't find the corresponding LID to GID %d", cgid);

      dense(full_rlid, full_clid) = rvals[i];
      csum += std::abs(rvals[i]);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::solveNonSymmEigenValueProblem(LINALG::SerialDenseMatrix& mat,
    LINALG::SerialDenseVector& reigenvalues, LINALG::SerialDenseVector& ieigenvalues) const
{
  // start debugging
  //  std::ofstream of("/home/hiermeier/Workspace/o/test/mat.csv",std::ios_base::out);
  //  const int num_cols = mat.M();
  //  for ( int r=0; r<mat.N(); ++r )
  //    for ( int c=0; c<num_cols; ++c )
  //      of << std::setprecision(16) << std::scientific << mat(r,c)
  //         << ( c != num_cols-1 ? "," : "\n" );
  //  of.close();
  // end debugging

  reigenvalues.Size(mat.N());
  ieigenvalues.Size(mat.N());

  callGEEV(mat, reigenvalues, ieigenvalues);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::callGGEV(LINALG::SerialDenseMatrix& mat,
    LINALG::SerialDenseVector& reigenvalues, LINALG::SerialDenseVector& ieigenvalues) const
{
  int info = 0;
  int lwork = 8.0 * mat.N();
  std::vector<double> work(lwork);

  if (mat.N() != mat.M()) dserror("A N-by-N real non-symmetric matrix is expected by GGEV!");

  // create dummy B matrix
  static bool first_execution = true;
  static LINALG::SerialDenseMatrix bmat;
  if (first_execution or bmat.N() != mat.N())
  {
    bmat.Reshape(mat.N(), mat.N());
    for (int r = 0; r < bmat.N(); ++r) bmat(r, r) = 1.0;
  }

  LINALG::SerialDenseVector beta(reigenvalues);

  Epetra_LAPACK lapack;
  lapack.GGEV('N', 'N', mat.N(), mat.A(), mat.LDA(), bmat.A(), bmat.LDA(), reigenvalues.A(),
      ieigenvalues.A(), beta.A(), NULL, 1, NULL, 1, work.data(), lwork, &info);

  // postprocess the actual eigenvalues
  for (int i = 0; i < reigenvalues.N(); ++i)
  {
    if (beta(i) == 0.0)
      dserror(
          "The BETA factor is equal to zero! The eigenvalues can not be"
          " computed.");
    reigenvalues(i) /= beta(i);
    ieigenvalues(i) /= beta(i);
  }

  if (info) dserror("GGEV failed! (info = %d)", info);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::NLN::LinearSystem::callGEEV(LINALG::SerialDenseMatrix& mat,
    LINALG::SerialDenseVector& reigenvalues, LINALG::SerialDenseVector& ieigenvalues) const
{
  if (mat.N() != mat.M()) dserror("A N-by-N real non-symmetric matrix is expected by GEEV!");

  int info = 0;
  int lwork = 3.0 * mat.N();
  std::vector<double> work(lwork);

  Epetra_LAPACK lapack;
  lapack.GEEV('N', 'N', mat.N(), mat.A(), mat.LDA(), reigenvalues.A(), ieigenvalues.A(), NULL, 1,
      NULL, 1, work.data(), lwork, &info);

  if (info) dserror("GEEV failed!");
}
