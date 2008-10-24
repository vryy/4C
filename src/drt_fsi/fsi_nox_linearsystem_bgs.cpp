#ifdef CCADISCRET

#include <vector>

#include <Epetra_CrsMatrix.h>
#include <Epetra_LinearProblem.h>
#include <Epetra_Operator.h>
#include <Epetra_RowMatrix.h>
#include <Epetra_SerialDenseMatrix.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_VbrMatrix.h>
#include <Epetra_Vector.h>

#include <blitz/array.h>

#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include "../drt_lib/linalg_solver.H"

#include "fsi_nox_aitken.H"
#include "fsi_nox_fixpoint.H"
#include "fsi_nox_jacobian.H"
#include "fsi_nox_linearsystem_bgs.H"
#include "fsi_nox_linearsystem_gcr.H"
#include "fsi_nox_mpe.H"
#include "fsi_nox_sd.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::LinearBGS::LinearBGS(const LINALG::BlockSparseMatrixBase& A,
                               const Epetra_Vector &y,
                               Teuchos::RCP<LINALG::Solver> structure_solver,
                               Teuchos::RCP<LINALG::Solver> fluid_solver,
                               Teuchos::RCP<LINALG::Solver> ale_solver)
  :
  A_(A),
  y_(y),
  callcount_(0),
  structureSolver_(structure_solver),
  fluidSolver_(fluid_solver),
  aleSolver_(ale_solver)
{
  sy_ = A_.RangeExtractor().ExtractVector(y_,0);
  fy_ = A_.RangeExtractor().ExtractVector(y_,1);
  ay_ = A_.RangeExtractor().ExtractVector(y_,2);

  tmpsx_ = Teuchos::rcp(new Epetra_Vector(A_.RangeMap(0)));
  tmpfx_ = Teuchos::rcp(new Epetra_Vector(A_.RangeMap(1)));
  tmpax_ = Teuchos::rcp(new Epetra_Vector(A_.RangeMap(2)));

  tmpsy_ = Teuchos::rcp(new Epetra_Vector(A_.RangeMap(0)));
  tmpfy_ = Teuchos::rcp(new Epetra_Vector(A_.RangeMap(1)));
  tmpay_ = Teuchos::rcp(new Epetra_Vector(A_.RangeMap(2)));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearBGS::computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  const LINALG::SparseMatrix& A00 = A_.Matrix(0,0);	// structureInnerOp
  const LINALG::SparseMatrix& A01 = A_.Matrix(0,1);
  const LINALG::SparseMatrix& A02 = A_.Matrix(0,2);
  const LINALG::SparseMatrix& A10 = A_.Matrix(1,0);
  const LINALG::SparseMatrix& A11 = A_.Matrix(1,1);	// fluidInnerOp
  const LINALG::SparseMatrix& A12 = A_.Matrix(1,2);
  const LINALG::SparseMatrix& A20 = A_.Matrix(2,0);
  const LINALG::SparseMatrix& A21 = A_.Matrix(2,1);
  const LINALG::SparseMatrix& A22 = A_.Matrix(2,2);	// aleInnerOp

  Teuchos::RCP<Epetra_Vector> sx = A_.DomainExtractor().ExtractVector(x,0);
  Teuchos::RCP<Epetra_Vector> fx = A_.DomainExtractor().ExtractVector(x,1);
  Teuchos::RCP<Epetra_Vector> ax = A_.DomainExtractor().ExtractVector(x,2);

  // Structure
  // sx_{k+1} = A00^(-1) (sy - A01 * fx_{k} - A02 * ax_{k})
  A01.Multiply(false,*fx,*tmpsx_);
  tmpsy_->Update(-1.0,*tmpsx_,1.0,*sy_,0.0);
  A02.Multiply(false,*ax,*tmpsx_);
  tmpsy_->Update(-1.0,*tmpsx_,1.0);
  structureSolver_->Solve(A00.EpetraMatrix(),sx,tmpsy_,true,callcount_==0);

  // ALE
  // ax_{k+1} = A22^(-1) (ay - A20 * sx_{k+1} - A21 * fx_{k+1})
  A20.Multiply(false,*sx,*tmpax_);
  tmpay_->Update(-1.0,*tmpax_,1.0,*ay_,0.0);
  A21.Multiply(false,*fx,*tmpax_);
  tmpay_->Update(-1.0,*tmpax_,1.0);
  aleSolver_->Solve(A22.EpetraMatrix(),ax,tmpay_,true,callcount_==0);

  // Fluid
  // fx_{k+1} = A11^(-1) (fy - A10 * sx_{k+1} - A12 * ax_{k})
  A10.Multiply(false,*sx,*tmpfx_);
  tmpfy_->Update(-1.0,*tmpfx_,1.0,*fy_,0.0);
  A12.Multiply(false,*ax,*tmpfx_);
  tmpfy_->Update(-1.0,*tmpfx_,1.0);
  fluidSolver_->Solve(A11.EpetraMatrix(),fx,tmpfy_,true,callcount_==0);

  A_.DomainExtractor().InsertVector(*sx,0,F);
  A_.DomainExtractor().InsertVector(*fx,1,F);
  A_.DomainExtractor().InsertVector(*ax,2,F);

  F.Update(-1.,x,1.);
  callcount_ += 1;
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::LinearBGSSolver::LinearBGSSolver(Teuchos::ParameterList& printParams,
                                           Teuchos::ParameterList& linearSolverParams,
                                           const Teuchos::RefCountPtr< NOX::Epetra::Interface::Jacobian>& iJac,
                                           const Teuchos::RefCountPtr<Epetra_Operator>& J,
                                           const NOX::Epetra::Vector& cloneVector,
                                           Teuchos::RCP < LINALG::Solver > structure_solver,
                                           Teuchos::RCP < LINALG::Solver > fluid_solver,
                                           Teuchos::RCP < LINALG::Solver > ale_solver,
                                           INPUTPARAMS::FSILinearBlockSolver linearsolverstrategy,
                                           const Teuchos::RefCountPtr< NOX::Epetra::Scaling> s)
  :
  utils_(printParams),
  jacInterfacePtr_(iJac),
  jacType_(EpetraOperator),
  precType_(EpetraOperator),
  jacPtr_(J),
  scaling_(s),
  conditionNumberEstimate_(0.0),
  structureSolver_(structure_solver),
  fluidSolver_(fluid_solver),
  aleSolver_(ale_solver),
  timer_(cloneVector.getEpetraVector().Comm()),
  timeApplyJacbianInverse_(0.0),
  linearsolverstrategy_(linearsolverstrategy)
{
  tmpVectorPtr_ = Teuchos::rcp(new NOX::Epetra::Vector(cloneVector));

  //cout << "STRUCTURE SOLVER: " << *structureSolver_ << " " << structureSolver_ << endl;
  //cout << "FLUID SOLVER: " << *fluidSolver_ << " " << fluidSolver_ << endl;
  //cout << "ALE SOLVER: " << *aleSolver_  << " " << aleSolver_ << endl;

  // Jacobian operator is supplied.
  // get type of it
  jacType_ = getOperatorType(*jacPtr_);

  reset(linearSolverParams);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::LinearBGSSolver::~LinearBGSSolver()
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::LinearBGSSolver::OperatorType NOX::FSI::LinearBGSSolver::getOperatorType(
  const Epetra_Operator& Op)
{
  // check per dynamik cast, welcher typ von jacobian operator �bergeben wurde

  const Epetra_Operator* testOperator = 0;

  testOperator
    = dynamic_cast<const LINALG::BlockSparseMatrix<LINALG::DefaultBlockMatrixStrategy>*>(&Op);
  if (testOperator != 0)
    return BlockSparseMatrix;

  testOperator = dynamic_cast<const LINALG::SparseMatrix*>(&Op);
  if (testOperator != 0)
    return SparseMatrix;

  testOperator = dynamic_cast<const Epetra_CrsMatrix*>(&Op);
  if (testOperator != 0)
    return EpetraCrsMatrix;

  testOperator = dynamic_cast<const Epetra_VbrMatrix*>(&Op);
  if (testOperator != 0)
    return EpetraVbrMatrix;

  testOperator = dynamic_cast<const Epetra_RowMatrix*>(&Op);
  if (testOperator != 0)
    return EpetraRowMatrix;

  return EpetraOperator;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearBGSSolver::reset(Teuchos::ParameterList& linearSolverParams)
{
  zeroInitialGuess_ = linearSolverParams.get("Zero Initial Guess", false);
  manualScaling_ = linearSolverParams.get("Compute Scaling Manually", true);
  outputSolveDetails_ = linearSolverParams.get("Output Solver Details", true);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearBGSSolver::applyJacobian(const NOX::Epetra::Vector& input,
                                              NOX::Epetra::Vector& result) const
{
  jacPtr_->SetUseTranspose(false);
  int status = jacPtr_->Apply(input.getEpetraVector(), result.getEpetraVector());
  return status == 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearBGSSolver::applyJacobianTranspose(const NOX::Epetra::Vector& input,
                                                       NOX::Epetra::Vector& result) const
{
  jacPtr_->SetUseTranspose(true);
  int status = jacPtr_->Apply(input.getEpetraVector(), result.getEpetraVector());
  jacPtr_->SetUseTranspose(false);
  return status == 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearBGSSolver::applyJacobianInverse(Teuchos::ParameterList &p,
                                                     const NOX::Epetra::Vector &input,
                                                     NOX::Epetra::Vector &result)
{
  double startTime = timer_.WallTime();

  // Zero out the delta X of the linear problem if requested by user.
  if (zeroInitialGuess_)
    result.init(0.0);

  int maxit = p.get("Max Iterations", 30);
  double tol = p.get("Tolerance", 1.0e-10);

  bgs(Teuchos::dyn_cast<const LINALG::BlockSparseMatrixBase>(*jacPtr_.get()),
      result,
      input,
      maxit,
      tol);

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
void NOX::FSI::LinearBGSSolver::bgs(const LINALG::BlockSparseMatrixBase& A,
                                    NOX::Epetra::Vector& result,
                                    const NOX::Epetra::Vector& input,
                                    int& maxit,
                                    double& tol)
{
  // solve linear BGS system via NOX

  Teuchos::ParameterList nlParams;
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& lineSearchParams = nlParams.sublist("Line Search");

  nlParams.set("Norm abs F", tol);
  nlParams.set("Max Iterations", maxit);

  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  ///////////////////////////////////////////////////////////////////

  // fixed relaxation parameter
  //lineSearchParams.set("Method", "Full Step");
  //lineSearchParams.sublist("Full Step").set("Full Step", fsidyn.get<double>("RELAX"));

  switch (linearsolverstrategy_)
  {
  case INPUTPARAMS::fsi_BGSAitken:
  {
    dirParams.set("Method","User Defined");
    Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
      Teuchos::rcp(new NOX::FSI::FixPointFactory());
    dirParams.set("User Defined Direction Factory",fixpointfactory);

    nlParams.set("Jacobian", "None");

    // Aitken relaxation
    Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> aitkenfactory =
      Teuchos::rcp(new NOX::FSI::AitkenFactory());
    lineSearchParams.set("Method","User Defined");
    lineSearchParams.set("User Defined Line Search Factory", aitkenfactory);

    lineSearchParams.sublist("Aitken").set("max step size", fsidyn.get<double>("MAXOMEGA"));
    break;
  }
  case INPUTPARAMS::fsi_BGSVectorExtrapolation:
  {
    nlParams.set("Jacobian", "None");
    dirParams.set("Method","User Defined");

    Teuchos::RCP<NOX::Direction::UserDefinedFactory> factory =
      Teuchos::rcp(new NOX::FSI::MinimalPolynomialFactory());
    dirParams.set("User Defined Direction Factory",factory);

    Teuchos::ParameterList& exParams = dirParams.sublist("Extrapolation");
    exParams.set("Tolerance", fsidyn.get<double>("BASETOL"));
    exParams.set("omega", fsidyn.get<double>("RELAX"));
    exParams.set("kmax", 10);
    exParams.set("Method", "RRE");

    lineSearchParams.set("Method", "Full Step");
    lineSearchParams.sublist("Full Step").set("Full Step", 1.0);
    break;
  }
  case INPUTPARAMS::fsi_BGSJacobianFreeNewtonKrylov:
    // Newton-Krylov on block Gauss-Seidel? Can there be any difference to
    // preconditioned block Newton-Krylov (real monolithic)?
  default:
    dserror("unsupported linear block solver strategy: %d", linearsolverstrategy_);
  }

  ///////////////////////////////////////////////////////////////////

  Teuchos::RCP<NOX::Epetra::Interface::Required> bgs =
    Teuchos::rcp(new LinearBGS(A,input.getEpetraVector(),structureSolver_,fluidSolver_,aleSolver_));

  // Create the linear system
  // actually, we do not need any
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys = Teuchos::null;

  // Create the Group
  Teuchos::RCP<NOX::Epetra::Group> grp =
    Teuchos::rcp(new NOX::Epetra::Group(printParams, bgs, result, linSys));

  // Convergence Tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo = CreateStatusTest(nlParams, grp);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp,combo,RCP<ParameterList>(&nlParams,false));

  // solve the whole thing
  NOX::StatusTest::StatusType status = solver->solve();

  if (status != NOX::StatusTest::Converged)
    utils_.out() << RED "linear BGS solver failed to converge!" END_COLOR << "\n";

  // Get the Epetra_Vector with the final solution from the solver
  const NOX::Epetra::Group& finalGroup = Teuchos::dyn_cast<const NOX::Epetra::Group>(solver->getSolutionGroup());
  result = Teuchos::dyn_cast<const NOX::Epetra::Vector>(finalGroup.getX());

  maxit = nlParams.sublist("Output").get("Nonlinear Iterations",0);
  tol = nlParams.sublist("Output").get("2-Norm of Residual", 0.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo>
NOX::FSI::LinearBGSSolver::CreateStatusTest(ParameterList& nlParams,
                                            Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo       = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::StatusTest::Combo> converged   = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv    = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // setup the real tests
  CreateStatusTest(nlParams,grp,converged);

  return combo;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void
NOX::FSI::LinearBGSSolver::CreateStatusTest(ParameterList& nlParams,
                                            Teuchos::RCP<NOX::Epetra::Group> grp,
                                            Teuchos::RCP<NOX::StatusTest::Combo> converged)
{
  if (nlParams.isParameter("Norm abs F"))
  {
    Teuchos::RCP<NOX::StatusTest::NormF> absresid =
      Teuchos::rcp(new NOX::StatusTest::NormF(nlParams.get("Norm abs F", 1.0e-6)));
    converged->addStatusTest(absresid);
  }

  if (nlParams.isParameter("Norm Update"))
  {
    Teuchos::RCP<NOX::StatusTest::NormUpdate> update =
      Teuchos::rcp(new NOX::StatusTest::NormUpdate(nlParams.get("Norm Update", 1.0e-5)));
    converged->addStatusTest(update);
  }

  if (nlParams.isParameter("Norm rel F"))
  {
    Teuchos::RCP<NOX::StatusTest::NormF> relresid =
      Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), nlParams.get("Norm rel F", 1.0e-2)));
    converged->addStatusTest(relresid);
  }

  //Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms     = Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  //converged->addStatusTest(wrms);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearBGSSolver::applyRightPreconditioning(bool useTranspose,
                                                          Teuchos::ParameterList& params,
                                                          const NOX::Epetra::Vector& input,
                                                          NOX::Epetra::Vector& result) const
{
  if (&result != &input)
    result = input;
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP< NOX::Epetra::Scaling> NOX::FSI::LinearBGSSolver::getScaling()
{
  return scaling_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearBGSSolver::resetScaling(const Teuchos::RCP< NOX::Epetra::Scaling>& scalingObject)
{
  scaling_ = scalingObject;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearBGSSolver::computeJacobian(const NOX::Epetra::Vector& x)
{
  bool success = jacInterfacePtr_->computeJacobian(x.getEpetraVector(),
                                                   *jacPtr_);
  return success;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearBGSSolver::createPreconditioner(const NOX::Epetra::Vector& x,
                                                     Teuchos::ParameterList& p,
                                                     bool recomputeGraph) const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearBGSSolver::destroyPreconditioner() const

{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearBGSSolver::recomputePreconditioner(const NOX::Epetra::Vector& x,
                                                        Teuchos::ParameterList& linearSolverParams) const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::LinearBGSSolver::PreconditionerReusePolicyType
NOX::FSI::LinearBGSSolver::getPreconditionerPolicy(bool advanceReuseCounter)
{
  return PRPT_REBUILD;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearBGSSolver::isPreconditionerConstructed() const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearBGSSolver::hasPreconditioner() const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Operator>
NOX::FSI::LinearBGSSolver::getJacobianOperator() const
{
  return jacPtr_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> NOX::FSI::LinearBGSSolver::getJacobianOperator()
{
  return jacPtr_;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Operator>
NOX::FSI::LinearBGSSolver::getGeneratedPrecOperator() const
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator> NOX::FSI::LinearBGSSolver::getGeneratedPrecOperator()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearBGSSolver::setJacobianOperatorForSolve(
  const Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
  jacPtr_ = Teuchos::rcp_const_cast<Epetra_Operator>(solveJacOp);
  jacType_ = getOperatorType(*solveJacOp);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearBGSSolver::setPrecOperatorForSolve(
  const Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
  throwError("setPrecOperatorForSolve", "no preconditioner supported");
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearBGSSolver::throwError(const string& functionName,
                                           const string& errorMsg) const
{
  if (utils_.isPrintType(NOX::Utils::Error))

  {
    utils_.out() << "NOX::FSI::LinearBGSSolver::" << functionName << " - "
                 << errorMsg << endl;
  }
  throw "NOX Error";
}

#endif
