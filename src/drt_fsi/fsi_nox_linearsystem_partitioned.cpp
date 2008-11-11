#ifdef CCADISCRET

#include "fsi_nox_linearsystem_partitioned.H"
#include "fsi_partitionedmonolithic.H"

#include "fsi_nox_aitken.H"
#include "fsi_nox_fixpoint.H"
#include "fsi_nox_jacobian.H"
#include "fsi_nox_mpe.H"


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::LinearPartitioned::LinearPartitioned(::FSI::PartitionedMonolithic& algorithm,
                                               ADAPTER::Structure& structurefield,
                                               ADAPTER::Fluid& fluidfield,
                                               ADAPTER::Ale& alefield)
  : callcount_(0),
    algorithm_(algorithm),
    structurefield_(structurefield),
    fluidfield_(fluidfield),
    alefield_(alefield)
{
  srhs_ = structurefield_.RHS();
  frhs_ = fluidfield_    .RHS();
  arhs_ = alefield_      .RHS();

  s_ = structurefield_.SystemMatrix();
  f_ = fluidfield_    .SystemMatrix();
  a_ = alefield_      .SystemMatrix();

  // extract interface Dirichlet lines from fluid matrix and remove them
  // afterwards
  fluiddirichlet_ = f_->ExtractDirichletLines(*fluidfield_.Interface().CondMap());
  f_->ApplyDirichlet(*fluidfield_.Interface().CondMap(),true);

  // get an idea of interface displacement
  idispn_ = structurefield_.ExtractInterfaceDispn();
  iveln_ = algorithm_.FluidToStruct(fluidfield_.ExtractInterfaceVeln());

  sx_ = Teuchos::rcp(new Epetra_Vector(s_->RowMap()));
  fx_ = Teuchos::rcp(new Epetra_Vector(f_->RowMap()));
  ax_ = Teuchos::rcp(new Epetra_Vector(a_->RowMap()));

  stmp_ = Teuchos::rcp(new Epetra_Vector(s_->RowMap()));
  ftmp_ = Teuchos::rcp(new Epetra_Vector(f_->RowMap()));
  atmp_ = Teuchos::rcp(new Epetra_Vector(a_->RowMap()));

  slin_ = structurefield_.LinearSolver();
  flin_ = fluidfield_    .LinearSolver();
  alin_ = alefield_      .LinearSolver();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearPartitioned::computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  // We do it all ourselves here. In a better world we would just call linear
  // solve methods on our fields.

  Teuchos::RCP<const Epetra_Vector> idisp = Teuchos::rcp(&x,false);

  Teuchos::RCP<Epetra_Vector> iforce  = FluidOp(idisp, fillFlag);
  Teuchos::RCP<Epetra_Vector> idispnp = StructOp(iforce, fillFlag);

  F.Update(1.0,*idispnp,-1.0,*idisp,0.0);

  callcount_ += 1;
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
NOX::FSI::LinearPartitioned::StructOp(Teuchos::RCP<const Epetra_Vector> iforce, const FillType fillFlag)
{
  // prepare RHS

  if (fillFlag==User)
  {
    stmp_->PutScalar(0.);
  }
  else
  {
    stmp_->Update(1.0,*srhs_,0.0);
  }

  structurefield_.Interface().AddCondVector(iforce,stmp_);

  // solve

  slin_->Solve(s_->EpetraOperator(),sx_,stmp_,true,callcount_==0);

  // extract interface displacement

  Teuchos::RCP<Epetra_Vector> idisp =
    structurefield_.Interface().ExtractCondVector(sx_);

  return idisp;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
NOX::FSI::LinearPartitioned::FluidOp(Teuchos::RCP<const Epetra_Vector> idisp, const FillType fillFlag)
{
  //Teuchos::RCP<Epetra_Vector> ivel = algorithm_.StructToFluid(InterfaceVelocity(idisp));
  Teuchos::RCP<Epetra_Vector> ivel = algorithm_.StructToFluid(idisp);
  ivel->Scale(fluidfield_.TimeScaling());

  // prepare RHS

  if (fillFlag==User)
  {
    ftmp_->PutScalar(0.);
  }
  else
  {
    ftmp_->Update(1.0,*frhs_,0.0);
  }

  LINALG::ApplyDirichlettoSystem(fx_,ftmp_,ivel,*fluidfield_.Interface().CondMap());

  // solve

  flin_->Solve(f_->EpetraOperator(),fx_,ftmp_,true,callcount_==0);

  // calculate interface forces

  Teuchos::RCP<Epetra_Vector> force = Teuchos::rcp(new Epetra_Vector(fx_->Map()));
  fluiddirichlet_->Apply(*fx_,*force);

  // Can we have fluid forces at the coupling interface?
  force->Update(-1.0,*frhs_,1.0);

  Teuchos::RCP<Epetra_Vector> iforce =
    fluidfield_.Interface().ExtractCondVector(force);

  // make the interface forces a true force
  iforce->Scale(-fluidfield_.ResidualScaling());

  return algorithm_.FluidToStruct(iforce);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void
NOX::FSI::LinearPartitioned::AleOp(Teuchos::RCP<const Epetra_Vector> idisp)
{
  atmp_->Update(1.0,*arhs_,0.0);

  // set just those Dirichlet conditions at the interface
  // The matrix has already been modified.
  LINALG::ApplyDirichlettoSystem(ax_,atmp_,algorithm_.StructToAle(idisp),*alefield_.Interface().CondMap());

  alin_->Solve(a_->EpetraOperator(),ax_,atmp_,true,true);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
NOX::FSI::LinearPartitioned::InterfaceVelocity(Teuchos::RCP<const Epetra_Vector> idispnp) const
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
  Teuchos::RCP<Epetra_Vector> ivel = Teuchos::null;
  double dt = algorithm_.Dt();

  if (Teuchos::getIntegralValue<int>(fsidyn,"SECONDORDER") == 1)
  {
    ivel = rcp(new Epetra_Vector(*iveln_));
    ivel->Update(2./dt, *idispnp, -2./dt, *idispn_, -1.);
  }
  else
  {
    ivel = rcp(new Epetra_Vector(*idispn_));
    ivel->Update(1./dt, *idispnp, -1./dt);
  }
  return ivel;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearPartitioned::ExtractResult(Teuchos::RCP<const Epetra_Vector> idisp,
                                                Epetra_Vector& result)
{
  AleOp(idisp);
  algorithm_.SetupVector(result,sx_,fx_,ax_,0.);
  result.Scale(-1.0);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::LinearPartitionedSolver::LinearPartitionedSolver(Teuchos::ParameterList& printParams,
                                                           Teuchos::ParameterList& linearSolverParams,
                                                           const LINALG::MultiMapExtractor& extractor,
                                                           ::FSI::PartitionedMonolithic& algorithm,
                                                           ADAPTER::Structure& structurefield,
                                                           ADAPTER::Fluid& fluidfield,
                                                           ADAPTER::Ale& alefield,
                                                           INPUTPARAMS::FSILinearBlockSolver linearsolverstrategy)
  : utils_(printParams),
    extractor_(extractor),
    algorithm_(algorithm),
    structurefield_(structurefield),
    fluidfield_(fluidfield),
    alefield_(alefield),
    linearsolverstrategy_(linearsolverstrategy)
{
  reset(linearSolverParams);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::LinearPartitionedSolver::~LinearPartitionedSolver()
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearPartitionedSolver::reset(Teuchos::ParameterList& linearSolverParams)
{
  zeroInitialGuess_ = linearSolverParams.get("Zero Initial Guess", false);
  outputSolveDetails_ = linearSolverParams.get("Output Solver Details", true);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearPartitionedSolver::applyJacobian(const NOX::Epetra::Vector& input,
                                                      NOX::Epetra::Vector& result) const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearPartitionedSolver::applyJacobianTranspose(const NOX::Epetra::Vector& input,
                                                               NOX::Epetra::Vector& result) const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearPartitionedSolver::applyJacobianInverse(Teuchos::ParameterList &params,
                                                             const NOX::Epetra::Vector &input,
                                                             NOX::Epetra::Vector &result)
{
  TEUCHOS_FUNC_TIME_MONITOR("NOX::FSI::LinearPartitionedSolver::applyJacobianInverse");

  // Zero out the delta X of the linear problem if requested by user.
  if (zeroInitialGuess_)
    result.init(0.0);

  int maxit = params.get("Max Iterations", 100);
  double tol = params.get("Tolerance", 1.0e-10);

  LinearPartitionedSolve(result,input,maxit,tol);

  // Set the output parameters in the "Output" sublist
  if (outputSolveDetails_)
  {
    Teuchos::ParameterList& outputList = params.sublist("Output");
    int prevLinIters = outputList.get("Total Number of Linear Iterations", 0);
    int curLinIters = maxit;
    double achievedTol = tol;

    outputList.set("Number of Linear Iterations", curLinIters);
    outputList.set("Total Number of Linear Iterations", (prevLinIters + curLinIters));
    outputList.set("Achieved Tolerance", achievedTol);
  }

  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearPartitionedSolver::applyRightPreconditioning(bool useTranspose,
                                                                  Teuchos::ParameterList& params,
                                                                  const NOX::Epetra::Vector& input,
                                                                  NOX::Epetra::Vector& result) const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::Scaling> NOX::FSI::LinearPartitionedSolver::getScaling()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearPartitionedSolver::resetScaling(const Teuchos::RCP<NOX::Epetra::Scaling>& s)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearPartitionedSolver::computeJacobian(const NOX::Epetra::Vector& x)
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearPartitionedSolver::createPreconditioner(const NOX::Epetra::Vector& x,
                                                             Teuchos::ParameterList& p,
                                                             bool recomputeGraph) const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearPartitionedSolver::destroyPreconditioner() const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearPartitionedSolver::recomputePreconditioner(const NOX::Epetra::Vector& x,
                                                                Teuchos::ParameterList& linearSolverParams) const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::LinearPartitionedSolver::PreconditionerReusePolicyType
NOX::FSI::LinearPartitionedSolver::getPreconditionerPolicy(bool advanceReuseCounter)
{
  return PRPT_REBUILD;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearPartitionedSolver::isPreconditionerConstructed() const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool NOX::FSI::LinearPartitionedSolver::hasPreconditioner() const
{
  return false;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Operator>
NOX::FSI::LinearPartitionedSolver::getJacobianOperator() const
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator>
NOX::FSI::LinearPartitionedSolver::getJacobianOperator()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Operator>
NOX::FSI::LinearPartitionedSolver::getGeneratedPrecOperator() const
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Operator>
NOX::FSI::LinearPartitionedSolver::getGeneratedPrecOperator()
{
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearPartitionedSolver::setJacobianOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solveJacOp)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearPartitionedSolver::setPrecOperatorForSolve(const Teuchos::RCP<const Epetra_Operator>& solvePrecOp)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::LinearPartitionedSolver::LinearPartitionedSolve(NOX::Epetra::Vector& result,
                                                               const NOX::Epetra::Vector& input,
                                                               int& maxit,
                                                               double& tol)
{
  // solve linear partitioned system via NOX

  Teuchos::ParameterList nlParams;
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& lineSearchParams = nlParams.sublist("Line Search");

  nlParams.set("Norm rel F", tol);
  nlParams.set("Norm abs F", tol);
  nlParams.set("Max Iterations", maxit);

  Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
  solverOptions.set<std::string>("Status Test Check Type","Complete");

  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

  ///////////////////////////////////////////////////////////////////

  // fixed relaxation parameter
  //lineSearchParams.set("Method", "Full Step");
  //lineSearchParams.sublist("Full Step").set("Full Step", fsidyn.get<double>("RELAX"));

  switch (linearsolverstrategy_)
  {
  case INPUTPARAMS::fsi_PartitionedAitken:
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
  case INPUTPARAMS::fsi_PartitionedVectorExtrapolation:
  {
    nlParams.set("Jacobian", "None");
    dirParams.set("Method","User Defined");

    Teuchos::RCP<NOX::Direction::UserDefinedFactory> factory =
      Teuchos::rcp(new NOX::FSI::MinimalPolynomialFactory());
    dirParams.set("User Defined Direction Factory",factory);

    Teuchos::ParameterList& exParams = dirParams.sublist("Extrapolation");
    //exParams.set("Tolerance", fsidyn.get<double>("BASETOL"));
    exParams.set("Tolerance", tol);
    exParams.set("omega", fsidyn.get<double>("RELAX"));
    exParams.set("kmax", 25);
    exParams.set("Method", "RRE");

    lineSearchParams.set("Method", "Full Step");
    lineSearchParams.sublist("Full Step").set("Full Step", 1.0);
    break;
  }
  case INPUTPARAMS::fsi_PartitionedJacobianFreeNewtonKrylov:
    dserror("todo!");
    break;
  default:
    dserror("unsupported linear block solver strategy: %d", linearsolverstrategy_);
  }

  ///////////////////////////////////////////////////////////////////

//   Teuchos::RCP<Epetra_Vector> disp =
//     extractor_.ExtractVector(result.getEpetraVector(),0);

//   Teuchos::RCP<Epetra_Vector> idisp =
//     structurefield_.Interface().ExtractCondVector(disp);

  Teuchos::RCP<Epetra_Vector> idisp = Teuchos::rcp(new Epetra_Vector(*structurefield_.Interface().CondMap()));

  Teuchos::RCP<LinearPartitioned> partitionedInterface =
    Teuchos::rcp(new LinearPartitioned(algorithm_,structurefield_,fluidfield_,alefield_));

  Teuchos::RCP<NOX::Epetra::Interface::Required> partitioned = partitionedInterface;

  // Create the linear system
  // actually, we do not need any
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys = Teuchos::null;

  // Create the Group
  Teuchos::RCP<NOX::Epetra::Group> grp =
    Teuchos::rcp(new NOX::Epetra::Group(printParams, partitioned, NOX::Epetra::Vector(*idisp), linSys));

  // Convergence Tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo = CreateStatusTest(nlParams, grp);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp,combo,RCP<ParameterList>(&nlParams,false));

  // solve the whole thing
  NOX::StatusTest::StatusType status = solver->solve();

  if (status != NOX::StatusTest::Converged)
    utils_.out() << RED "linear partitioned solver failed to converge!" END_COLOR << "\n";

  const NOX::Epetra::Group& finalGroup = Teuchos::dyn_cast<const NOX::Epetra::Group>(solver->getSolutionGroup());

  partitionedInterface->ExtractResult(
    Teuchos::rcp(&Teuchos::dyn_cast<const NOX::Epetra::Vector>(finalGroup.getX()).getEpetraVector(),false),
    result.getEpetraVector());

  maxit = nlParams.sublist("Output").get("Nonlinear Iterations",0);
  tol = nlParams.sublist("Output").get("2-Norm of Residual", 0.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo>
NOX::FSI::LinearPartitionedSolver::CreateStatusTest(ParameterList& nlParams,
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
NOX::FSI::LinearPartitionedSolver::CreateStatusTest(ParameterList& nlParams,
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
void NOX::FSI::LinearPartitionedSolver::throwError(const string& functionName,
                                                   const string& errorMsg) const
{
  if (utils_.isPrintType(NOX::Utils::Error))

  {
    utils_.out() << "NOX::FSI::LinearPartitionedSolver::" << functionName << " - "
                 << errorMsg << endl;
  }
  throw "NOX Error";
}


#endif
