
#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

#include "fsi_debugwriter.H"
#include "fsi_monolithic.H"
#include "fsi_nox_aitken.H"
#include "fsi_nox_group.H"
#include "fsi_nox_newton.H"
#include "fsi_statustest.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_lib/drt_colors.H"

#include "../drt_io/io_control.H"

#ifdef PARALLEL
#include <mpi.h>
#endif



/*----------------------------------------------------------------------*/
// Note: The order of calling the three BaseAlgorithm-constructors is
// important here! In here control file entries are written. And these
// entries define the order in which the filters handle the
// Discretizations, which in turn defines the dof number ordering of the
// Discretizations.
/*----------------------------------------------------------------------*/
FSI::MonolithicBase::MonolithicBase(Epetra_Comm& comm)
  : AlgorithmBase(comm,DRT::Problem::Instance()->FSIDynamicParams()),
    StructureBaseAlgorithm(DRT::Problem::Instance()->FSIDynamicParams()),
    FluidBaseAlgorithm(DRT::Problem::Instance()->FSIDynamicParams(),true),
    AleBaseAlgorithm(DRT::Problem::Instance()->FSIDynamicParams())
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicBase::~MonolithicBase()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBase::ReadRestart(int step)
{
  StructureField().ReadRestart(step);
  FluidField()    .ReadRestart(step);
  AleField()      .ReadRestart(step);

  SetTimeStep(FluidField().Time(),FluidField().Step());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBase::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  StructureField().PrepareTimeStep();
  FluidField().    PrepareTimeStep();
  AleField().      PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBase::Update()
{
  StructureField().Update();
  FluidField().    Update();
  AleField().      Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBase::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  StructureField().Output();
  FluidField().    Output();
  AleField().      Output();

  FluidField().LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::StructToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::StructToAle(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToStruct(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::StructToFluid(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::FluidToStruct(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToFluid(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_.SlaveToMaster(iv);
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Monolithic::Monolithic(Epetra_Comm& comm)
  : MonolithicBase(comm)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  // enable debugging
  if (Teuchos::getIntegralValue<int>(fsidyn,"DEBUGOUTPUT"))
  {
    sdbg_ = Teuchos::rcp(new UTILS::DebugWriter(StructureField().Discretization()));
    //fdbg_ = Teuchos::rcp(new UTILS::DebugWriter(FluidField().Discretization()));
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
  // make sure we didn't destroy the maps before we entered the timeloop
  Extractor().CheckForValidMapExtractor();

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = NOXParameterList();

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  //Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  //Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm().MyPID());

#if 0
  // there is a strange NOX bug...
  Teuchos::ParameterList& oo = printParams.sublist("Output Information");

  oo.set<bool>("Error",false);
  oo.set<bool>("Warning",false);
  oo.set<bool>("Outer Iteration",false);
  oo.set<bool>("Inner Iteration",false);
  oo.set<bool>("Parameters",false);
  oo.set<bool>("Details",false);
  oo.set<bool>("Outer Iteration StatusTest",false);
  oo.set<bool>("Linear Solver Details",false);
  oo.set<bool>("Test Details",false);
  oo.set<bool>("Stepper Iteration",false);
  oo.set<bool>("Stepper Details",false);
  oo.set<bool>("Stepper Parameters",false);
  oo.set<bool>("Debug",false);
#else
  printParams.set("Output Information",
                  NOX::Utils::Error |
                  NOX::Utils::Warning |
                  NOX::Utils::OuterIteration |
                  NOX::Utils::InnerIteration |
                  //NOX::Utils::Parameters |
                  NOX::Utils::Details |
                  NOX::Utils::OuterIterationStatusTest |
                  NOX::Utils::LinearSolverDetails |
                  NOX::Utils::TestDetails |
                  NOX::Utils::StepperIteration |
                  NOX::Utils::StepperDetails |
                  NOX::Utils::StepperParameters |
                  NOX::Utils::Debug |
                  0);
#endif

  // Create printing utilities
  utils_ = Teuchos::rcp(new NOX::Utils(printParams));

  Teuchos::RefCountPtr<std::ofstream> log;
  if (Comm().MyPID()==0)
  {
    std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
    s.append(".iteration");
    log = Teuchos::rcp(new std::ofstream(s.c_str()));
    (*log) << "# num procs      = " << Comm().NumProc() << "\n"
           << "# Method         = " << nlParams.sublist("Direction").get("Method","Newton") << "\n"
           << "#\n"
      ;
  }

  Teuchos::Time timer("time step timer");

  while (NotFinished())
  {
    PrepareTimeStep();

    if (sdbg_!=Teuchos::null)
      sdbg_->NewTimeStep(Step(),"struct");
    if (fdbg_!=Teuchos::null)
      fdbg_->NewTimeStep(Step(),"fluid");

    // start time measurement
    Teuchos::RefCountPtr<Teuchos::TimeMonitor> timemonitor = rcp(new Teuchos::TimeMonitor(timer,true));

    // calculate initial linear system at current position
    // (no increment)
    // This initializes the field algorithms and creates the first linear
    // systems. And this is the reason we know the initial linear system is
    // there when we create the NOX::Group.
    Evaluate(Teuchos::null);

    // Get initial guess.
    // The initial system is there, so we can happily extract the
    // initial guess. (The Dirichlet conditions are already build in!)
    Teuchos::RCP<Epetra_Vector> initial_guess = Teuchos::rcp(new Epetra_Vector(*DofRowMap()));
    InitialGuess(initial_guess);

    NOX::Epetra::Vector noxSoln(initial_guess, NOX::Epetra::Vector::CreateView);

    // Create the linear system
    Teuchos::RCP<NOX::Epetra::LinearSystem> linSys =
      CreateLinearSystem(nlParams, noxSoln, utils_);

    // Create the Group
    Teuchos::RCP<NOX::FSI::Group> grp =
      Teuchos::rcp(new NOX::FSI::Group(*this, printParams, interface, noxSoln, linSys));

    // Convergence Tests
    Teuchos::RCP<NOX::StatusTest::Combo> combo = CreateStatusTest(nlParams, grp);

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp,combo,RCP<ParameterList>(&nlParams,false));

    // we know we already have the first linear system calculated
    grp->CaptureSystemState();

    // solve the whole thing
    NOX::StatusTest::StatusType status = solver->solve();

    if (status != NOX::StatusTest::Converged)
      if (Comm().MyPID()==0)
        utils_->out() << RED "Nonlinear solver failed to converge!" END_COLOR << endl;

    // cleanup
    //mat_->Zero();

    // stop time measurement
    timemonitor = Teuchos::null;

    if (Comm().MyPID()==0)
    {
      (*log) << Step()
             << " " << timer.totalElapsedTime()
             << " " << nlParams.sublist("Output").get("Nonlinear Iterations",0)
             << " " << nlParams.sublist("Output").get("2-Norm of Residual", 0.)
             << " " << lsParams.sublist("Output").get("Total Number of Linear Iterations",0)
        ;
      (*log) << std::endl;
      lsParams.sublist("Output").set("Total Number of Linear Iterations",0);
    }

    Update();
    Output();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::Evaluate");

  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  Teuchos::RCP<const Epetra_Vector> ax;

  if (x!=Teuchos::null)
  {
    // This does not seem very reasonable.
    // It is a premature optimization.
#if 0
    double norm;
    int err = x->Norm2(&norm);
    if (err)
      dserror("failed to calculate norm");

    if (norm==0.)
      return;
#endif

    ExtractFieldVectors(x,sx,fx,ax);

    if (sdbg_!=Teuchos::null)
    {
      sdbg_->NewIteration();
      sdbg_->WriteVector("x",*StructureField().Interface().ExtractFSICondVector(sx));
    }
  }

  // Call all elements and assemble rhs and matrices
  // We only need the rhs here because NOX will ask for the rhs
  // only. But the Jacobian is stored internally and will be returnd
  // later on without looking at x again!

  Utils()->out() << "\nEvaluate elements\n";

  {
    Epetra_Time ts(Comm());
    StructureField().Evaluate(sx);
    Utils()->out() << "structure: " << ts.ElapsedTime() << "\n";
  }

  {
    Epetra_Time ta(Comm());
    AleField()      .Evaluate(ax);
    Utils()->out() << "ale      : " << ta.ElapsedTime() << "\n";
  }

  // transfer the current ale mesh positions to the fluid field
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluid(AleField().ExtractDisplacement());
  FluidField().ApplyMeshDisplacement(fluiddisp);

  {
    Epetra_Time tf(Comm());
    FluidField().Evaluate(fx);
    Utils()->out() << "fluid    : " << tf.ElapsedTime() << "\n";
  }

  Utils()->out() << "\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps)
{
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  blockrowdofmap_.Setup(*fullmap,maps);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::SetDefaultParameters(const Teuchos::ParameterList& fsidyn,
                                           Teuchos::ParameterList& list)
{
  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = list;

  nlParams.set<std::string>("Nonlinear Solver", "Line Search Based");
  //nlParams.set("Preconditioner", "None");
  //nlParams.set("Norm abs F", fsidyn.get<double>("CONVTOL"));

  nlParams.set("Max Iterations", fsidyn.get<int>("ITEMAX"));
  //nlParams.set("Max Iterations", 1);

  nlParams.set("Norm abs pres", fsidyn.get<double>("CONVTOL"));
  nlParams.set("Norm abs vel",  fsidyn.get<double>("CONVTOL"));
  nlParams.set("Norm abs disp", fsidyn.get<double>("CONVTOL"));

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  //Teuchos::ParameterList& lineSearchParams = nlParams.sublist("Line Search");



#if 0 // NonlinearCG

  dirParams.set("Method","NonlinearCG");
  lineSearchParams.set("Method","NonlinearCG");
  Teuchos::ParameterList& nlcgParams = dirParams.sublist("Nonlinear CG");
  ParameterList& lsParams = nlcgParams.sublist("Linear Solver");
  nlcgParams.set("Restart Frequency",10);
  nlcgParams.set("Precondition", "On");
  nlcgParams.set("Orthogonalize", "Fletcher-Reeves");
  lsParams.set("Preconditioning", "User Supplied Preconditioner");
  lsParams.set("Preconditioner","User Defined");

#endif // NonlinearCG

#if 0 // Steepest Decent

  dirParams.set("Method","Steepest Descent");
  //lineSearchParams.set("Method","Full Step");
  Teuchos::ParameterList& SDParams = dirParams.sublist("Steepest Descent");
  ParameterList& lsParams = SDParams.sublist("Linear Solver");
  dirParams.set("Precondition","On");
  SDParams.set("Precondition","On"); // not recognized
  // Quadratic Model Min / 2-Norm / F 2-Norm / None
  SDParams.set("Scaling Type","2-Norm");
  lsParams.set("Preconditioning","User Supplied Preconditioner");
  lsParams.set("Preconditioner","User Defined");

#endif // Steepest Decent


#if 1 // Newton

  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  dirParams.set<std::string>("Method","User Defined");
  Teuchos::RCP<NOX::Direction::UserDefinedFactory> newtonfactory = Teuchos::rcp(this,false);
  dirParams.set("User Defined Direction Factory",newtonfactory);


#if 0 // polynominal line search
  lineSearchParams.set("Method","Polynomial");
  Teuchos::ParameterList& polyparams = lineSearchParams.sublist("Polynomial");
  polyparams.set("Default Step",1.0);
  polyparams.set("Max Iters",5);
  polyparams.set("Minimum Step",0.01);
  polyparams.set("Maximum Step",1.5);
  polyparams.set("Recovery Step Type","Constant");
  polyparams.set("Recovery Step",0.1);
  //polyparams.set("Interpolation Type","Cubic");
  polyparams.set("Interpolation Type","Quadratic");
  //polyparams.set("Interpolation Type","Quadratic3");
  //polyparams.set("Sufficient Decrease Condition","Armijo-Goldstein");
  polyparams.set("Sufficient Decrease Condition","Ared/Pred");
  //polyparams.set("Sufficient Decrease Condition","None");
  //polyparams.set("Sufficient Decrease",0.001);
  polyparams.set("Alpha Factor",0.1);
  polyparams.set("Use Counters",false);
  polyparams.set("Maximum Iteration for Increase",0);
  polyparams.set("Optimize Slope Calculation",false);
  polyparams.set("Allowed Relative Increase",1.0);
  polyparams.set("Min Bounds Factor",0.1);
  polyparams.set("Max Bounds Factor",0.95);
#endif

#if 0 // backtracking
  lineSearchParams.set("Method","Backtrack");
  Teuchos::ParameterList& backtrack = lineSearchParams.sublist("Backtrack");
  backtrack.set("Default Step",1.0);
  backtrack.set("Minimum Step",0.01);
  backtrack.set("Recovery Step",0.1);
  backtrack.set("Reduction Factor",0.5);
  backtrack.set("Max Iters",10);
  backtrack.set("Maximum Iteration for Increase",0);
  backtrack.set("Optimize Slope Calculation",false);
  polyparams.set("Min Bounds Factor",0.1);
  polyparams.set("Max Bounds Factor",1.5);
#endif

#if 0 // More'-Thuente
  lineSearchParams.set("Method","More'-Thuente");
  Teuchos::ParameterList& mt = lineSearchParams.sublist("More'-Thuente");
  mt.set("Default Step",1.0);
  mt.set("Minimum Step",0.01);
  mt.set("Maximum Step",1.5);
  mt.set("Max Iters",5);
  mt.set("Recovery Step",0.1);
  //mt.set("Sufficient Decrease Condition","Armijo-Goldstein");
  mt.set("Sufficient Decrease Condition","Ared/Pred");
  mt.set("Alpha Factor",0.1);
  mt.set("Optimize Slope Calculation",false);
  mt.set("Maximum Iteration for Increase",0);
  mt.set("Curvature Condition",0.9999);
  mt.set("Min Bounds Factor",0.1);
  mt.set("Max Bounds Factor",1.1);
#endif


#endif // Newton

  // status tests are expensive, but instructive
  //solverOptions.set<std::string>("Status Test Check Type","Minimal");
  solverOptions.set<std::string>("Status Test Check Type","Complete");

  // be explicit about linear solver parameters
  lsParams.set<std::string>("Aztec Solver","GMRES");
  //lsParams.set<std::string>("BiCGStab","GMRES");
  lsParams.set<std::string>("Orthogonalization","Modified");

  // "r0", "rhs", "norm", "no scaling", "sol"
  lsParams.set<std::string>("Convergence Test","r0");

  lsParams.set<int>("Size of Krylov Subspace",50);
  lsParams.set<int>("Max Iterations",2500);
  lsParams.set<std::string>("Preconditioner","User Defined");
  lsParams.set<int>("Output Frequency",10);
  lsParams.set<bool>("Output Solver Details",true);

  // adaptive tolerance settings
  lsParams.set<double>("base tolerance",fsidyn.get<double>("BASETOL"));
  lsParams.set<double>("adaptive distance",fsidyn.get<double>("ADAPTIVEDIST"));

#if 0
  // add Aitken relaxation to Newton step
  // there is nothing to be gained...
  Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> aitkenfactory =
    Teuchos::rcp(new NOX::FSI::AitkenFactory());
  lineSearchParams.set("Method","User Defined");
  lineSearchParams.set("User Defined Line Search Factory", aitkenfactory);

  //lineSearchParams.sublist("Aitken").set("max step size",
  //fsidyn.get<double>("MAXOMEGA"));
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Direction::Generic>
FSI::Monolithic::buildDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
                                Teuchos::ParameterList& params) const
{
  Teuchos::RCP<NOX::FSI::Newton> newton = Teuchos::rcp(new NOX::FSI::Newton(gd,params));
  for (unsigned i=0; i<statustests_.size(); ++i)
  {
    statustests_[i]->SetNewton(newton);
  }
  return newton;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::Monolithic::computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::computeF");
  Evaluate(Teuchos::rcp(&x,false));
  SetupRHS(F);
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::Monolithic::computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac)
{
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::Monolithic::computePreconditioner(const Epetra_Vector &x,
                                            Epetra_Operator &M,
                                            Teuchos::ParameterList *precParams)
{
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::BlockMonolithic::BlockMonolithic(Epetra_Comm& comm)
  : Monolithic(comm),
    precondreusecount_(0)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::BlockMonolithic::computeJacobian(const Epetra_Vector &x, Epetra_Operator &Jac)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::BlockMonolithic::computeJacobian");
  Evaluate(Teuchos::rcp(&x,false));
  LINALG::BlockSparseMatrixBase& mat = Teuchos::dyn_cast<LINALG::BlockSparseMatrixBase>(Jac);
  SetupSystemMatrix(mat);
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::BlockMonolithic::computePreconditioner(const Epetra_Vector &x,
                                            Epetra_Operator &M,
                                            Teuchos::ParameterList *precParams)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::BlockMonolithic::computePreconditioner");

  if (precondreusecount_<=0)
  {
    // Create preconditioner operator. The blocks are already there. This is
    // the perfect place to initialize the block preconditioners.
    SystemMatrix()->SetupPreconditioner();

    const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();
    precondreusecount_ = fsidyn.get<int>("PRECONDREUSE");
  }

  precondreusecount_ -= 1;

  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::BlockMonolithic::PrepareTimeStep()
{
  FSI::Monolithic::PrepareTimeStep();

  // new time step, rebuild preconditioner
  precondreusecount_ = 0;
}


#endif
