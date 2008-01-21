
#ifdef CCADISCRET

#include "mfsi_algorithm.H"
#include "mfsi_nox_thyra_group.H"
#include "mfsi_preconditionerfactory.H"
#include "mfsi_statustest.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_colors.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include <NOX.H>

#include <Thyra_EpetraThyraWrappers.hpp>
#include <Thyra_EpetraLinearOp.hpp>
#include <Thyra_get_Epetra_Operator.hpp>

#include <Thyra_VectorStdOps.hpp>
#include <Thyra_DefaultIdentityLinearOp.hpp>

// fix clashes between ccarat and Thyra::AmesosLinearOpWithSolveFactory
#ifdef UMFPACK
#undef UMFPACK
#endif

#include <Thyra_DefaultRealLinearSolverBuilder.hpp>
#include <Thyra_AmesosLinearOpWithSolveFactory.hpp>
#include <Thyra_AztecOOLinearOpWithSolveFactory.hpp>
#include <Thyra_IfpackPreconditionerFactory.hpp>
#include <Thyra_MLPreconditionerFactory.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C"
{
#include "../headers/standardtypes.h"
}


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;

/*----------------------------------------------------------------------*
 | global variable *solv, vector of lenght numfld of structures SOLVAR  |
 | defined in solver_control.c                                          |
 |                                                                      |
 |                                                       m.gee 11/00    |
 *----------------------------------------------------------------------*/
extern struct _SOLVAR  *solv;

/*!----------------------------------------------------------------------
\brief file pointers

<pre>                                                         m.gee 8/00
This structure struct _FILES allfiles is defined in input_control_global.c
and the type is in standardtypes.h
It holds all file pointers and some variables needed for the FRSYSTEM
</pre>
*----------------------------------------------------------------------*/
extern struct _FILES  allfiles;

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | structure of flags to control output                                 |
 | defined in out_global.c                                              |
 *----------------------------------------------------------------------*/
extern struct _IO_FLAGS     ioflags;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
extern Teuchos::RCP<Teuchos::ParameterList> globalparameterlist;


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
MFSI::Algorithm::Algorithm(Epetra_Comm& comm)
  : comm_(comm)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  step_ = 0;
  time_ = 0.;
  dt_ = fsidyn.get<double>("TIMESTEP");
  nstep_ = fsidyn.get<int>("NUMSTEP");
  maxtime_ = fsidyn.get<double>("MAXTIME");

  SetupStructure();
  SetupFluid();
  SetupAle();

  // Create the system matrix. Right now it is empty. It is filled in
  // create_W_op().
  mat_ = Teuchos::rcp(new Thyra::DefaultBlockedLinearOp<double>());

  //--------------------------------------------------
  evaluatetimer_ = Teuchos::TimeMonitor::getNewTimer("MFSI::Algorithm::evalModelImpl");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::SetupStructure()
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("MFSI::Algorithm::SetupStructure");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf,0);

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RefCountPtr<IO::DiscretizationWriter> output =
    rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR*         actsolv  = &solv[genprob.numsf];

  const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& sdyn     = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  RefCountPtr<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // create a generalized alpha time integrator
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> genalphaparams = rcp(new ParameterList());
  StruGenAlpha::SetDefaults(*genalphaparams);

  genalphaparams->set<bool>  ("damping",Teuchos::getIntegralValue<int>(sdyn,"DAMPING"));
  genalphaparams->set<double>("damping factor K",sdyn.get<double>("K_DAMP"));
  genalphaparams->set<double>("damping factor M",sdyn.get<double>("M_DAMP"));

  genalphaparams->set<double>("beta",sdyn.get<double>("BETA"));
  genalphaparams->set<double>("gamma",sdyn.get<double>("GAMMA"));
  genalphaparams->set<double>("alpha m",sdyn.get<double>("ALPHA_M"));
  genalphaparams->set<double>("alpha f",sdyn.get<double>("ALPHA_F"));

  genalphaparams->set<double>("total time",0.0);
  genalphaparams->set<double>("delta time",fsidyn.get<double>("TIMESTEP"));
  genalphaparams->set<int>   ("step",0);
  genalphaparams->set<int>   ("nstep",fsidyn.get<int>("NUMSTEP"));
  genalphaparams->set<int>   ("max iterations",sdyn.get<int>("MAXITER"));
  genalphaparams->set<int>   ("num iterations",-1);
  genalphaparams->set<double>("tolerance displacements",sdyn.get<double>("TOLDISP"));

  genalphaparams->set<bool>  ("io structural disp",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_DISP"));
  genalphaparams->set<int>   ("io disp every nstep",fsidyn.get<int>("UPRES"));
  genalphaparams->set<bool>  ("io structural stress",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_STRESS"));
  genalphaparams->set<int>   ("io stress every nstep",sdyn.get<int>("RESEVRYSTRS"));

  genalphaparams->set<int>   ("restart",probtype.get<int>("RESTART"));
  genalphaparams->set<int>   ("write restart every",fsidyn.get<int>("RESTARTEVRY"));

  genalphaparams->set<bool>  ("print to screen",true);
  genalphaparams->set<bool>  ("print to err",true);
  genalphaparams->set<FILE*> ("err file",allfiles.out_err);

  // takes values "full newton" , "modified newton" , "nonlinear cg"
  genalphaparams->set<string>("equilibrium iteration","full newton");

  structure_ = rcp(new StructureAdapter(genalphaparams,actdis,solver,output));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::SetupFluid()
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("MFSI::Algorithm::SetupFluid");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numff,0);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RefCountPtr<IO::DiscretizationWriter> output =
    rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR        *actsolv  = &solv[genprob.numff];

  const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  RefCountPtr<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // create a fluid nonlinear time integrator
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> fluidtimeparams = rcp(new ParameterList());
  FluidImplicitTimeInt::SetDefaults(*fluidtimeparams);

  FLUID_TIMEINTTYPE iop = Teuchos::getIntegralValue<FLUID_TIMEINTTYPE>(fdyn,"TIMEINTEGR");

  // number of degrees of freedom
  fluidtimeparams->set<int>              ("number of velocity degrees of freedom" ,probsize.get<int>("DIM"));
  // the default time step size
  fluidtimeparams->set<double>           ("time step size"           ,fsidyn.get<double>("TIMESTEP"));
  // max. sim. time
  fluidtimeparams->set<double>           ("total time"               ,fsidyn.get<double>("MAXTIME"));
  // parameter for time-integration
  fluidtimeparams->set<double>           ("theta"                    ,fdyn.get<double>("THETA"));
  // which kind of time-integration
  fluidtimeparams->set<FLUID_TIMEINTTYPE>("time int algo"            ,iop);
  // bound for the number of timesteps
  fluidtimeparams->set<int>              ("max number timesteps"     ,fsidyn.get<int>("NUMSTEP"));
  // number of steps with start algorithm
  fluidtimeparams->set<int>              ("number of start steps"    ,fdyn.get<int>("NUMSTASTEPS"));
  // parameter for start algo
  fluidtimeparams->set<double>           ("start theta"              ,fdyn.get<double>("START_THETA"));


  // ---------------------------------------------- nonlinear iteration
  // set linearisation scheme
  fluidtimeparams->set<bool>("Use reaction terms for linearisation",
                            Teuchos::getIntegralValue<int>(fdyn,"NONLINITER")==2);
  // maximum number of nonlinear iteration steps
  fluidtimeparams->set<int>             ("max nonlin iter steps"     ,fdyn.get<int>("ITEMAX"));
  // stop nonlinear iteration when both incr-norms are below this bound
  fluidtimeparams->set<double>          ("tolerance for nonlin iter" ,fdyn.get<double>("CONVTOL"));

  // ----------------------------------------------- restart and output
  // restart
  fluidtimeparams->set                 ("write restart every"       ,fsidyn.get<int>("RESTARTEVRY"));
  // solution output
  fluidtimeparams->set                 ("write solution every"      ,fsidyn.get<int>("UPRES"));
  // flag for writing stresses
  fluidtimeparams->set                 ("write stresses"            ,Teuchos::getIntegralValue<int>(ioflags,"FLUID_STRESS"));

  //--------------------------------------------------
  // evaluate error for test flows with analytical solutions
  int init = Teuchos::getIntegralValue<int>(fdyn,"INITIALFIELD");
  fluidtimeparams->set                  ("eval err for analyt sol"   ,init);


  //--------------------------------------------------
  // create all vectors and variables associated with the time
  // integration (call the constructor)
  // the only parameter from the list required here is the number of
  // velocity degrees of freedom
  fluid_ = rcp(new FluidAdapter(actdis, solver, fluidtimeparams, output));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::SetupAle()
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("MFSI::Algorithm::SetupAle");
  Teuchos::TimeMonitor monitor(*t);

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numaf,0);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RefCountPtr<IO::DiscretizationWriter> output =
    rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR        *actsolv  = &solv[genprob.numaf];

  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  RefCountPtr<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  RefCountPtr<ParameterList> params = rcp(new ParameterList());
  params->set<int>("numstep",    fsidyn.get<int>("NUMSTEP"));
  params->set<double>("maxtime", fsidyn.get<double>("MAXTIME"));
  params->set<double>("dt",      fsidyn.get<double>("TIMESTEP"));

  ale_ = rcp(new FSI::AleLinear(actdis, solver, params, output, false));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::Timeloop()
{
  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = *globalparameterlist;

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method","Newton"));
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  //Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm().MyPID());

  // turn on output
  printParams.set("Output Information", 0xffff);

  // Create printing utilities
  utils_ = Teuchos::rcp(new NOX::Utils(printParams));

  Teuchos::RefCountPtr<std::ofstream> log;
  if (Comm().MyPID()==0)
  {
    std::string s = allfiles.outputfile_kenner;
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

    // start time measurement
    Teuchos::RefCountPtr<Teuchos::TimeMonitor> timemonitor = rcp(new Teuchos::TimeMonitor(timer,true));

    // calculate initial linear system at current position
    // (no increment)
    Evaluate(Teuchos::null);

    // Get initial guess.
    // The initial system is there, so we can happily extract the
    // initial guess. (The Dirichlet conditions are already build in!)
    Thyra::DefaultProductVector<double> initial_guess(dofrowmap_);
    InitialGuess(initial_guess);

    // Create the NOX Group
    // The initial system is already known. It must not be reevaluated!
    Teuchos::RCP<NOX::Thyra::Group> grp = Teuchos::rcp(new NOX_Thyra_Group(initial_guess, this));

    // Convergence Tests
    Teuchos::RCP<NOX::StatusTest::Combo> combo = CreateStatusTest(nlParams, grp);

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp,combo,RCP<ParameterList>(&nlParams,false));

    // solve the whole thing
    NOX::StatusTest::StatusType status = solver->solve();

    if (status != NOX::StatusTest::Converged)
      if (Comm().MyPID()==0)
        utils_->out() << RED "Nonlinear solver failed to converge!" END_COLOR << endl;

    // cleanup
    // Do not keep the block matrices. They are too heavy! And the
    // next Evaluate() call will replace them anyway.
    mat_->uninitialize();

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
    }

    Update();
    Output();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Thyra::ModelEvaluatorBase::InArgs<double> MFSI::Algorithm::createInArgs() const
{
  // Here we define what kinds of input arguments
  // MFSI::OverlapAlgorithm::evalModelImpl() supports
  Thyra::ModelEvaluatorBase::InArgsSetup<double> inArgs;
  inArgs.setModelEvalDescription(this->description());
  inArgs.setSupports(Thyra::ModelEvaluatorBase::IN_ARG_x);
  return inArgs;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Thyra::ModelEvaluatorBase::OutArgs<double> MFSI::Algorithm::createOutArgsImpl() const
{
  // Here we define what result types
  // MFSI::OverlapAlgorithm::evalModelImpl() supports
  Thyra::ModelEvaluatorBase::OutArgsSetup<double> outArgs;
  outArgs.setModelEvalDescription(this->description());
  outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_f);
  outArgs.setSupports(Thyra::ModelEvaluatorBase::OUT_ARG_W);
  return outArgs;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<double> &inArgs,
                                    const Thyra::ModelEvaluatorBase::OutArgs<double> &outArgs) const
{
  Teuchos::TimeMonitor monitor(*evaluatetimer_);

  const Thyra::VectorBase<double> &x_bar = *inArgs.get_x();
  const Thyra::DefaultProductVector<double> &x = dynamic_cast<const Thyra::DefaultProductVector<double>&>(x_bar);

  Evaluate(Teuchos::rcp(&x,false));

  Teuchos::RCP<Thyra::VectorBase<double> > f_bar = outArgs.get_f();
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > W_bar = outArgs.get_W();

  if (f_bar!=Teuchos::null)
  {
    Thyra::DefaultProductVector<double> &f = dynamic_cast<Thyra::DefaultProductVector<double>&>(*f_bar);
    SetupRHS(f);
  }

  if (W_bar!=Teuchos::null)
  {
    // We know that W_bar contains mat_. We cannot extract it,
    // though. We could set it again.
    SetupSysMat(*SysMat());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > MFSI::Algorithm::create_W() const
{
  Teuchos::RCP<Thyra::LinearOpWithSolveBase<double> > solver = solverfactory_->createOp();

  // initialize the solver with the composite matrix
  solverfactory_->initializeOp(Thyra::defaultLinearOpSource<double,double>(create_W_op()),
                               &*solver,
                               Thyra::SUPPORT_SOLVE_FORWARD_ONLY);
  return solver;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Thyra::LinearOpBase<double> > MFSI::Algorithm::create_W_op() const
{
  SetupSysMat(*mat_);
  return mat_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::Evaluate(Teuchos::RCP<const Thyra::DefaultProductVector<double> > x) const
{
  if (x!=Teuchos::null)
  {
    if (Thyra::norm(*x)!=0)
    {
      Utils()->out() << YELLOW_LIGHT "element call with new x" END_COLOR << endl;

      Teuchos::RCP<const Epetra_Vector> sx;
      Teuchos::RCP<const Epetra_Vector> fx;
      Teuchos::RCP<const Epetra_Vector> ax;

      ExtractFieldVectors(x,sx,fx,ax);

      // debug
      debug_.DumpVector("sx",*StructureField()->Discretization(),*sx);
      debug_.DumpVector("fx",*FluidField()->Discretization(),*fx);
      debug_.DumpVector("ax",*AleField()->Discretization(),*ax);

      // Call all elements and assemble rhs and matrices
      // We only need the rhs here because NOX will ask for the rhs
      // only. But the Jacobian is stored internally and will be returnd
      // later on without looking at x again!
      StructureField()->Evaluate(sx);
      AleField()      ->Evaluate(ax);

      // transfer the current ale mesh positions to the fluid field
      Teuchos::RefCountPtr<Epetra_Vector> fluiddisp = AleToFluid(AleField()->ExtractDisplacement());
      FluidField()->ApplyMeshDisplacement(fluiddisp);

      FluidField()    ->Evaluate(fx);
    }
  }
  else
  {
    Utils()->out() << YELLOW_LIGHT "element call at current x" END_COLOR << endl;

    StructureField()->Evaluate(Teuchos::null);
    AleField()      ->Evaluate(Teuchos::null);

    // transfer the current ale mesh positions to the fluid field
    Teuchos::RefCountPtr<Epetra_Vector> fluiddisp = AleToFluid(AleField()->ExtractDisplacement());
    FluidField()->ApplyMeshDisplacement(fluiddisp);

    FluidField()    ->Evaluate(Teuchos::null);
  }

  //debug.DumpVector("sres",*StructureField()->Discretization(),*StructureField()->RHS());
  //debug.DumpVector("fres",*FluidField()->Discretization(),*FluidField()->RHS());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > MFSI::Algorithm::get_x_space() const
{
  return dofrowmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > MFSI::Algorithm::get_f_space() const
{
  return dofrowmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;

  structure_->PrepareTimeStep();
  fluid_->    PrepareTimeStep();
  ale_->      PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::Update()
{
  // Just go on. Tell each field to update itself. The fields already
  // have all the information they need because each residual
  // evaluation updates the current variable in each field.

  structure_->Update();
  fluid_    ->Update();
  ale_      ->Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::Output()
{
  structure_->Output();
  fluid_    ->Output();
  ale_      ->Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::Algorithm::StructToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::Algorithm::AleToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::Algorithm::StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::Algorithm::FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::Algorithm::AleToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::Algorithm::StructToAle(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::Algorithm::AleToStruct(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::Algorithm::StructToFluid(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_.MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::Algorithm::FluidToStruct(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> MFSI::Algorithm::AleToFluid(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> >
MFSI::Algorithm::CreateSolverFactory(struct _SOLVAR* actsolv)
{
  //ParameterList params;
  //LINALG::Solver::TranslateSolverParameters(params, actsolv);

  Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<double> > solverfactory;

  switch (actsolv->solvertyp)
  {
  case amesos_klu_sym://====================================== Tim Davis' KLU
  case amesos_klu_nonsym://=================================== Tim Davis' KLU
    solverfactory = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory(Thyra::Amesos::KLU));
    break;
  case umfpack://========================================= Tim Davis' Umfpack
    solverfactory = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory(Thyra::Amesos::UMFPACK));
    break;
  case aztec_msr://================================================= AztecOO
  {
    Teuchos::RCP<Thyra::AztecOOLinearOpWithSolveFactory> aztecfactory =
      Teuchos::rcp(new Thyra::AztecOOLinearOpWithSolveFactory());

    AZVAR* azvar = actsolv->azvar;
    Teuchos::RCP<ParameterList> aztecparams = Teuchos::rcp(new ParameterList("Aztec Parameters"));
    //aztecparams.set("Output Every RHS",);
    //aztecparams.set("Adjoint Solve",);

    ParameterList& solveparams = aztecparams->sublist("Forward Solve");
    solveparams.set("Max Iterations",azvar->aziter);
    solveparams.set("Tolerance",azvar->aztol);

    ParameterList& params = solveparams.sublist("AztecOO Settings");

    switch (azvar->azsolvertyp)
    {
    case azsolv_CG:       params.set("Aztec Solver","CG");       break;
    case azsolv_GMRES:    params.set("Aztec Solver","GMRES");    break;
    case azsolv_CGS:      params.set("Aztec Solver","CGS");      break;
    case azsolv_BiCGSTAB: params.set("Aztec Solver","BiCGStab"); break;
    case azsolv_LU:       params.set("Aztec Solver","LU");       break;
    case azsolv_TFQMR:    params.set("Aztec Solver","TFQMR");    break;
    default: dserror("Unknown solver for AztecOO");               break;
    }

    // build preconditioner factory

    switch (azvar->azprectyp)
    {
    case azprec_none:
      params.set("Aztec Preconditioner","none");
      break;
    case azprec_ILUT:
    case azprec_ILU:
    case azprec_LU:
    case azprec_ICC:
    {
      params.set("Aztec Preconditioner","none");
      Teuchos::RCP<ParameterList> pcparams = Teuchos::rcp(new ParameterList("IFPACK Parameters"));
      switch (azvar->azprectyp)
      {
      case azprec_ILUT:  pcparams->set("Prec Type","ILUT");   break;
      case azprec_ILU:   pcparams->set("Prec Type","ILU");    break;
      case azprec_LU:    pcparams->set("Prec Type","Amesos"); break;
      case azprec_ICC:   pcparams->set("Prec Type","IC");     break;
      default:
        dserror("Logical breakdown. Panic.");
      }
      ParameterList& ifpackparams = pcparams->sublist("Ifpack Settings");
      ifpackparams.set("fact: drop tolerance",azvar->azdrop);
      ifpackparams.set("fact: level-of-fill",azvar->azgfill);
      ifpackparams.set("fact: ilut level-of-fill",azvar->azfill);
      ifpackparams.set("schwarz: combine mode","Add"); // can be "Add", "Zero", "Insert", "InsertAdd", "Average", "AbsMax"
      ifpackparams.set("schwarz: reordering type","rcm"); // "rcm" or "metis"
      ifpackparams.set("amesos: solver type", "Amesos_Klu"); // can be "Amesos_Klu", "Amesos_Umfpack", "Amesos_Superlu"

      Teuchos::RCP<Thyra::IfpackPreconditionerFactory > precondfactory =
        Teuchos::rcp(new Thyra::IfpackPreconditionerFactory());
      precondfactory->setParameterList(pcparams);
      aztecfactory->setPreconditionerFactory(precondfactory, "IFPACK Preconditioner");
      break;
    }
    case azprec_Jacobi:
      params.set("Aztec Preconditioner","Jacobi");
      break;
    case azprec_Neumann:
      params.set("Aztec Preconditioner","Polynomial");
      break;
    case azprec_Least_Squares:
      params.set("Aztec Preconditioner","Least-squares Polynomial");
      break;
    case azprec_SymmGaussSeidel:
      params.set("Aztec Preconditioner","Symmetric Gauss-Seidel");
      break;
#if 0
    case azprec_RILU:
      azlist.set("AZ_precond",AZ_dom_decomp);
      azlist.set("AZ_subdomain_solve",AZ_rilu);
      azlist.set("AZ_graph_fill",azvar->azgfill);
      break;
#endif
    case azprec_ML:
    case azprec_MLfluid:
    case azprec_MLAPI:
    case azprec_MLfluid2:
    {
      params.set("Aztec Preconditioner","none");
      Teuchos::RCP<Thyra::MLPreconditionerFactory> mlfactory =
        Teuchos::rcp(new Thyra::MLPreconditionerFactory());
      Teuchos::RCP<ParameterList> mlparams = Teuchos::rcp(new ParameterList("ML Parameters"));
      ML_Epetra::SetDefaults("SA",*mlparams);

      switch (azvar->azprectyp)
      {
      case azprec_ML: // do nothing, this is standard
        break;
      case azprec_MLAPI: // set flag to use mlapi operator
        mlparams->set<bool>("LINALG::AMG_Operator",true);
        break;
      case azprec_MLfluid: // unsymmetric, unsmoothed restruction
        mlparams->set("aggregation: use tentative restriction",true);
        break;
      case azprec_MLfluid2: // full Pretrov-Galerkin unsymmetric smoothed
        mlparams->set("energy minimization: enable",true);
        mlparams->set("energy minimization: type",3); // 1,2,3 cheap -> expensive
        mlparams->set("aggregation: block scaling",false);
        break;
      default: dserror("Unknown type of ml preconditioner");
      }
      mlparams->set("output"                          ,azvar->mlprint);
      if (azvar->mlprint==10)
        mlparams->set("print unused"                  ,1);
      else
        mlparams->set("print unused"                  ,-2);
      mlparams->set("increasing or decreasing"        ,"increasing");
      mlparams->set("coarse: max size"                ,azvar->mlcsize);
      mlparams->set("max levels"                      ,azvar->mlmaxlevel);
      mlparams->set("smoother: pre or post"           ,"both");
      mlparams->set("aggregation: threshold"          ,azvar->ml_threshold);
      mlparams->set("aggregation: damping factor"     ,azvar->mldamp_prolong);
      mlparams->set("aggregation: nodes per aggregate",azvar->mlaggsize);
      switch (azvar->mlcoarsentype)
      {
      case 0:  mlparams->set("aggregation: type","Uncoupled");  break;
      case 1:  mlparams->set("aggregation: type","METIS");      break;
      case 2:  mlparams->set("aggregation: type","VBMETIS");    break;
      case 3:  mlparams->set("aggregation: type","MIS");        break;
      default: dserror("Unknown type of coarsening for ML"); break;
      }
      // set ml smoothers
      for (int i=0; i<azvar->mlmaxlevel-1; ++i)
      {
        char levelstr[11];
        sprintf(levelstr,"(level %d)",i);
        int type;
        double damp;
        if (i==0)
        {
          type = azvar->mlsmotype_fine;
          damp = azvar->mldamp_fine;
        }
        else if (i < azvar->mlmaxlevel-1)
        {
          type = azvar->mlsmotype_med;
          damp = azvar->mldamp_med;
        }
        else
        {
          type = azvar->mlsmotype_coarse;
          damp = azvar->mldamp_coarse;
        }
        switch (type)
        {
        case 0:
          mlparams->set("smoother: type "+(string)levelstr                    ,"symmetric Gauss-Seidel");
          mlparams->set("smoother: sweeps "+(string)levelstr                  ,azvar->mlsmotimes[i]);
          mlparams->set("smoother: damping factor "+(string)levelstr          ,damp);
          break;
        case 1:
          mlparams->set("smoother: type "+(string)levelstr                    ,"Jacobi");
          mlparams->set("smoother: sweeps "+(string)levelstr                  ,azvar->mlsmotimes[i]);
          mlparams->set("smoother: damping factor "+(string)levelstr          ,damp);
          break;
        case 2:
          mlparams->set("smoother: type "+(string)levelstr                    ,"MLS");
          mlparams->set("smoother: MLS polynomial order "+(string)levelstr    ,azvar->mlsmotimes[i]);
          break;
        case 3:
          mlparams->set("smoother: type (level 0)"                            ,"MLS");
          mlparams->set("smoother: MLS polynomial order "+(string)levelstr    ,-azvar->mlsmotimes[i]);
          break;
        case 4:
        {
          mlparams->set("smoother: type "+(string)levelstr                    ,"IFPACK");
          mlparams->set("smoother: ifpack type "+(string)levelstr             ,"ILU");
          mlparams->set("smoother: ifpack overlap "+(string)levelstr          ,0);
          ParameterList& ifpacklist = mlparams->sublist("smoother: ifpack list");
          ifpacklist.set<int>("fact: level-of-fill",azvar->mlsmotimes[i]);
          ifpacklist.set("schwarz: reordering type","rcm");
        }
        break;
        case 5:
          mlparams->set("smoother: type "+(string)levelstr,"Amesos-KLU");
          break;
#ifdef PARALLEL
        case 6:
          mlparams->set("smoother: type "+(string)levelstr,"Amesos-Superludist");
          break;
#endif
        default: dserror("Unknown type of smoother for ML"); break;
        } // switch (type)
      } // for (int i=0; i<azvar->mlmaxlevel-1; ++i)
      // set coarse grid solver
      const int coarse = azvar->mlmaxlevel-1;
      switch (azvar->mlsmotype_coarse)
      {
      case 0:
        mlparams->set("coarse: type"          ,"symmetric Gauss-Seidel");
        mlparams->set("coarse: sweeps"        , azvar->mlsmotimes[coarse]);
        mlparams->set("coarse: damping factor",azvar->mldamp_coarse);
        break;
      case 1:
        mlparams->set("coarse: type"          ,"Jacobi");
        mlparams->set("coarse: sweeps"        , azvar->mlsmotimes[coarse]);
        mlparams->set("coarse: damping factor",azvar->mldamp_coarse);
        break;
      case 2:
        mlparams->set("coarse: type"                ,"MLS");
        mlparams->set("coarse: MLS polynomial order",azvar->mlsmotimes[coarse]);
        break;
      case 3:
        mlparams->set("coarse: type"                ,"MLS");
        mlparams->set("coarse: MLS polynomial order",-azvar->mlsmotimes[coarse]);
        break;
      case 4:
      {
        mlparams->set("coarse: type"          ,"IFPACK");
        mlparams->set("coarse: ifpack type"   ,"ILU");
        mlparams->set("coarse: ifpack overlap",0);
        ParameterList& ifpacklist = mlparams->sublist("coarse: ifpack list");
        ifpacklist.set<int>("fact: level-of-fill",azvar->mlsmotimes[coarse]);
        ifpacklist.set("schwarz: reordering type","rcm");
      }
      break;
      case 5:
        mlparams->set("coarse: type","Amesos-KLU");
        break;
      case 6:
        mlparams->set("coarse: type","Amesos-Superludist");
        break;
      default: dserror("Unknown type of coarse solver for ML"); break;
      } // switch (azvar->mlsmotype_coarse)
      // default values for nullspace
      mlparams->set("PDE equations",1);
      mlparams->set("null space: dimension",1);
      mlparams->set("null space: type","pre-computed");
      mlparams->set("null space: add default vectors",false);
      mlparams->set<double*>("null space: vectors",NULL);

      mlfactory->setParameterList(mlparams);
      aztecfactory->setPreconditionerFactory(mlfactory, "ML Preconditioner");
      break;
    }
    default:
      dserror("Unknown preconditioner for AztecOO");
    break;
    }

    //------------------------------------- set other aztec parameters
    params.set("Overlap",0);
    //params.set("Graph Fill",);
    params.set("Drop Tolerance",azvar->azdrop);
    //params.set("Fill Factor",);
    //params.set("Steps",);
    params.set("Polynomial Order",azvar->azpoly);
    //params.set("RCM Reordering",);
    //params.set("Orthogonalization",AZ_modified); // ???
    params.set("Size of Krylov Subspace",azvar->azsub);

    switch (azvar->azconv)
    {
    case AZ_r0:       params.set("Convergence Test","r0");         break;
    case AZ_rhs:      params.set("Convergence Test","rhs");        break;
    case AZ_Anorm:    params.set("Convergence Test","Anorm");      break;
    case AZ_noscaled: params.set("Convergence Test","no scaling"); break;
    case AZ_sol:      params.set("Convergence Test","sol");        break;
    case AZ_weighted:
    case AZ_expected_values:
    case AZTECOO_conv_test:
    case AZ_inf_noscaled:
    default:
      dserror("unsupported convergance criteria %d",azvar->azconv);
    }

    //params.set("Ill-Conditioning Threshold",);
    params.set("Output Frequency",azvar->azoutput);

    aztecfactory->setParameterList(aztecparams);

    solverfactory = aztecfactory;
    break;
  }
  default:
    dserror("Unsupported type %d of solver",actsolv->solvertyp);
  }

  return solverfactory;
}

#endif
