
#ifdef CCADISCRET

#include "fsi_algorithm.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"

#include <Teuchos_StandardParameterEntryValidators.hpp>


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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Algorithm::Algorithm(Epetra_Comm& comm)
  : comm_(comm)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  if (comm_.MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, fsidyn);

  step_ = 0;
  time_ = 0.;
  dt_ = fsidyn.get<double>("TIMESTEP");
  nstep_ = fsidyn.get<int>("NUMSTEP");
  maxtime_ = fsidyn.get<double>("MAXTIME");

  SetupStructure();
  SetupFluid();
  SetupAle();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Algorithm::~Algorithm()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::ReadRestart(int step)
{
  StructureField().ReadRestart(step);
  FluidField().ReadRestart(step);
  AleField().ReadRestart(step);

  time_ = FluidField().Time();
  step_ = FluidField().Step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::SetupStructure()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numsf,0);

  // set degrees of freedom in the discretization
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter> output =
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

  if (comm_.MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, sdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<ParameterList> solveparams = rcp(new ParameterList());
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // create a generalized alpha time integrator
  // -------------------------------------------------------------------
  RCP<ParameterList> genalphaparams = rcp(new ParameterList());
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
  genalphaparams->set<string>("convcheck",sdyn.get<string>("CONV_CHECK"));

  genalphaparams->set<bool>  ("io structural disp",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_DISP"));
  genalphaparams->set<int>   ("io disp every nstep",fsidyn.get<int>("UPRES"));
  genalphaparams->set<bool>  ("io structural stress",Teuchos::getIntegralValue<int>(ioflags,"STRUCT_STRESS"));
  genalphaparams->set<int>   ("io stress every nstep",sdyn.get<int>("RESEVRYSTRS"));

  genalphaparams->set<int>   ("restart",probtype.get<int>("RESTART"));
  genalphaparams->set<int>   ("write restart every",fsidyn.get<int>("RESTARTEVRY"));

  genalphaparams->set<bool>  ("print to screen",true);
  genalphaparams->set<bool>  ("print to err",true);
  genalphaparams->set<FILE*> ("err file",allfiles.out_err);

  switch (Teuchos::getIntegralValue<int>(sdyn,"NLNSOL"))
  {
  case STRUCT_DYNAMIC::fullnewton:
    genalphaparams->set<string>("equilibrium iteration","full newton");
    break;
  case STRUCT_DYNAMIC::modnewton:
    genalphaparams->set<string>("equilibrium iteration","modified newton");
    break;
  case STRUCT_DYNAMIC::matfreenewton:
    genalphaparams->set<string>("equilibrium iteration","matrixfree newton");
    break;
  case STRUCT_DYNAMIC::nlncg:
    genalphaparams->set<string>("equilibrium iteration","nonlinear cg");
    break;
  case STRUCT_DYNAMIC::ptc:
    genalphaparams->set<string>("equilibrium iteration","ptc");
    break;
  default:
    genalphaparams->set<string>("equilibrium iteration","full newton");
    break;
  }

  // set predictor (takes values "constant" "consistent")
  switch (Teuchos::getIntegralValue<int>(sdyn,"PREDICT"))
  {
  case STRUCT_DYNAMIC::pred_vague:
    dserror("You have to define the predictor");
    break;
  case STRUCT_DYNAMIC::pred_constdis:
    genalphaparams->set<string>("predictor","consistent");
    break;
  case STRUCT_DYNAMIC::pred_constdisvelacc:
    genalphaparams->set<string>("predictor","constant");
    break;
  default:
    dserror("Cannot cope with choice of predictor");
    break;
  }

  structure_ = rcp(new StructureAdapter(genalphaparams,actdis,solver,output));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::SetupFluid()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numff,0);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter> output =
    rcp(new IO::DiscretizationWriter(actdis));
  output->WriteMesh(0,0.0);

  // -------------------------------------------------------------------
  // set some pointers and variables
  // -------------------------------------------------------------------
  SOLVAR        *actsolv  = &solv[genprob.numff];

  //const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& probsize = DRT::Problem::Instance()->ProblemSizeParams();
  const Teuchos::ParameterList& ioflags  = DRT::Problem::Instance()->IOParams();
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  if (comm_.MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, fdyn);

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RCP<ParameterList> solveparams = rcp(new ParameterList());
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  // -------------------------------------------------------------------
  // create a fluid nonlinear time integrator
  // -------------------------------------------------------------------
  RCP<ParameterList> fluidtimeparams = rcp(new ParameterList());

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
  // parameter for grid velocity interpolation
  fluidtimeparams->set<int>              ("order gridvel"            ,fdyn.get<int>("GRIDVEL"));


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

  fluidtimeparams->set<FILE*>("err file",allfiles.out_err);

  //--------------------------------------------------
  // create all vectors and variables associated with the time
  // integration (call the constructor)
  // the only parameter from the list required here is the number of
  // velocity degrees of freedom
  fluid_ = rcp(new FluidAdapter(actdis, solver, fluidtimeparams, output));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::SetupAle()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RCP<DRT::Discretization> actdis = null;
  actdis = DRT::Problem::Instance()->Dis(genprob.numaf,0);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  if (!actdis->Filled()) actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  RCP<IO::DiscretizationWriter> output =
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
  RCP<ParameterList> solveparams = rcp(new ParameterList());
  RCP<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  RCP<ParameterList> params = rcp(new ParameterList());
  params->set<int>("numstep",    fsidyn.get<int>("NUMSTEP"));
  params->set<double>("maxtime", fsidyn.get<double>("MAXTIME"));
  params->set<double>("dt",      fsidyn.get<double>("TIMESTEP"));

  // ----------------------------------------------- restart and output
  // restart
  params->set<int>("write restart every", fsidyn.get<int>("RESTARTEVRY"));

  ale_ = rcp(new AleLinear(actdis, solver, params, output));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;

  if (Comm().MyPID()==0)
    std::cout << "\n"
              << method_ << "\n"
              << "TIME:  "    << std::scientific << time_ << "/" << std::scientific << maxtime_
              << "     DT = " << std::scientific << dt_
              << "     STEP = " YELLOW_LIGHT << setw(4) << step_ << END_COLOR "/" << setw(4) << nstep_
              << "\n"
              << NOX::Utils::fill(82)
              << "\n\n";

  StructureField().PrepareTimeStep();
  FluidField().    PrepareTimeStep();
  AleField().      PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::Update()
{
  StructureField().UpdateandOutput();
  FluidField().    Update();
  AleField().      Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Algorithm::Output()
{
  FluidField().Output();
  AleField().  Output();

  FluidField().LiftDrag();
}


#endif
