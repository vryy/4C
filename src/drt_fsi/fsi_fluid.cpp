#ifdef CCADISCRET

#include "fsi_fluid.H"

// further includes for FluidBaseAlgorithm:
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"
#include <Teuchos_StandardParameterEntryValidators.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>

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
FSI::Fluid::~Fluid()
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidAdapter::FluidAdapter(Teuchos::RCP<DRT::Discretization> dis,
                                 Teuchos::RCP<LINALG::Solver> solver,
                                 Teuchos::RCP<ParameterList> params,
                                 Teuchos::RCP<IO::DiscretizationWriter> output,
                                 bool isale)
  : fluid_(dis, *solver, *params, *output, isale),
    dis_(dis),
    solver_(solver),
    params_(params),
    output_(output)
{
  FSI::UTILS::SetupInterfaceExtractor(*dis,"FSICoupling",interface_);

  // build inner velocity map
  // dofs at the interface are excluded
  // we use only velocity dofs and only those without Dirichlet constraint

  Teuchos::RCP<const Epetra_Map> velmap = fluid_.VelocityRowMap();
  Teuchos::RCP<Epetra_Vector> dirichtoggle = fluid_.Dirichlet();
  Teuchos::RCP<const Epetra_Map> fullmap = DofRowMap();

  int numvelids = velmap->NumMyElements();
  std::vector<int> velids;
  velids.reserve(numvelids);
  for (int i=0; i<numvelids; ++i)
  {
    int gid = velmap->GID(i);
    if (not interface_.CondMap()->MyGID(gid) and (*dirichtoggle)[fullmap->LID(gid)]==0.)
    {
      velids.push_back(gid);
    }
  }

  innervelmap_ = Teuchos::rcp(new Epetra_Map(-1,velids.size(), &velids[0], 0, velmap->Comm()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FSI::FluidAdapter::InitialGuess() const
{
  return fluid_.InitialGuess();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FSI::FluidAdapter::RHS() const
{
  return fluid_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FSI::FluidAdapter::Velnp() const
{
  return fluid_.Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FSI::FluidAdapter::Veln() const
{
  return fluid_.Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FSI::FluidAdapter::DofRowMap() const
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> FSI::FluidAdapter::SystemMatrix() const
{
  return fluid_.SystemMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> FSI::FluidAdapter::Discretization()
{
  return fluid_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
// Teuchos::RCP<Epetra_Vector> FSI::FluidAdapter::StructCondRHS() const
// {
//   return interface_.ExtractCondVector(Velnp());
// }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::PrepareTimeStep()
{
  fluid_.PrepareTimeStep();

  // we add the whole fluid mesh displacement later on?
  //fluid_.Dispnp()->PutScalar(0.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::Evaluate(Teuchos::RCP<const Epetra_Vector> vel) const
{
  // Yes, this is complicated. But we have to be very careful
  // here. The field solver always expects an increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest
  // increment only.
  if (vel!=Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> incvel = Teuchos::rcp(new Epetra_Vector(*vel));
    incvel->Update(-1.0,*fluid_.Velnp(),1.0);
    fluid_.Evaluate(incvel);
  }
  else
  {
    fluid_.Evaluate(Teuchos::null);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::Update()
{
  fluid_.TimeUpdate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::Output()
{
  fluid_.Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::NonlinearSolve()
{
  fluid_.NonlinearSolve();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FSI::FluidAdapter::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FSI::FluidAdapter::PressureRowMap()
{
  return fluid_.PressureRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::FluidAdapter::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  meshmap_.Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::FluidAdapter::ResidualScaling() const
{
  return fluid_.ResidualScaling();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::ReadRestart(int step)
{
  fluid_.ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::FluidAdapter::Time()
{
  return fluid_.Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int FSI::FluidAdapter::Step()
{
  return fluid_.Step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::LiftDrag()
{
  fluid_.LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAdapter::ExtractInterfaceForces()
{
  return interface_.ExtractCondVector(fluid_.TrueResidual());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  interface_.InsertCondVector(ivel,fluid_.Velnp());

  // mark all interface velocities as dirichlet values
  // this is very easy, but there are two dangers:
  // - We change ivel here. It must not be used afterwards.
  // - The algorithm must support the sudden change of dirichtoggle_
  ivel->PutScalar(1.0);
  interface_.InsertCondVector(ivel,fluid_.Dirichlet());

  //----------------------- compute an inverse of the dirichtoggle vector
  fluid_.InvDirichlet()->PutScalar(1.0);
  fluid_.InvDirichlet()->Update(-1.0,*fluid_.Dirichlet(),1.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp) const
{
  meshmap_.InsertCondVector(fluiddisp,fluid_.Dispnp());

  // new grid velocity
  fluid_.UpdateGridv();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel) const
{
  meshmap_.InsertCondVector(gridvel,fluid_.GridVel());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int FSI::FluidAdapter::Itemax() const
{
  return fluid_.Itemax();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidAdapter::SetItemax(int itemax)
{
  fluid_.SetItemax(itemax);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAdapter::IntegrateInterfaceShape()
{
  return interface_.ExtractCondVector(fluid_.IntegrateInterfaceShape("FSICoupling"));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidAdapter::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
{
  const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
  Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap,true);
  interface_.InsertCondVector(ivel,relax);
  fluid_.LinearRelaxationSolve(relax);
  return ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> FSI::FluidAdapter::CreateFieldTest()
{
  return Teuchos::rcp(new FluidResultTest(fluid_));
}
//----------------------------------------------------------------------
//----------------------------------------------------------------------
//----------------------------------------------------------------------
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidGenAlphaAdapter::FluidGenAlphaAdapter(
  Teuchos::RCP<DRT::Discretization>      dis,
  Teuchos::RCP<LINALG::Solver>           solver,
  Teuchos::RCP<ParameterList>            params,
  Teuchos::RCP<IO::DiscretizationWriter> output,
  bool                                   isale)
  : fluid_ (dis, *solver, *params, *output, isale),
    dis_   (dis),
    solver_(solver),
    params_(params),
    output_(output)
{
  FSI::UTILS::SetupInterfaceExtractor(*dis,"FSICoupling",interface_);

  // build inner velocity map
  // dofs at the interface are excluded
  // we use only velocity dofs and only those without Dirichlet constraint

  Teuchos::RCP<const Epetra_Map>    velmap       = fluid_.VelocityRowMap();
  Teuchos::RCP<const Epetra_Vector> dirichtoggle = fluid_.Dirichlet();
  Teuchos::RCP<const Epetra_Map>    fullmap      = DofRowMap();

  int numvelids = velmap->NumMyElements();
  std::vector<int> velids;
  velids.reserve(numvelids);
  for (int i=0; i<numvelids; ++i)
  {
    int gid = velmap->GID(i);
    if (not interface_.CondMap()->MyGID(gid) and (*dirichtoggle)[fullmap->LID(gid)]==0.)
    {
      velids.push_back(gid);
    }
  }

  innervelmap_ = Teuchos::rcp(new Epetra_Map(-1,velids.size(), &velids[0], 0, velmap->Comm()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FSI::FluidGenAlphaAdapter::InitialGuess() const
{
  return fluid_.InitialGuess();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FSI::FluidGenAlphaAdapter::RHS() const
{
  return fluid_.Residual();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FSI::FluidGenAlphaAdapter::Velnp() const
{
  return fluid_.Velnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> FSI::FluidGenAlphaAdapter::Veln() const
{
  return fluid_.Veln();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FSI::FluidGenAlphaAdapter::DofRowMap() const
{
  const Epetra_Map* dofrowmap = dis_->DofRowMap();
  return Teuchos::rcp(dofrowmap, false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> FSI::FluidGenAlphaAdapter::SystemMatrix() const
{
  return fluid_.SysMat();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::Discretization> FSI::FluidGenAlphaAdapter::Discretization()
{
  return fluid_.Discretization();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidGenAlphaAdapter::PrepareTimeStep()
{
  fluid_.GenAlphaIncreaseTimeAndStep();

  fluid_.GenAlphaEchoToScreen("print time algorithm info");
  fluid_.GenAlphaPredictNewSolutionValues();
  fluid_.GenAlphaApplyDirichletAndNeumann();
  fluid_.GenAlphaCalcInitialAccelerations();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidGenAlphaAdapter::Evaluate(Teuchos::RCP<const Epetra_Vector> vel) const
{
  // Yes, this is complicated. But we have to be very careful
  // here. The field solver always expects an increment only. And
  // there are Dirichlet conditions that need to be preserved. So take
  // the sum of increments we get from NOX and apply the latest
  // increment only.
  if (vel!=Teuchos::null)
  {
    Teuchos::RCP<Epetra_Vector> incvel = Teuchos::rcp(new Epetra_Vector(*vel));
    incvel->Update(-1.0,*fluid_.Velnp(),1.0);
    fluid_.ExternIncrementOfVelnp(incvel);
  }

  fluid_.GenAlphaComputeIntermediateSol();

  fluid_.GenAlphaAssembleResidualAndMatrix();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidGenAlphaAdapter::Update()
{
  fluid_.GenAlphaTimeUpdate();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidGenAlphaAdapter::Output()
{
  fluid_.GenAlphaOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidGenAlphaAdapter::NonlinearSolve()
{
  fluid_.DoGenAlphaPredictorCorrectorIteration();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FSI::FluidGenAlphaAdapter::InnerVelocityRowMap()
{
  return innervelmap_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> FSI::FluidGenAlphaAdapter::PressureRowMap()
{
  return fluid_.PressureRowMap();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::FluidGenAlphaAdapter::SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
{
  meshmap_.Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::FluidGenAlphaAdapter::ResidualScaling() const
{
  return fluid_.ResidualScaling();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidGenAlphaAdapter::ReadRestart(int step)
{
  fluid_.ReadRestart(step);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
double FSI::FluidGenAlphaAdapter::Time()
{
  return fluid_.Time();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int FSI::FluidGenAlphaAdapter::Step()
{
  return fluid_.Step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidGenAlphaAdapter::LiftDrag()
{
  fluid_.LiftDrag();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidGenAlphaAdapter::ExtractInterfaceForces()
{
  return interface_.ExtractCondVector(fluid_.TrueResidual());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidGenAlphaAdapter::ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
{
  interface_.InsertCondVector(ivel,fluid_.Velnp());

  // mark all interface velocities as dirichlet values
  // this is very easy, but there are two dangers:
  // - We change ivel here. It must not be used afterwards.
  // - The algorithm must support the sudden change of dirichtoggle_
  ivel->PutScalar(1.0);
  interface_.InsertCondVector(ivel,fluid_.Dirichlet());

  //----------------------- compute an inverse of the dirichtoggle vector
  fluid_.InvDirichlet()->PutScalar(1.0);
  fluid_.InvDirichlet()->Update(-1.0,*fluid_.Dirichlet(),1.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidGenAlphaAdapter::ApplyMeshDisplacement(Teuchos::RCP<Epetra_Vector> fluiddisp) const
{
  meshmap_.InsertCondVector(fluiddisp,fluid_.Dispnp());

  // new grid velocity
  fluid_.UpdateGridv();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidGenAlphaAdapter::ApplyMeshVelocity(Teuchos::RCP<Epetra_Vector> gridvel) const
{
  meshmap_.InsertCondVector(gridvel,fluid_.GridVel());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int FSI::FluidGenAlphaAdapter::Itemax() const
{
  return fluid_.Itemax();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidGenAlphaAdapter::SetItemax(int itemax)
{
  fluid_.SetItemax(itemax);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidGenAlphaAdapter::IntegrateInterfaceShape()
{
  return interface_.ExtractCondVector(fluid_.IntegrateInterfaceShape("FSICoupling"));
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::FluidGenAlphaAdapter::RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
{
  const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
  Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap,true);
  interface_.InsertCondVector(ivel,relax);
  fluid_.LinearRelaxationSolve(relax);
  return ExtractInterfaceForces();
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<DRT::ResultTest> FSI::FluidGenAlphaAdapter::CreateFieldTest()
{
  return Teuchos::rcp(new FluidResultTest(fluid_));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidBaseAlgorithm::FluidBaseAlgorithm(const Teuchos::ParameterList& prbdyn, bool isale)
{
  SetupFluid(prbdyn, isale);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::FluidBaseAlgorithm::~FluidBaseAlgorithm()
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::FluidBaseAlgorithm::SetupFluid(const Teuchos::ParameterList& prbdyn, bool& isale)
{
  Teuchos::RCP<Teuchos::Time> t = Teuchos::TimeMonitor::getNewTimer("FSI::FluidBaseAlgorithm::SetupFluid");
  Teuchos::TimeMonitor monitor(*t);

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

  if ((actdis->Comm()).MyPID()==0)
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

  // -----------------------sublist containing stabilization parameters
  fluidtimeparams->sublist("STABILIZATION")=fdyn.sublist("STABILIZATION");

  // --------------------------sublist containing turbulence parameters
  {
    fluidtimeparams->sublist("TURBULENCE MODEL")=fdyn.sublist("TURBULENCE MODEL");

    fluidtimeparams->sublist("TURBULENCE MODEL").set<string>("statistics outfile",allfiles.outputfile_kenner);
  }

  if (iop == timeint_one_step_theta || iop == timeint_bdf2)
  {
    // number of degrees of freedom
    fluidtimeparams->set<int>              ("number of velocity degrees of freedom" ,probsize.get<int>("DIM"));
    // the default time step size
    fluidtimeparams->set<double>           ("time step size"           ,prbdyn.get<double>("TIMESTEP"));
    // max. sim. time
    fluidtimeparams->set<double>           ("total time"               ,prbdyn.get<double>("MAXTIME"));
    // parameter for time-integration
    fluidtimeparams->set<double>           ("theta"                    ,fdyn.get<double>("THETA"));
    // which kind of time-integration
    fluidtimeparams->set<FLUID_TIMEINTTYPE>("time int algo"            ,iop);
    // bound for the number of timesteps
    fluidtimeparams->set<int>              ("max number timesteps"     ,prbdyn.get<int>("NUMSTEP"));
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
    fluidtimeparams->set                 ("write restart every"       ,prbdyn.get<int>("RESTARTEVRY"));
    // solution output
    fluidtimeparams->set                 ("write solution every"      ,prbdyn.get<int>("UPRES"));
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
    fluid_ = rcp(new FluidAdapter(actdis, solver, fluidtimeparams, output, isale));
  }
  else if (iop == timeint_gen_alpha)
  {
    // -------------------------------------- number of degrees of freedom
    // number of degrees of freedom
    fluidtimeparams->set<int>              ("number of velocity degrees of freedom" ,probsize.get<int>("DIM"));

    // ------------------------------------------------ basic scheme, i.e.
    // --------------------- solving nonlinear or linearised flow equation
    fluidtimeparams->set<int>("type of nonlinear solve" ,
                             Teuchos::getIntegralValue<int>(fdyn,"DYNAMICTYP"));

    // -------------------------------------------------- time integration
    // the default time step size
    fluidtimeparams->set<double>           ("time step size"           ,prbdyn.get<double>("TIMESTEP"));
    // maximum simulation time
    fluidtimeparams->set<double>           ("total time"               ,prbdyn.get<double>("MAXTIME"));
    // maximum number of timesteps
    fluidtimeparams->set<int>              ("max number timesteps"     ,prbdyn.get<int>("NUMSTEP"));

    // ---------------------------------------------- nonlinear iteration
    // set linearisation scheme
    fluidtimeparams->set<bool>("Use reaction terms for linearisation",
                              Teuchos::getIntegralValue<int>(fdyn,"NONLINITER")==2);
    // maximum number of nonlinear iteration steps
    fluidtimeparams->set<int>             ("max nonlin iter steps"     ,fdyn.get<int>("ITEMAX"));
    // stop nonlinear iteration when both incr-norms are below this bound
    fluidtimeparams->set<double>          ("tolerance for nonlin iter" ,fdyn.get<double>("CONVTOL"));
    // set convergence check
    fluidtimeparams->set<string>          ("CONVCHECK"  ,fdyn.get<string>("CONVCHECK"));

    // ----------------------------------------------- restart and output
    // restart
    fluidtimeparams->set                  ("write restart every"       ,prbdyn.get<int>("RESTARTEVRY"));
    // solution output
    fluidtimeparams->set                  ("write solution every"      ,prbdyn.get<int>("UPRES"));
    // flag for writing stresses
    fluidtimeparams->set                  ("write stresses"            ,Teuchos::getIntegralValue<int>(ioflags,"FLUID_STRESS"));
    // ---------------------------------------------------- lift and drag
    fluidtimeparams->set<int>("liftdrag",Teuchos::getIntegralValue<int>(fdyn,"LIFTDRAG"));

    // -----------evaluate error for test flows with analytical solutions
    int init = Teuchos::getIntegralValue<int>(fdyn,"INITIALFIELD");
    fluidtimeparams->set                  ("eval err for analyt sol"   ,init);

    // -------------------------------------------------------------------
    // parameter alpha_M for for generalized-alpha scheme
    fluidtimeparams->set<double>           ("alpha_M"                  ,fdyn.get<double>("ALPHA_M"));
    // parameter alpha_F for for generalized-alpha scheme
    fluidtimeparams->set<double>           ("alpha_F"                  ,fdyn.get<double>("ALPHA_F"));

    //--------------------------------------------------
    // create all vectors and variables associated with the time
    // integration (call the constructor)
    // the only parameter from the list required here is the number of
    // velocity degrees of freedom
    fluid_ = rcp(new FluidGenAlphaAdapter(actdis, solver, fluidtimeparams, output, isale));
  }
  else
  {
    dserror("Unknown time integration for FSI fluid\n");
  }
  return;
}

#endif
