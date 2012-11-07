/*----------------------------------------------------------------------*/
/*!
\file fsi_monolithic.cpp

\brief General framework for monolithic fsi solution schemes

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-15262
</pre>
*/

/*----------------------------------------------------------------------*/

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <NOX_Epetra_Interface_Preconditioner.H>
#include <NOX_Direction_UserDefinedFactory.H>

#include "fsi_monolithic.H"
#include "fsi_debugwriter.H"
#include "fsi_nox_aitken.H"
#include "fsi_nox_group.H"
#include "fsi_nox_newton.H"
#include "fsi_statustest.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../linalg/linalg_blocksparsematrix.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "../drt_constraint/constraint_manager.H"

#include "../drt_io/io_control.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_ale/ale_utils_mapextractor.H"


/*----------------------------------------------------------------------*/
// Note: The order of calling the three BaseAlgorithm-constructors is
// important here! In here control file entries are written. And these
// entries define the order in which the filters handle the
// Discretizations, which in turn defines the dof number ordering of the
// Discretizations.
/*----------------------------------------------------------------------*/
FSI::MonolithicBase::MonolithicBase(const Epetra_Comm& comm,
                                    const Teuchos::ParameterList& timeparams)
  : AlgorithmBase(comm,timeparams)
{
  // ask base algorithm for the structural time integrator
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(timeparams));
  structure_ = rcp_dynamic_cast<ADAPTER::FSIStructureWrapper>(structure->StructureFieldrcp());

  if(structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");

  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(timeparams,true));
  fluid_ = fluid->FluidFieldrcp();

  // ask base algorithm for the ale time integrator
  Teuchos::RCP<ALE::AleBaseAlgorithm> ale = Teuchos::rcp(new ALE::AleBaseAlgorithm(timeparams));
  ale_ = ale->AleFieldrcp();

  coupsf_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupsa_ = Teuchos::rcp(new ADAPTER::Coupling());
  coupfa_ = Teuchos::rcp(new ADAPTER::Coupling());
  icoupfa_ = Teuchos::rcp(new ADAPTER::Coupling());
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
  StructureField()->ReadRestart(step);
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

  StructureField()->PrepareTimeStep();
  FluidField().    PrepareTimeStep();
  AleField().      PrepareTimeStep();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBase::Update()
{
  StructureField()->Update();
  FluidField().    Update();
  AleField().      Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBase::PrepareOutput()
{
  StructureField()->PrepareOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBase::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  StructureField()->Output();
  FluidField().    Output();
  AleField().      Output();

  FluidField().LiftDrag();

  if (StructureField()->GetConstraintManager()->HaveMonitor())
  {
    StructureField()->GetConstraintManager()->ComputeMonitorValues(StructureField()->Dispnp());
    if(Comm().MyPID() == 0)
      StructureField()->GetConstraintManager()->PrintMonitorValues();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::StructToAle(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsa_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::StructToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::FluidToStruct(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupsf_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::FluidToAleInterface(Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::StructToAle(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToStruct(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::StructToFluid(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_->MasterToSlave(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::FluidToStruct(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToFluid(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::FluidToAleInterface(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfa_->MasterToSlave(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Monolithic::Monolithic(const Epetra_Comm& comm,
                            const Teuchos::ParameterList& timeparams)
  : MonolithicBase(comm,timeparams)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  // enable debugging
  if (DRT::INPUT::IntegralValue<int>(fsidyn,"DEBUGOUTPUT")==1)
  {
    sdbg_ = Teuchos::rcp(new UTILS::DebugWriter(StructureField()->Discretization()));
    //fdbg_ = Teuchos::rcp(new UTILS::DebugWriter(FluidField().Discretization()));
  }

  std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
  s.append(".iteration");
  log_ = Teuchos::rcp(new std::ofstream(s.c_str()));

  firstcall_ = true;

  ddgpre_ = Teuchos::null;
  dugpre_ = Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::SetupSystem()
{
  // right now we use matching meshes at the interface

  const int ndim = DRT::Problem::Instance()->NDim();

  ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  ADAPTER::Coupling& coupsa = StructureAleCoupling();
  ADAPTER::Coupling& coupfa = FluidAleCoupling();
  ADAPTER::Coupling& icoupfa = InterfaceFluidAleCoupling();

  // structure to fluid

  coupsf.SetupConditionCoupling(*StructureField()->Discretization(),
                                 StructureField()->Interface()->FSICondMap(),
                                *FluidField().Discretization(),
                                 FluidField().Interface()->FSICondMap(),
                                "FSICoupling",
                                 ndim);

  // structure to ale

  coupsa.SetupConditionCoupling(*StructureField()->Discretization(),
                                 StructureField()->Interface()->FSICondMap(),
                                *AleField().Discretization(),
                                 AleField().Interface()->FSICondMap(),
                                "FSICoupling",
                                 ndim);

  // fluid to ale at the interface

  icoupfa.SetupConditionCoupling(*FluidField().Discretization(),
                                   FluidField().Interface()->FSICondMap(),
                                   *AleField().Discretization(),
                                   AleField().Interface()->FSICondMap(),
                                   "FSICoupling",
                                   ndim);

  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf.MasterDofMap()->SameAs(*coupsa.MasterDofMap()))
    dserror("structure interface dof maps do not match");

  if (coupsf.MasterDofMap()->NumGlobalElements()==0)
    dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField().Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap   = AleField().Discretization()->NodeRowMap();

  coupfa.SetupCoupling(*FluidField().Discretization(),
                       *AleField().Discretization(),
                       *fluidnodemap,
                       *alenodemap,
                        ndim);

  FluidField().SetMeshMap(coupfa.MasterDofMap());

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
  PrepareTimeloop();

  while (NotFinished())
  {
    PrepareTimeStep();
    TimeStep(interface);
    PrepareOutput();
    Update();
    Output();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::PrepareTimeloop()
{
  // make sure we didn't destroy the maps before we entered the timeloop
  Extractor().CheckForValidMapExtractor();

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = NOXParameterList();

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm().MyPID());

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

  // Create printing utilities
  utils_ = Teuchos::rcp(new NOX::Utils(printParams));

  if (Comm().MyPID()==0)
  {
    (*log_) << "# num procs      = " << Comm().NumProc() << "\n"
            << "# Method         = " << nlParams.sublist("Direction").get("Method","Newton") << "\n"
            << "# step | time | time/step | #nliter | res-norm | #liter\n"
            << "#\n"
      ;
  }

  // check for prestressing,
  // do not allow monolithic in the pre-phase
  // allow monolithic in the post-phase
  {
    double time = 0.0;
    double dt = 0.0;
    double pstime = -1.0;
    const ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn,"PRESTRESS");
    if (pstype != INPAR::STR::prestress_none)
    {
      time   = StructureField()->GetTime();
      dt     = StructureField()->GetTimeStepSize();
      pstime = sdyn.get<double>("PRESTRESSTIME");
      if (time+dt <= pstime) dserror("No monolithic FSI in the pre-phase of prestressing, use Aitken!");
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::TimeStep(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
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

  Teuchos::Time timer("time step timer");

  if (sdbg_!=Teuchos::null)
    sdbg_->NewTimeStep(Step(),"struct");
  if (fdbg_!=Teuchos::null)
    fdbg_->NewTimeStep(Step(),"fluid");

  // start time measurement
  Teuchos::RCP<Teuchos::TimeMonitor> timemonitor = rcp(new Teuchos::TimeMonitor(timer,true));

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
  noxstatus_ = solver->solve();

  if (noxstatus_ != NOX::StatusTest::Converged)
    dserror("Nonlinear solver failed to converge!");

  // recover Lagrange multiplier \lambda_\Gamma at the interface at the end of each time step
  // (i.e. condensed forces onto the structure) needed for rhs in next time step
  RecoverLagrangeMultiplier();

  // cleanup
  //mat_->Zero();

  // stop time measurement
  timemonitor = Teuchos::null;

  if (Comm().MyPID()==0)
  {
    (*log_) << Step()
            << "\t" << Time()
            << "\t" << timer.totalElapsedTime()
            << "\t" << nlParams.sublist("Output").get("Nonlinear Iterations",0)
            << "\t" << nlParams.sublist("Output").get("2-Norm of Residual", 0.)
            << "\t" << lsParams.sublist("Output").get("Total Number of Linear Iterations",0)
      ;
    (*log_) << std::endl;
    lsParams.sublist("Output").set("Total Number of Linear Iterations",0);
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
    ExtractFieldVectors(x,sx,fx,ax);

    if (sdbg_!=Teuchos::null)
    {
      sdbg_->NewIteration();
      sdbg_->WriteVector("x",*StructureField()->Interface()->ExtractFSICondVector(sx));
    }
  }

  // Call all elements and assemble rhs and matrices
  // We only need the rhs here because NOX will ask for the rhs
  // only. But the Jacobian is stored internally and will be returnd
  // later on without looking at x again!

  Utils()->out() << "\nEvaluate elements\n";

  {
    Epetra_Time ts(Comm());
    StructureField()->Evaluate(sx);
    Utils()->out() << "structure: " << ts.ElapsedTime() << " sec\n";
  }

  {
    Epetra_Time ta(Comm());
    AleField()      .Evaluate(ax);
    Utils()->out() << "ale      : " << ta.ElapsedTime() << " sec\n";
  }

  // transfer the current ale mesh positions to the fluid field
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluid(AleField().ExtractDisplacement());
  FluidField().ApplyMeshDisplacement(fluiddisp);

  {
    Epetra_Time tf(Comm());
    FluidField().Evaluate(fx);
    Utils()->out() << "fluid    : " << tf.ElapsedTime() << " sec\n";
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



  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  dirParams.set<std::string>("Method","User Defined");
  Teuchos::RCP<NOX::Direction::UserDefinedFactory> newtonfactory = Teuchos::rcp(this,false);
  dirParams.set("User Defined Direction Factory",newtonfactory);


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
  lsParams.set<int>("Max Iterations",1000);
  lsParams.set<std::string>("Preconditioner","User Defined");
  lsParams.set<int>("Output Frequency",10);
  lsParams.set<bool>("Output Solver Details",true);

  // adaptive tolerance settings
  lsParams.set<double>("base tolerance",fsidyn.get<double>("BASETOL"));
  lsParams.set<double>("adaptive distance",fsidyn.get<double>("ADAPTIVEDIST"));

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
FSI::BlockMonolithic::BlockMonolithic(const Epetra_Comm& comm,
                                      const Teuchos::ParameterList& timeparams)
  : Monolithic(comm,timeparams),
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

  if (pcdbg_!=Teuchos::null)
  {
    pcdbg_->NewLinearSystem();
  }

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

