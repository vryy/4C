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
#include "fsi_nox_group.H"
#include "fsi_nox_newton.H"
#include "fsi_statustest.H"

#include "../drt_inpar/inpar_fsi.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_blocksparsematrix.H"
#include "../linalg/linalg_utils.H"

#include "../drt_ale/ale.H"

#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include "../drt_constraint/constraint_manager.H"

#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"

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
  // access the structural discretization
  Teuchos::RCP<DRT::Discretization> structdis = DRT::Problem::Instance()->GetDis("structure");

  // access structural dynamic params list which will be possibly modified while creating the time integrator
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // ask base algorithm for the structural time integrator
  Teuchos::RCP<ADAPTER::StructureBaseAlgorithm> structure = Teuchos::rcp(new ADAPTER::StructureBaseAlgorithm(timeparams, const_cast<Teuchos::ParameterList&>(sdyn), structdis));
  structure_ = Teuchos::rcp_dynamic_cast<ADAPTER::FSIStructureWrapper>(structure->StructureFieldrcp());

  if(structure_ == Teuchos::null)
    dserror("cast from ADAPTER::Structure to ADAPTER::FSIStructureWrapper failed");

  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(timeparams,DRT::Problem::Instance()->FluidDynamicParams(),true));
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
  FluidField().ReadRestart(step);
  AleField().ReadRestart(step);

  SetTimeStep(FluidField().Time(), FluidField().Step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBase::PrepareTimeStep()
{
  IncrementTimeAndStep();

  PrintHeader();

  StructureField()->PrepareTimeStep();
  FluidField().PrepareTimeStep();
  AleField().PrepareTimeStep();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBase::Update()
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  bool timeadapton = DRT::INPUT::IntegralValue<int>(fsidyn.sublist("TIMEADAPTIVITY"),"TIMEADAPTON");

  if (not timeadapton)
    StructureField()->Update(); // constant dt
  else
    StructureField()->Update(Time()); // variable/adaptive dt

  FluidField().Update();
  AleField().Update();
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
  FluidField().Output();
  AleField().Output();

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
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToFluidInterface(Teuchos::RCP<Epetra_Vector> iv) const
{
  return icoupfa_->SlaveToMaster(iv);
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
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBase::AleToFluidInterface(Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfa_->SlaveToMaster(iv);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Monolithic::Monolithic(const Epetra_Comm& comm,
                            const Teuchos::ParameterList& timeparams)
  : MonolithicBase(comm,timeparams),
    firstcall_(true),
    noxiter_(0),
    log_(Teuchos::null),
    logada_(Teuchos::null)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

  // enable debugging
  if (DRT::INPUT::IntegralValue<int>(fsidyn, "DEBUGOUTPUT") == 1)
  {
    sdbg_ = Teuchos::rcp(new UTILS::DebugWriter(StructureField()->Discretization()));
    //fdbg_ = Teuchos::rcp(new UTILS::DebugWriter(FluidField().Discretization()));
  }

  // write iterations-file
  std::string fileiter = DRT::Problem::Instance()->OutputControlFile()->FileName();
  fileiter.append(".iteration");
  log_ = Teuchos::rcp(new std::ofstream(fileiter.c_str()));

  // "Initialize" interface solution increments due to structural predictor
  ddgpred_ = Teuchos::null;

  erroraction_ = FSI::Monolithic::erroraction_stop;
  numhalvestep_ = 0;

  //-------------------------------------------------------------------------
  // time step size adaptivity
  //-------------------------------------------------------------------------
  const bool timeadapton = DRT::INPUT::IntegralValue<bool>(fsidyn.sublist("TIMEADAPTIVITY"),"TIMEADAPTON");

  if (timeadapton)
  {
    InitTimIntAda(fsidyn);
  }
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
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  const bool timeadapton = DRT::INPUT::IntegralValue<bool>(fsidyn.sublist("TIMEADAPTIVITY"),"TIMEADAPTON");

  // Run time loop with constant or adaptive time step size (depending on the user's will)
  if (not timeadapton)
  {
    // call time loop with constant time step size
    TimeloopConstDt(interface);
  }
  else
  {
    // call time loop with adaptive time step size
    TimeloopAdaDt(interface);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::TimeloopConstDt(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
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

  // write header of log-file
  if (Comm().MyPID() == 0)
  {
    (*log_) << "# num procs      = " << Comm().NumProc() << "\n"
            << "# Method         = " << nlParams.sublist("Direction").get<std::string>("Method") << "\n"
            << "# step | time | time/step  |  #nliter  | res-norm |   #liter  |   dt  | ";

    (*log_) << "#\n\n";
  }

  WriteAdaFileHeader();

  // check for prestressing,
  // do not allow monolithic in the pre-phase
  // allow monolithic in the post-phase
  {
    const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
    INPAR::STR::PreStress pstype = DRT::INPUT::IntegralValue<INPAR::STR::PreStress>(sdyn,"PRESTRESS");
    if (pstype != INPAR::STR::prestress_none)
    {
      const double pstime = sdyn.get<double>("PRESTRESSTIME");

      if (Time() + Dt() <= pstime)
        dserror("No monolithic FSI in the pre-phase of prestressing, use Aitken!");
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

  if (sdbg_ != Teuchos::null)
    sdbg_->NewTimeStep(Step(), "struct");
  if (fdbg_ != Teuchos::null)
    fdbg_->NewTimeStep(Step(), "fluid");

  // Single field predictors have been applied, so store the structural
  // interface displacement increment due to predictor or inhomogeneous
  // Dirichlet boundary conditions
  ddgpred_ = Teuchos::rcp(new Epetra_Vector(*StructureField()->ExtractInterfaceDispnp()));
  ddgpred_->Update(-1.0, *StructureField()->ExtractInterfaceDispn(), 1.0);

  // start time measurement
  Teuchos::RCP<Teuchos::TimeMonitor> timemonitor = Teuchos::rcp(new Teuchos::TimeMonitor(timer,true));

  // calculate initial linear system at current position (no increment)
  // This initializes the field algorithms and creates the first linear
  // systems. And this is the reason we know the initial linear system is
  // there when we create the NOX::FSI::Group.
  Evaluate(Teuchos::null);

  // Get initial guess.
  // The initial system is there, so we can happily extract the
  // initial guess. (The Dirichlet conditions are already build in!)
  Teuchos::RCP<Epetra_Vector> initial_guess = Teuchos::rcp(new Epetra_Vector(*DofRowMap(),true));
  InitialGuess(initial_guess);

  NOX::Epetra::Vector noxSoln(initial_guess, NOX::Epetra::Vector::CreateView);

  // Create the linear system
  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys = CreateLinearSystem(nlParams, noxSoln, utils_);

  // Create the Group
  Teuchos::RCP<NOX::FSI::Group> grp =
    Teuchos::rcp(new NOX::FSI::Group(*this, printParams, interface, noxSoln, linSys));

  // Convergence Tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo = CreateStatusTest(nlParams, grp);

  // Create the solver
  Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp,combo,Teuchos::RCP<Teuchos::ParameterList>(&nlParams,false));

  // we know we already have the first linear system calculated
  grp->CaptureSystemState();

  // solve the whole thing
  noxstatus_ = solver->solve();
  noxiter_ = solver->getNumIterations();

  // recover Lagrange multiplier \lambda_{\Gamma} at the interface at the end of each time step
  // (i.e. condensed traction/forces onto the structure) needed for rhs in next time step
  RecoverLagrangeMultiplier();

  // stop time measurement
  timemonitor = Teuchos::null;

  if (Comm().MyPID()==0)
  {
    (*log_) << Step()
            << "\t" << Time()
            << "\t" << timer.totalElapsedTime()
            << "\t" << nlParams.sublist("Output").get<int>("Nonlinear Iterations")
            << "\t" << nlParams.sublist("Output").get<double>("2-Norm of Residual")
            << "\t" << lsParams.sublist("Output").get<int>("Total Number of Linear Iterations")
            << "\t" << Dt();
  }
  (*log_) << std::endl;
  lsParams.sublist("Output").set("Total Number of Linear Iterations",0);

  // perform the error check to determine the error action to be performed
  NonLinErrorCheck();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::NonLinErrorCheck()
{
  // if everything is fine, then return right now
  if (NoxStatus() == NOX::StatusTest::Converged)
  {
    erroraction_ = FSI::Monolithic::erroraction_none;
    return;
  }

  // The nonlinear solver did not converge. Thus, we have to take some action
  // that depends on the user's will given in the input file

  // get the FSI parameter list
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

  // get the user's will
  const INPAR::FSI::DivContAct divcontype
    = DRT::INPUT::IntegralValue<INPAR::FSI::DivContAct>(fsidyn.sublist("TIMEADAPTIVITY"),("DIVERCONT"));

  if (NoxStatus() != NOX::StatusTest::Converged)
  {
    switch (divcontype)
    {
      case INPAR::FSI::divcont_stop:
      {
        // set the corresponding error action
        erroraction_ = FSI::Monolithic::erroraction_stop;

        // stop the simulation
        dserror("Nonlinear solver did not converge in %i iterations.", noxiter_);
        break;
      }
      case INPAR::FSI::divcont_continue:
      {
        // set the corresponding error action
        erroraction_ = FSI::Monolithic::erroraction_continue;

        // Notify user about non-converged nonlinear solver, but do not abort the simulation
        if (Comm().MyPID() == 0)
        {
          IO::cout << "\n*** Nonlinear solver did not converge in " << noxiter_ << " iterations. Continue ...\n";
        }
        break;
      }
      case INPAR::FSI::divcont_halve_step:
      {
        // set the corresponding error action
        erroraction_ = FSI::Monolithic::erroraction_halve_step;

        numhalvestep_++;

        const bool timeadapton = DRT::INPUT::IntegralValue<bool>(fsidyn.sublist("TIMEADAPTIVITY"),"TIMEADAPTON");
        if (not timeadapton)
        {
          dserror("Nonlinear solver wants to halve the time step size. This is not possible in a time integrator with constant Delta t.");
        }

        // Notify user about non-converged nonlinear solver, but do not abort the simulation
        if (Comm().MyPID() == 0)
        {
          IO::cout << "\n*** Nonlinear solver did not converge in " << noxiter_ << " iterations. Halve the time step size.\n";
        }
        break;
      }
      default:
      {
        dserror("Unknown action to cope with non-converged nonlinear solver.");
        break;
      }
    }
  }
  else
  {
    dserror("Unknown NOX::StatusTest::StatusType.");
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::Evaluate(Teuchos::RCP<const Epetra_Vector> x)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::Evaluate");

  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  Teuchos::RCP<const Epetra_Vector> ax;

  if (x != Teuchos::null)
  {
    ExtractFieldVectors(x,sx,fx,ax);

    if (sdbg_!=Teuchos::null)
    {
      sdbg_->NewIteration();
      sdbg_->WriteVector("x", *StructureField()->Interface()->ExtractFSICondVector(sx));
    }
  }

  // Call all elements and assemble rhs and matrices
  // We only need the rhs here because NOX will ask for the rhs
  // only. But the Jacobian is stored internally and will be returned
  // later on without looking at x again!

  Utils()->out() << "\nEvaluate elements\n";

  {
    Epetra_Time ts(Comm());
    StructureField()->Evaluate(sx);
    Utils()->out() << "structure: " << ts.ElapsedTime() << " sec\n";
  }

  {
    Epetra_Time ta(Comm());
    AleField().Evaluate(ax);
    Utils()->out() << "ale      : " << ta.ElapsedTime() << " sec\n";
  }

  // transfer the current ale mesh positions to the fluid field
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluid(AleField().Dispnp());
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
  blockrowdofmap_.Setup(*fullmap, maps);
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

  // ToDo: Remove the CONVTOL-tolerances and replace them by the single field
  //       tolerances given just below also for lung- and constraint-FSI
  //
  // Currently, we have to keep the parameter CONVTOL for lung and constraint FSI
  nlParams.set("Norm abs pres", fsidyn.get<double>("CONVTOL"));
  nlParams.set("Norm abs vel",  fsidyn.get<double>("CONVTOL"));
  nlParams.set("Norm abs disp", fsidyn.get<double>("CONVTOL"));

  // set tolerances for nonlinear solver
  nlParams.set("Tol dis res L2",  fsidyn.get<double>("TOL_DIS_RES_L2"));
  nlParams.set("Tol dis res Inf", fsidyn.get<double>("TOL_DIS_RES_INF"));
  nlParams.set("Tol dis inc L2",  fsidyn.get<double>("TOL_DIS_INC_L2"));
  nlParams.set("Tol dis inc Inf", fsidyn.get<double>("TOL_DIS_INC_INF"));
  nlParams.set("Tol fsi res L2",  fsidyn.get<double>("TOL_FSI_RES_L2"));
  nlParams.set("Tol fsi res Inf", fsidyn.get<double>("TOL_FSI_RES_INF"));
  nlParams.set("Tol fsi inc L2",  fsidyn.get<double>("TOL_FSI_INC_L2"));
  nlParams.set("Tol fsi inc Inf", fsidyn.get<double>("TOL_FSI_INC_INF"));
  nlParams.set("Tol pre res L2",  fsidyn.get<double>("TOL_PRE_RES_L2"));
  nlParams.set("Tol pre res Inf", fsidyn.get<double>("TOL_PRE_RES_INF"));
  nlParams.set("Tol pre inc L2",  fsidyn.get<double>("TOL_PRE_INC_L2"));
  nlParams.set("Tol pre inc Inf", fsidyn.get<double>("TOL_PRE_INC_INF"));
  nlParams.set("Tol vel res L2",  fsidyn.get<double>("TOL_VEL_RES_L2"));
  nlParams.set("Tol vel res Inf", fsidyn.get<double>("TOL_VEL_RES_INF"));
  nlParams.set("Tol vel inc L2",  fsidyn.get<double>("TOL_VEL_INC_L2"));
  nlParams.set("Tol vel inc Inf", fsidyn.get<double>("TOL_VEL_INC_INF"));

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

  // adaptive tolerance settings for linear solver
  lsParams.set<double>("base tolerance",fsidyn.get<double>("BASETOL")); // relative tolerance
  lsParams.set<double>("adaptive distance",fsidyn.get<double>("ADAPTIVEDIST")); // adaptive distance
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Direction::Generic>
FSI::Monolithic::buildDirection(const Teuchos::RCP<NOX::GlobalData>& gd,
                                Teuchos::ParameterList& params) const
{
  Teuchos::RCP<NOX::FSI::Newton> newton = Teuchos::rcp(new NOX::FSI::Newton(gd, params));
  for (unsigned i = 0; i < statustests_.size(); ++i)
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
void FSI::Monolithic::SetupRHS(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::SetupRHS");

  firstcall_ = firstcall;

  // We want to add into a zero vector
  f.PutScalar(0.0);

  // contributions of single field residuals
  SetupRHSResidual(f);

  // contributions of Lagrange multiplier from last time step
  SetupRHSLambda(f);

  // contributions of special "first nonlinear iteration" terms
  if (firstcall)
    SetupRHSFirstiter(f);

  if (dbcmaps_ != Teuchos::null) // ToDo: Remove if-"Abfrage" after 'dbcmaps_' has been introduced in lung fsi
  {
    // Finally, we take care of Dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(f));
    Teuchos::RCP<const Epetra_Vector> zeros = Teuchos::rcp(new const Epetra_Vector(f.Map(), true));
    LINALG::ApplyDirichlettoSystem(rhs, zeros, *(dbcmaps_->CondMap()));
    f.Update(1.0, *rhs, 0.0);
  }

  // NOX expects the 'positive' residual. The negative sign for the
  // linearized Newton system J*dx=-r is done internally by NOX.
  // Since we assembled the right hand side, we have to invert the sign here.
  f.Scale(-1.);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::InitialGuess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::Monolithic::InitialGuess");

  CombineFieldVectors(*ig,
                      StructureField()->InitialGuess(),
                      FluidField().InitialGuess(),
                      AleField().InitialGuess(),
                      true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::CombineFieldVectors(Epetra_Vector& v,
                                          Teuchos::RCP<const Epetra_Vector> sv,
                                          Teuchos::RCP<const Epetra_Vector> fv,
                                          Teuchos::RCP<const Epetra_Vector> av)
{
  Extractor().AddVector(*sv, 0, v);
  Extractor().AddVector(*fv, 1, v);
  Extractor().AddVector(*av, 2, v);
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

  if (precondreusecount_ <= 0)
  {
    // Create preconditioner operator. The blocks are already there. This is
    // the perfect place to initialize the block preconditioners.
    SystemMatrix()->SetupPreconditioner();

    const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
    precondreusecount_ = fsidyn.get<int>("PRECONDREUSE");
  }

  precondreusecount_ -= 1;

  if (pcdbg_ != Teuchos::null)
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

