
#ifdef CCADISCRET

#include "fsi_monolithic.H"
#include "fsi_nox_group.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"

#include "../drt_lib/drt_colors.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

extern "C"
{
#include "../headers/standardtypes.h"
}


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
// Note: The order of calling the three BaseAlgorithm-constructors is
// important here! In here control file entries are written. And these
// entries define the order in which the filters handle the
// Discretizations, which in turn defines the dof number ordering of the
// Discretizations.
/*----------------------------------------------------------------------*/
FSI::MonolithicBase::MonolithicBase(Epetra_Comm& comm)
  : StructureBaseAlgorithm(),
    FluidBaseAlgorithm(DRT::Problem::Instance()->FSIDynamicParams(),true),
    AleBaseAlgorithm(),
    comm_(comm)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  if (comm_.MyPID()==0)
    DRT::INPUT::PrintDefaultParameters(std::cout, fsidyn);

  step_ = 0;
  time_ = 0.;
  dt_ = fsidyn.get<double>("TIMESTEP");
  nstep_ = fsidyn.get<int>("NUMSTEP");
  maxtime_ = fsidyn.get<double>("MAXTIME");
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
  FluidField().ReadRestart(step);
  AleField().ReadRestart(step);

  time_ = FluidField().Time();
  step_ = FluidField().Step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBase::PrepareTimeStep()
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
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = NOXParameterList();

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
    Teuchos::RCP<NOXGroup> grp =
      Teuchos::rcp(new NOXGroup(*this, printParams, interface, noxSoln, linSys));

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
void FSI::Monolithic::Evaluate(Teuchos::RCP<const Epetra_Vector> x) const
{
  if (x!=Teuchos::null)
  {
    double norm;
    int err = x->Norm2(&norm);
    if (err)
      dserror("failed to calculate norm");

    if (norm!=0)
    {
      Utils()->out() << YELLOW_LIGHT "element call with new x" END_COLOR << endl;

      Teuchos::RCP<const Epetra_Vector> sx;
      Teuchos::RCP<const Epetra_Vector> fx;
      Teuchos::RCP<const Epetra_Vector> ax;

      ExtractFieldVectors(x,sx,fx,ax);

      // debug
      //debug_.DumpVector("sx",*StructureField()->Discretization(),*sx);
      //debug_.DumpVector("fx",*FluidField()->Discretization(),*fx);
      //debug_.DumpVector("ax",*AleField()->Discretization(),*ax);

      // Call all elements and assemble rhs and matrices
      // We only need the rhs here because NOX will ask for the rhs
      // only. But the Jacobian is stored internally and will be returnd
      // later on without looking at x again!
      StructureField().Evaluate(sx);
      AleField()      .Evaluate(ax);

      // transfer the current ale mesh positions to the fluid field
      Teuchos::RefCountPtr<Epetra_Vector> fluiddisp = AleToFluid(AleField().ExtractDisplacement());
      FluidField().ApplyMeshDisplacement(fluiddisp);

      FluidField().Evaluate(fx);
    }
  }
  else
  {
    Utils()->out() << YELLOW_LIGHT "element call at current x" END_COLOR << endl;

    StructureField().Evaluate(Teuchos::null);
    AleField()      .Evaluate(Teuchos::null);

    // transfer the current ale mesh positions to the fluid field
    Teuchos::RefCountPtr<Epetra_Vector> fluiddisp = AleToFluid(AleField().ExtractDisplacement());
    FluidField().ApplyMeshDisplacement(fluiddisp);

    FluidField().Evaluate(Teuchos::null);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Monolithic::SetDofRowMaps(const std::vector<Teuchos::RCP<const Epetra_Map> >& maps)
{
  Teuchos::RCP<Epetra_Map> fullmap = LINALG::MultiMapExtractor::MergeMaps(maps);
  blockrowdofmap_.Setup(*fullmap,maps);
}


#endif
