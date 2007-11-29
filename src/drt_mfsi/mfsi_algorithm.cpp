
#ifdef CCADISCRET

#include "mfsi_algorithm.H"
#include "mfsi_nox_thyra_group.H"
#include "mfsi_preconditionerfactory.H"
#include "mfsi_statustest.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_colors.H"

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
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;

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
  FSI_DYNAMIC *fsidyn = alldyn[3].fsidyn;
  step_ = 0;
  time_ = 0.;
  dt_ = fsidyn->dt;
  nstep_ = fsidyn->nstep;
  maxtime_ = fsidyn->maxtime;

  SetupStructure();
  SetupFluid();
  SetupAle();

  // right now we use matching meshes at the interface

  // structure to fluid

  coupsf_.SetupConditionCoupling(*structure_->Discretization(),
                                 *fluid_->Discretization(),
                                 "FSICoupling");

  // setup the very simple structure to fluid coupling
  // u(n+1)*dt = d(n+1) - d(n)
  // but we solve both fields for the increments
  // the structure is even solved for the middle increment
  //
  // du(n+1)*dt = 1/(1-alpha_f)*dd(n+m)

  coupsf_.SetupCouplingMatrices(*structure_->Discretization()->DofRowMap(),
                                *fluid_->Discretization()->DofRowMap());

//   coupsf_.MasterToMasterMat()->Scale(structure_->DispIncrFactor());
//   coupsf_.MasterToMasterMatTrans()->Scale(structure_->DispIncrFactor());
  coupsf_.SlaveToMasterMat()->Scale(-dt_);
  coupsf_.SlaveToMasterMatTrans()->Scale(-dt_);

  // structure to ale

  coupsa_.SetupConditionCoupling(*structure_->Discretization(),
                                 *ale_->Discretization(),
                                 "FSICoupling");

  // setup structure to ale coupling
  //
  // dd(G,n+1) = 1/(1-alpha_f)*dd(n+m)

  coupsa_.SetupCouplingMatrices(*structure_->Discretization()->DofRowMap(),
                                *ale_->Discretization()->DofRowMap());

//   coupsa_.MasterToMasterMat()->Scale(structure_->DispIncrFactor());
//   coupsa_.MasterToMasterMatTrans()->Scale(structure_->DispIncrFactor());
  coupsa_.SlaveToMasterMat()->Scale(-1.);
  coupsa_.SlaveToMasterMatTrans()->Scale(-1.);

  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf_.MasterDofMap()->SameAs(*coupsa_.MasterDofMap()))
    dserror("structure interface dof maps do not match");

  if (coupsf_.MasterDofMap()->NumGlobalElements()==0)
    dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");

  // init transfer from interface to field
  structure_->SetInterfaceMap(coupsf_.MasterDofMap());
  fluid_    ->SetInterfaceMap(coupsf_.SlaveDofMap());
  ale_      ->SetInterfaceMap(coupsa_.SlaveDofMap());

  ifstruct_ = Teuchos::rcp(new Epetra_Vector(*structure_->InterfaceMap()));
  iastruct_ = Teuchos::rcp(new Epetra_Vector(*structure_->InterfaceMap()));

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = fluid_->Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap   = ale_->Discretization()->NodeRowMap();

  coupfa_.SetupCoupling(*fluid_->Discretization(),
                        *ale_->Discretization(),
                        *fluidnodemap,
                        *alenodemap);

  fluid_->SetMeshMap(coupfa_.MasterDofMap());

  // create Thyra vector spaces from Epetra maps for all fields
  smap_ = Thyra::create_VectorSpace(structure_->DofRowMap());
  fmap_ = Thyra::create_VectorSpace(fluid_->DofRowMap());
  amap_ = Thyra::create_VectorSpace(ale_->DofRowMap());

  sfmap_ = Thyra::create_VectorSpace(coupsf_.MasterDofMap());

  // stack vector spaces to build a composite that contains them all
  int numBlocks = 5;
  std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double> > > vecSpaces(numBlocks);
  vecSpaces[0] = smap_;
  vecSpaces[1] = fmap_;
  vecSpaces[2] = amap_;
  vecSpaces[3] = sfmap_;
  vecSpaces[4] = sfmap_;
  dofrowmap_ = Teuchos::rcp(new Thyra::DefaultProductVectorSpace<double>(numBlocks, &vecSpaces[0]));

#if 0
  Teuchos::RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  // use a solver builder for the standard field solvers so we get all the default behaviour
  Thyra::DefaultRealLinearSolverBuilder linearSolverBuilder;

  if (comm_.MyPID()==0)
  {
    std::cout << *linearSolverBuilder.getValidParameters() << std::endl;

    //linearSolverBuilder.getValidParameters()->print(std::cout,
    //                                                Teuchos::ParameterList::PrintOptions().showDoc(true).indent(2).showTypes(true));
  }

  Teuchos::RCP< Teuchos::ParameterList > paramList = Teuchos::rcp(new ParameterList());
  linearSolverBuilder.setParameterList(paramList);

  if (comm_.MyPID()==0)
  {
    paramList->print(std::cout);
  }

  structsolverfactory_ = linearSolverBuilder.createLinearSolveStrategy("");
  structsolverfactory_->setOStream(out);
  structsolverfactory_->setVerbLevel(Teuchos::VERB_LOW);

  fluidsolverfactory_ = linearSolverBuilder.createLinearSolveStrategy("");
  fluidsolverfactory_->setOStream(out);
  fluidsolverfactory_->setVerbLevel(Teuchos::VERB_LOW);

  alesolverfactory_ = linearSolverBuilder.createLinearSolveStrategy("");
  alesolverfactory_->setOStream(out);
  alesolverfactory_->setVerbLevel(Teuchos::VERB_LOW);
#else
  // field solvers used within the block preconditioner
  structsolverfactory_ = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory(Thyra::Amesos::KLU));
  fluidsolverfactory_ = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory(Thyra::Amesos::UMFPACK));
  alesolverfactory_ = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory(Thyra::Amesos::KLU));
#endif

  sfidentity_ = Teuchos::rcp(new Thyra::DefaultIdentityLinearOp<double>(Thyra::create_VectorSpace(structure_->InterfaceMap())));

  //Thyra::ConstLinearOperator<double> sfihandle = sfidentity;

  // the factory to create the special block preconditioner
  preconditionerfactory_ =
    Teuchos::rcp(new PreconditionerFactory(structsolverfactory_,
                                           fluidsolverfactory_,
                                           alesolverfactory_,
                                           sfidentity_));

  // Lets use aztec for now. This a about the only choice we have got.
  solverfactory_ = Teuchos::rcp(new Thyra::AztecOOLinearOpWithSolveFactory());
  solverfactory_->setPreconditionerFactory(preconditionerfactory_, "FSI block preconditioner");

  // We cannot use Amesos, since it expects the unterlying matrix to
  // be a Epetra_Operator.
  //solverfactory_ = Teuchos::rcp(new Thyra::AmesosLinearOpWithSolveFactory());

  // Create the system matrix. Right now it is empty. It is filled in
  // create_W_op().
  mat_ = Teuchos::rcp(new Thyra::DefaultBlockedLinearOp<double>());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::SetupStructure()
{
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
  FSI_DYNAMIC*    fsidyn   = alldyn[3].fsidyn;
  STRUCT_DYNAMIC* sdyn     = alldyn[genprob.numsf].sdyn;

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

  genalphaparams->set<bool>  ("damping",sdyn->damp);
  genalphaparams->set<double>("damping factor K",sdyn->k_damp);
  genalphaparams->set<double>("damping factor M",sdyn->m_damp);

  genalphaparams->set<double>("beta",sdyn->beta);
  genalphaparams->set<double>("gamma",sdyn->gamma);
  genalphaparams->set<double>("alpha m",sdyn->alpha_m);
  genalphaparams->set<double>("alpha f",sdyn->alpha_f);

  genalphaparams->set<double>("total time",0.0);
  genalphaparams->set<double>("delta time",fsidyn->dt);
  genalphaparams->set<int>   ("step",0);
  genalphaparams->set<int>   ("nstep",fsidyn->nstep);
  genalphaparams->set<int>   ("max iterations",sdyn->maxiter);
  genalphaparams->set<int>   ("num iterations",-1);
  genalphaparams->set<double>("tolerance displacements",sdyn->toldisp);

  genalphaparams->set<bool>  ("io structural disp",ioflags.struct_disp);
  genalphaparams->set<int>   ("io disp every nstep",sdyn->updevry_disp);
  genalphaparams->set<bool>  ("io structural stress",ioflags.struct_stress);
  genalphaparams->set<int>   ("io stress every nstep",sdyn->updevry_stress);
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

  FSI_DYNAMIC *fsidyn     = alldyn[3].fsidyn;
  FLUID_DYNAMIC *fdyn     = alldyn[genprob.numff].fdyn;
  fdyn->step              =   0;
  fdyn->acttime           = 0.0;

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

  // number of degrees of freedom
  fluidtimeparams->set<int>              ("number of velocity degrees of freedom" ,genprob.ndim);
  // the default time step size
  fluidtimeparams->set<double>           ("time step size"           ,fsidyn->dt);
  // max. sim. time
  fluidtimeparams->set<double>           ("total time"               ,fsidyn->maxtime);
  // parameter for time-integration
  fluidtimeparams->set<double>           ("theta"                    ,fdyn->theta);
  // which kind of time-integration
  fluidtimeparams->set<FLUID_TIMEINTTYPE>("time int algo"            ,fdyn->iop);
  // bound for the number of timesteps
  fluidtimeparams->set<int>              ("max number timesteps"     ,fsidyn->nstep);
  // number of steps with start algorithm
  fluidtimeparams->set<int>              ("number of start steps"    ,fdyn->nums);
  // parameter for start algo
  fluidtimeparams->set<double>           ("start theta"              ,fdyn->thetas);


  // ---------------------------------------------- nonlinear iteration
  // maximum number of nonlinear iteration steps
  fluidtimeparams->set<int>             ("max nonlin iter steps"     ,fdyn->itemax);
  // stop nonlinear iteration when both incr-norms are below this bound
  fluidtimeparams->set<double>          ("tolerance for nonlin iter" ,fdyn->ittol);

  // ----------------------------------------------- restart and output
  // restart
  fluidtimeparams->set                 ("write restart every"       ,fdyn->uprestart);
  // solution output
  fluidtimeparams->set                 ("write solution every"      ,fdyn->upres);
  // flag for writing stresses
  fluidtimeparams->set                 ("write stresses"            ,ioflags.fluid_stress);

  //--------------------------------------------------
  // evaluate error for test flows with analytical solutions
  fluidtimeparams->set                  ("eval err for analyt sol"   ,fdyn->init);


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

  FSI_DYNAMIC *fsidyn     = alldyn[3].fsidyn;
  ALE_DYNAMIC *adyn     = alldyn[genprob.numaf].adyn;
  adyn->step            =   0;
  adyn->time            = 0.0;

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  RefCountPtr<ParameterList> solveparams = rcp(new ParameterList());
  RefCountPtr<LINALG::Solver> solver =
    rcp(new LINALG::Solver(solveparams,actdis->Comm(),allfiles.out_err));
  solver->TranslateSolverParameters(*solveparams,actsolv);
  actdis->ComputeNullSpaceIfNecessary(*solveparams);

  RefCountPtr<ParameterList> params = rcp(new ParameterList());
  params->set<int>("nstep", fsidyn->nstep);
  params->set<double>("maxtime", fsidyn->maxtime);
  params->set<double>("dt", fsidyn->dt);

  ale_ = rcp(new FSI::AleLinear(actdis, solver, params, output));
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
  printParams.set("MyPID", comm_.MyPID());

  // turn on output
  printParams.set("Output Information", 0xffff);

  // Create printing utilities
  utils_ = Teuchos::rcp(new NOX::Utils(printParams));

  Teuchos::RefCountPtr<std::ofstream> log;
  if (comm_.MyPID()==0)
  {
    std::string s = allfiles.outputfile_kenner;
    s.append(".iteration");
    log = Teuchos::rcp(new std::ofstream(s.c_str()));
    (*log) << "# num procs      = " << comm_.NumProc() << "\n"
           << "# Method         = " << nlParams.sublist("Direction").get("Method","Newton") << "\n"
           << "#\n"
      ;
  }

  Teuchos::Time timer("time step timer");

  while (step_ < nstep_ and time_ <= maxtime_)
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
      if (comm_.MyPID()==0)
        utils_->out() << RED "Nonlinear solver failed to converge!" END_COLOR << endl;

    // cleanup
    // Do not keep the block matrices. They are too heavy! And the
    // next Evaluate() call will replace them anyway.
    mat_->uninitialize();

    // stop time measurement
    timemonitor = Teuchos::null;

    if (comm_.MyPID()==0)
    {
      (*log) << step_
             << " " << timer.totalElapsedTime()
             << " " << nlParams.sublist("Output").get("Nonlinear Iterations",0)
             << " " << nlParams.sublist("Output").get("2-Norm of Residual", 0.)
             << " " << lsParams.sublist("Output").get("Total Number of Linear Iterations",0)
        ;
      //for (std::vector<int>::size_type i=0; i<counter_.size(); ++i)
      //{
      //  (*log) << " " << counter_[i];
      //}
      (*log) << std::endl;
    }

    Update();
    Output();
  }
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
Thyra::ModelEvaluatorBase::InArgs<double> MFSI::Algorithm::createInArgs() const
{
  // Here we define what kinds of input arguments
  // MFSI::Algorithm::evalModelImpl() supports
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
  // MFSI::Algorithm::evalModelImpl() supports
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
    SetupSysMat(*mat_);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::Evaluate(Teuchos::RCP<const Thyra::DefaultProductVector<double> > x) const
{
  if (x!=Teuchos::null)
  {
    if (Thyra::norm(*x)!=0)
    {
      utils_->out() << YELLOW_LIGHT "element call with new x" END_COLOR << endl;

      Teuchos::RCP<const Epetra_Vector> sx = Thyra::get_Epetra_Vector(*structure_->DofRowMap(), x->getVectorBlock(0));
      Teuchos::RCP<const Epetra_Vector> fx = Thyra::get_Epetra_Vector(*fluid_    ->DofRowMap(), x->getVectorBlock(1));
      Teuchos::RCP<const Epetra_Vector> ax = Thyra::get_Epetra_Vector(*ale_      ->DofRowMap(), x->getVectorBlock(2));

//       debug.DumpVector("sx",*structure_->Discretization(),*sx);
//       debug.DumpVector("fx",*fluid_->Discretization(),*fx);
//       debug.DumpVector("ax",*ale_->Discretization(),*ax);

      // Call all elements and assemble rhs and matrices
      // We only need the rhs here because NOX will ask for the rhs
      // only. But the Jacobian is stored internally and will be returnd
      // later on without looking at x again!
      structure_->Evaluate(sx);
      ale_      ->Evaluate(ax);

      // transfer the current ale mesh positions to the fluid field
      Teuchos::RefCountPtr<Epetra_Vector> fluiddisp = coupfa_.SlaveToMaster(ale_->ExtractDisplacement());
      fluid_->ApplyMeshDisplacement(fluiddisp);

      fluid_    ->Evaluate(fx);
    }
  }
  else
  {
    utils_->out() << YELLOW_LIGHT "element call at current x" END_COLOR << endl;

    structure_->Evaluate(Teuchos::null);
    ale_      ->Evaluate(Teuchos::null);

    // transfer the current ale mesh positions to the fluid field
    Teuchos::RefCountPtr<Epetra_Vector> fluiddisp = coupfa_.SlaveToMaster(ale_->ExtractDisplacement());
    fluid_->ApplyMeshDisplacement(fluiddisp);

    fluid_    ->Evaluate(Teuchos::null);
  }

  //debug.DumpVector("sres",*structure_->Discretization(),*structure_->RHS());
  //debug.DumpVector("fres",*fluid_->Discretization(),*fluid_->RHS());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
#if 0
void MFSI::Algorithm::CurrentX(Thyra::DefaultProductVector<double>& x)
{
  Teuchos::RCP< Thyra::VectorBase< double > > sx = Thyra::create_Vector(structure_->Dispm(), smap_);
  Teuchos::RCP< Thyra::VectorBase< double > > fx = Thyra::create_Vector(fluid_->Vel(), fmap_);
  Teuchos::RCP< Thyra::VectorBase< double > > ax = Thyra::create_Vector(ale_->Disp(), amap_);

  Teuchos::RCP< Thyra::VectorBase< double > > sfx =
    Thyra::create_Vector(ifstruct_, Thyra::create_VectorSpace(structure_->InterfaceMap()));

  Teuchos::RCP< Thyra::VectorBase< double > > sax =
    Thyra::create_Vector(iastruct_, Thyra::create_VectorSpace(structure_->InterfaceMap()));

  int numBlocks = 5;
  std::vector<Teuchos::RCP<Thyra::VectorBase<double> > > vec(numBlocks);
  vec[0] = sx;
  vec[1] = fx;
  vec[2] = ax;
  vec[3] = sfx;
  vec[4] = sax;
  x.initialize(dofrowmap_,&vec[0]);
}
#endif


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::InitialGuess(Thyra::DefaultProductVector<double>& ig)
{
  // the linear field systems must be setup before the initial guess
  // is known
  Teuchos::RCP< const Thyra::VectorBase< double > > sig = Thyra::create_Vector(structure_->InitialGuess(), smap_);
  Teuchos::RCP< const Thyra::VectorBase< double > > fig = Thyra::create_Vector(fluid_->InitialGuess(), fmap_);
  Teuchos::RCP< const Thyra::VectorBase< double > > aig = Thyra::create_Vector(ale_->InitialGuess(), amap_);

  Teuchos::RCP< const Thyra::VectorBase< double > > sfig =
    Thyra::create_Vector(ifstruct_, Thyra::create_VectorSpace(structure_->InterfaceMap()));

  Teuchos::RCP< const Thyra::VectorBase< double > > saig =
    Thyra::create_Vector(iastruct_, Thyra::create_VectorSpace(structure_->InterfaceMap()));

  int numBlocks = 5;
  std::vector<Teuchos::RCP<const Thyra::VectorBase<double> > > vec(numBlocks);
  vec[0] = sig;
  vec[1] = fig;
  vec[2] = aig;
  vec[3] = sfig;
  vec[4] = saig;
  ig.initialize(dofrowmap_,&vec[0]);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::SetupRHS(Thyra::DefaultProductVector<double> &f) const
{
  // Extract RHS and put it into f

  // make a local copy so that we can modify the rhs vectors
  Teuchos::RCP<Epetra_Vector> sv = rcp(new Epetra_Vector(*structure_->RHS()));
  Teuchos::RCP<Epetra_Vector> fv = rcp(new Epetra_Vector(*fluid_->RHS()));
  Teuchos::RCP<Epetra_Vector> av = rcp(new Epetra_Vector(*ale_->RHS()));

//   debug.DumpVector("sf",*structure_->Discretization(),*sv);
//   debug.DumpVector("ff",*fluid_->Discretization(),*fv);
//   debug.DumpVector("af",*ale_->Discretization(),*av);

  // wrap epetra vectors in thyra
  Teuchos::RCP< Thyra::VectorBase< double > > srhs = Thyra::create_Vector(sv, smap_);
  Teuchos::RCP< Thyra::VectorBase< double > > frhs = Thyra::create_Vector(fv, fmap_);
  Teuchos::RCP< Thyra::VectorBase< double > > arhs = Thyra::create_Vector(av, amap_);

  // We couple absolute vectors, no increments. So we have a nonzero rhs.

  Teuchos::RCP< Epetra_Vector > ifstruct = structure_->FluidCondRHS();
  ifstruct->Update(1.0*dt_,*coupsf_.SlaveToMaster(fluid_->StructCondRHS()),-1.0);
  Teuchos::RCP< Thyra::VectorBase< double > > sfrhs =
    Thyra::create_Vector(ifstruct, Thyra::create_VectorSpace(structure_->InterfaceMap()));

  Teuchos::RCP< Epetra_Vector > iastruct = structure_->MeshCondRHS();
  iastruct->Update(1.0,*coupsa_.SlaveToMaster(ale_->StructCondRHS()),-1.0);
  Teuchos::RCP< Thyra::VectorBase< double > > sarhs =
    Thyra::create_Vector(iastruct, Thyra::create_VectorSpace(structure_->InterfaceMap()));

  // create block vector
  int numBlocks = 5;
  std::vector<Teuchos::RCP<Thyra::VectorBase<double> > > vec(numBlocks);
  vec[0] = srhs;
  vec[1] = frhs;
  vec[2] = arhs;
  vec[3] = sfrhs;
  vec[4] = sarhs;
  f.initialize(dofrowmap_,&vec[0]);

  // NOX expects a different sign here.
  Thyra::scale(-1., &f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void MFSI::Algorithm::SetupSysMat(Thyra::DefaultBlockedLinearOp<double>& mat) const
{
  // extract Jacobian matrices and put them into composite system
  // matrix W
  Teuchos::RCP<Thyra::LinearOpBase<double> > smat = Teuchos::rcp(new Thyra::EpetraLinearOp(structure_->SysMat()));
  Teuchos::RCP<Thyra::LinearOpBase<double> > fmat = Teuchos::rcp(new Thyra::EpetraLinearOp(fluid_->SysMat()));
  Teuchos::RCP<Thyra::LinearOpBase<double> > amat = Teuchos::rcp(new Thyra::EpetraLinearOp(ale_->SysMat()));

  mat.beginBlockFill(dofrowmap_, dofrowmap_);
  mat.setBlock(0,0,smat);
  mat.setBlock(1,1,fmat);
  mat.setBlock(2,2,amat);

  // structure to fluid coupling
  mat.setBlock(3,0,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsf_.MasterToMasterMat())));
  mat.setBlock(3,1,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsf_.SlaveToMasterMat())));
  mat.setBlock(0,3,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsf_.MasterToMasterMatTrans())));
  mat.setBlock(1,3,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsf_.SlaveToMasterMatTrans())));
  //mat.setBlock(3,3,sfidentity_);

  // structure to ale coupling
  // note there is no ale effect on the structural equations
  mat.setBlock(4,0,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsa_.MasterToMasterMat())));
  mat.setBlock(4,2,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsa_.SlaveToMasterMat())));
  //mat.setBlock(0,4,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsa_.MasterToMasterMatTrans())));
  mat.setBlock(2,4,Teuchos::rcp(new Thyra::EpetraLinearOp(coupsa_.SlaveToMasterMatTrans())));
  //mat.setBlock(4,4,sfidentity_);

  mat.endBlockFill();
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
Teuchos::RCP<NOX::StatusTest::Combo>
MFSI::Algorithm::CreateStatusTest(Teuchos::ParameterList& nlParams,
                                  Teuchos::RCP<NOX::Thyra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo       = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::StatusTest::Combo> converged   = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv    = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  Teuchos::RCP<PartialNormF> structureDisp =
    Teuchos::rcp(new PartialNormF("displacement",
                                  0,
                                  *structure_->DofRowMap(),
                                  *structure_->InnerDisplacementRowMap(),
                                  nlParams.get("Norm abs disp", 1.0e-6),
                                  PartialNormF::Scaled));
  converged->addStatusTest(structureDisp);

  Teuchos::RCP<PartialNormF> innerFluidVel =
    Teuchos::rcp(new PartialNormF("velocity",
                                  1,
                                  *fluid_->DofRowMap(),
                                  *fluid_->InnerVelocityRowMap(),
                                  nlParams.get("Norm abs vel", 1.0e-6),
                                  PartialNormF::Scaled));
  converged->addStatusTest(innerFluidVel);

  Teuchos::RCP<PartialNormF> fluidPress =
    Teuchos::rcp(new PartialNormF("pressure",
                                  1,
                                  *fluid_->DofRowMap(),
                                  *fluid_->PressureRowMap(),
                                  nlParams.get("Norm abs pres", 1.0e-6),
                                  PartialNormF::Scaled));
  converged->addStatusTest(fluidPress);

  Teuchos::RCP<InterfaceNormF> interface =
    Teuchos::rcp(new InterfaceNormF(1.,
                                    *structure_->DofRowMap(),
                                    *coupsf_.MasterDofMap(),
                                    fluid_->ResidualScaling(),
                                    *fluid_->DofRowMap(),
                                    *coupsf_.SlaveDofMap(),
                                    coupsf_,
                                    nlParams.get("Norm abs interface", 1.0e-6),
                                    PartialNormF::Scaled));
  converged->addStatusTest(interface);

#if 0
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
#endif

  //Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms     = Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  //converged->addStatusTest(wrms);

  return combo;
}


#endif
