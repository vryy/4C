
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE

#include "fsi_dirichletneumann.H"
#include "fsi_utils.H"
#include "../drt_lib/drt_globalproblem.H"

#include "../drt_lib/drt_colors.H"

#include "fsi_nox_aitken.H"
#include "fsi_nox_extrapolate.H"
#include "fsi_nox_michler.H"
#include "fsi_nox_fixpoint.H"
#include "fsi_nox_jacobian.H"
#include "fsi_nox_sd.H"


/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | vector of numfld FIELDs, defined in global_control.c                 |
 *----------------------------------------------------------------------*/
extern struct _FIELD      *field;

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

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | pointer to allocate dynamic variables if needed                      |
 | dedfined in global_control.c                                         |
 | ALLDYNA               *alldyn;                                       |
 *----------------------------------------------------------------------*/
extern ALLDYNA      *alldyn;

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
extern Teuchos::RefCountPtr<Teuchos::ParameterList> globalparameterlist;


FSI::DirichletNeumannCoupling::DirichletNeumannCoupling(Epetra_Comm& comm)
  : comm_(comm),
    counter_(7)
{
  FSI_DYNAMIC *fsidyn = alldyn[3].fsidyn;
  step_ = 0;
  time_ = 0.;
  dt_ = fsidyn->dt;
  nstep_ = fsidyn->nstep;
  maxtime_ = fsidyn->maxtime;

  displacementcoupling_ = globalparameterlist->get("Displacement Coupling", true);

  SetupStructure();
  SetupFluid();
  SetupAle();

  //cout << structure_->Discretization();
  //cout << fluid_->Discretization();
  //cout << ale_->Discretization();

  coupsf_.SetupConditionCoupling(structure_->Discretization(),
                                 fluid_->Discretization(),
                                 "FSICoupling");

  coupsa_.SetupConditionCoupling(structure_->Discretization(),
                                 ale_->Discretization(),
                                 "FSICoupling");

  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf_.MasterDofMap()->SameAs(*coupsa_.MasterDofMap()))
    dserror("structure interface dof maps do not match");

  // init transfer from interface to field
  structure_->SetInterfaceMap(coupsf_.MasterDofMap());
  fluid_    ->SetInterfaceMap(coupsf_.SlaveDofMap());
  ale_      ->SetInterfaceMap(coupsa_.SlaveDofMap());

  const Epetra_Map* fluidnodemap = fluid_->Discretization().NodeRowMap();
  const Epetra_Map* alenodemap   = ale_->Discretization().NodeRowMap();

  coupfa_.SetupCoupling(fluid_->Discretization(),
                        ale_->Discretization(),
                        *fluidnodemap,
                        *alenodemap);

  fluid_->SetMeshMap(coupfa_.MasterDofMap());
}


void FSI::DirichletNeumannCoupling::SetupStructure()
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
  Structure::SetDefaults(*genalphaparams);

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
  genalphaparams->set<int>   ("io disp every nstep",sdyn->updevry_stress);
  genalphaparams->set<bool>  ("print to screen",true);
  genalphaparams->set<bool>  ("print to err",true);
  genalphaparams->set<FILE*> ("err file",allfiles.out_err);

  // takes values "full newton" , "modified newton" , "nonlinear cg"
  genalphaparams->set<string>("equilibrium iteration","full newton");

  structure_ = rcp(new Structure(genalphaparams,actdis,solver,output));
}


void FSI::DirichletNeumannCoupling::SetupFluid()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
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
  Fluid::SetDefaults(*fluidtimeparams);

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

  // restart
  fluidtimeparams->set                  ("write restart every"       ,fdyn->uprestart);

  //--------------------------------------------------
  // evaluate error for test flows with analytical solutions
  fluidtimeparams->set                  ("eval err for analyt sol"   ,fdyn->init);


  //--------------------------------------------------
  // create all vectors and variables associated with the time
  // integration (call the constructor)
  // the only parameter from the list required here is the number of
  // velocity degrees of freedom
  fluid_ = rcp(new Fluid(actdis, solver, fluidtimeparams, output));
}


void FSI::DirichletNeumannCoupling::SetupAle()
{
  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  RefCountPtr<DRT::Discretization> actdis = null;
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

  ale_ = rcp(new AleLinear(actdis, solver, params, output));
}


void FSI::DirichletNeumannCoupling::Timeloop(const Teuchos::RefCountPtr<NOX::Epetra::Interface::Required>& interface)
{
  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = *globalparameterlist;

  // sublists

  //Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method","Newton"));
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  printParams.set("MyPID", comm_.MyPID());

  // Create printing utilities
  Teuchos::RefCountPtr<NOX::Utils> utils = Teuchos::rcp(new NOX::Utils(printParams));

  // Set user defined aitken line search object.
  if (nlParams.sublist("Line Search").get("Method","Aitken")=="Aitken")
  {
    // insert user defined aitken relaxation
    Teuchos::ParameterList& linesearch = nlParams.sublist("Line Search");
    Teuchos::RefCountPtr<NOX::LineSearch::Generic> aitken = Teuchos::rcp(new NOX::FSI::AitkenRelaxation(utils,linesearch));

    // We change the method here.
    linesearch.set("Method","User Defined");
    linesearch.set("User Defined Line Search",aitken);
  }

  // Set user defined steepest descent line search object.
  else if (nlParams.sublist("Line Search").get("Method","Aitken")=="SD")
  {
    // insert user defined aitken relaxation
    Teuchos::ParameterList& linesearch = nlParams.sublist("Line Search");
    Teuchos::RefCountPtr<NOX::LineSearch::Generic> sd = Teuchos::rcp(new NOX::FSI::SDRelaxation(utils,linesearch));

    // We change the method here.
    linesearch.set("Method","User Defined");
    linesearch.set("User Defined Line Search",sd);
  }

  // the very special experimental extrapolation
  else if (nlParams.sublist("Line Search").get("Method","Full Step")=="Extrapolate")
  {
    // insert user defined extrapolate relaxation
    Teuchos::ParameterList& linesearch = nlParams.sublist("Line Search");
    Teuchos::RefCountPtr<NOX::LineSearch::Generic> extrapolate = Teuchos::rcp(new NOX::FSI::Extrapolate(utils,linesearch));

    // We change the method here.
    linesearch.set("Method","User Defined");
    linesearch.set("User Defined Line Search",extrapolate);
  }

  // ==================================================================

  // log solver iterations

  Teuchos::RefCountPtr<std::ofstream> log;
  if (comm_.MyPID()==0)
  {
    std::string s = allfiles.outputfile_kenner;
    s.append(".iteration");
    log = Teuchos::rcp(new std::ofstream(s.c_str()));
    (*log) << "# num procs      = " << comm_.NumProc() << "\n"
           << "# Method         = " << dirParams.get("Method","Newton") << "\n"
           << "# Jacobian       = " << nlParams.get("Jacobian", "Matrix Free") << "\n"
           << "# Preconditioner = " << nlParams.get("Preconditioner","None") << "\n"
           << "# Line Search    = " << nlParams.sublist("Line Search").get("Method","Full Step") << "\n"
           << "#\n"
           << "# step  time/step  #nliter  #liter  Residual  Jac  Prec  FD_Res  MF_Res  MF_Jac  User\n"
      ;
  }

  idispn_ = InterfaceDisp();

  Teuchos::Time timer("time step timer");

  // ==================================================================

  while (step_ < nstep_ and time_ <= maxtime_)
  {
    PrepareTimeStep();

    // reset all counters
    std::fill(counter_.begin(),counter_.end(),0);
    lsParams.sublist("Output").set("Total Number of Linear Iterations",0);
    linsolvcount_.resize(0);

    // start time measurement
    Teuchos::RefCountPtr<Teuchos::TimeMonitor> timemonitor = rcp(new Teuchos::TimeMonitor(timer,true));

    /*----------------- CSD - predictor for itnum==0 --------------------*/

    // Begin Nonlinear Solver ************************************

    NOX::Epetra::Vector noxSoln(InterfaceDisp(), NOX::Epetra::Vector::CreateView);

    // Create the linear system
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Required> iReq = interface;

    // Ok. The variables.

    Teuchos::RefCountPtr<NOX::FSI::FSIMatrixFree> FSIMF;
    Teuchos::RefCountPtr<NOX::Epetra::MatrixFree> MF;
    Teuchos::RefCountPtr<NOX::Epetra::FiniteDifference> FD;
    Teuchos::RefCountPtr<NOX::Epetra::FiniteDifferenceColoring> FDC;
    Teuchos::RefCountPtr<NOX::Epetra::FiniteDifferenceColoring> FDC1;
    Teuchos::RefCountPtr<NOX::Epetra::BroydenOperator> B;

    Teuchos::RefCountPtr<NOX::Epetra::Interface::Jacobian> iJac;
    Teuchos::RefCountPtr<NOX::Epetra::Interface::Preconditioner> iPrec;

    Teuchos::RefCountPtr<Epetra_Operator> J;
    Teuchos::RefCountPtr<Epetra_Operator> M;

    Teuchos::RefCountPtr<NOX::Epetra::LinearSystemAztecOO> linSys;

    // ==================================================================
    // decide on Jacobian and preconditioner
    // We migh want to use no preconditioner at all. Some kind of
    // Jacobian has to be provided, otherwise the linear system uses
    // plain finite differences.

    std::string jacobian = nlParams.get("Jacobian", "None");
    std::string preconditioner = nlParams.get("Preconditioner", "None");

    // Special FSI based matrix free method
    if (jacobian=="FSI Matrix Free")
    {
      //Teuchos::ParameterList& mfParams = nlParams.sublist("FSI Matrix Free");

      // MatrixFree seems to be the most interessting choice. This
      // version builds on our steepest descent relaxation
      // implementation to approximate the Jacobian times x.

      // This is the default method.

      FSIMF = Teuchos::rcp(new NOX::FSI::FSIMatrixFree(printParams, interface, noxSoln));
      iJac = FSIMF;
      J = FSIMF;
    }

    // Matrix Free Newton Krylov. This is supposed to be the most
    // appropiate choice.
    else if (jacobian=="Matrix Free")
    {
      Teuchos::ParameterList& mfParams = nlParams.sublist("Matrix Free");
      double lambda = mfParams.get("lambda", 1.0e-6);
      //mfresitemax_ = mfParams.get("itemax", -1);

      // MatrixFree seems to be the most interessting choice. But you
      // must set a rather low tolerance for the linear solver.

      MF = Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface, noxSoln));
      MF->setLambda(lambda);
      iJac = MF;
      J = MF;
    }

    // No Jacobian at all. Do a fix point iteration. This is a user
    // extension, so we have to modify the parameter list here.
    else if (jacobian=="None")
    {
      Teuchos::RefCountPtr<NOX::Direction::Generic> fixpoint = Teuchos::rcp(new NOX::FSI::FixPoint(utils,nlParams));
      dirParams.set("Method","User Defined");
      dirParams.set("User Defined Direction",fixpoint);
      lsParams.set("Preconditioner","None");
      preconditioner="None";
    }

    // the strange coupling proposed by Michler
    else if (jacobian=="Michler")
    {
      Teuchos::RefCountPtr<NOX::Direction::Generic> michler = Teuchos::rcp(new NOX::FSI::Michler(utils,nlParams));
      dirParams.set("Method","User Defined");
      dirParams.set("User Defined Direction",michler);
      lsParams.set("Preconditioner","None");
      preconditioner="None";
    }

#if 0
    else if (jacobian=="Dumb Finite Difference")
    {
      Teuchos::ParameterList& fdParams = nlParams.sublist("Finite Difference");
      //fdresitemax_ = fdParams.get("itemax", -1);
      double alpha = fdParams.get("alpha", 1.0e-4);
      double beta  = fdParams.get("beta",  1.0e-6);
      std::string dt = fdParams.get("Difference Type","Forward");
      NOX::Epetra::FiniteDifference::DifferenceType dtype = NOX::Epetra::FiniteDifference::Forward;
      if (dt=="Forward")
        dtype = NOX::Epetra::FiniteDifference::Forward;
      else if (dt=="Backward")
        dtype = NOX::Epetra::FiniteDifference::Backward;
      else if (dt=="Centered")
        dtype = NOX::Epetra::FiniteDifference::Centered;
      else
        dserror("unsupported difference type '%s'",dt.c_str());

      FD = Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams, interface, noxSoln, rawGraph_, beta, alpha));
      FD->setDifferenceMethod(dtype);

      iJac = FD;
      J = FD;
    }
#endif

#if 0
    // Finite Difference with coloring. Build a Jacobian as good as it
    // gets. This is the closest we get to the true Jacobian.
    else if (jacobian=="Finite Difference")
    {
      Teuchos::ParameterList& fdParams = nlParams.sublist("Finite Difference");
      //fdresitemax_ = fdParams.get("itemax", -1);
      double alpha = fdParams.get("alpha", 1.0e-4);
      double beta  = fdParams.get("beta",  1.0e-6);
      std::string dt = fdParams.get("Difference Type","Forward");
      NOX::Epetra::FiniteDifferenceColoring::DifferenceType dtype = NOX::Epetra::FiniteDifferenceColoring::Forward;
      if (dt=="Forward")
        dtype = NOX::Epetra::FiniteDifferenceColoring::Forward;
      else if (dt=="Backward")
        dtype = NOX::Epetra::FiniteDifferenceColoring::Backward;
      else if (dt=="Centered")
        dtype = NOX::Epetra::FiniteDifferenceColoring::Centered;
      else
        dserror("unsupported difference type '%s'",dt.c_str());

      FDC = Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, colorMap, columns, true, false, beta, alpha));
      FDC->setDifferenceMethod(dtype);

      iJac = FDC;
      J = FDC;
    }
#endif

#if 0
    // Finite Difference with distance 1 coloring. Build just a
    // diagonal Jacobian.
    else if (jacobian=="Finite Difference 1")
    {
      Teuchos::ParameterList& fdParams = nlParams.sublist("Finite Difference");
      double alpha = fdParams.get("alpha", 1.0e-4);
      double beta  = fdParams.get("beta",  1.0e-6);
      std::string dt = fdParams.get("Difference Type","Forward");
      NOX::Epetra::FiniteDifferenceColoring::DifferenceType dtype = NOX::Epetra::FiniteDifferenceColoring::Forward;
      if (dt=="Forward")
        dtype = NOX::Epetra::FiniteDifferenceColoring::Forward;
      else if (dt=="Backward")
        dtype = NOX::Epetra::FiniteDifferenceColoring::Backward;
      else if (dt=="Centered")
        dtype = NOX::Epetra::FiniteDifferenceColoring::Centered;
      else
        dserror("unsupported difference type '%s'",dt.c_str());

      FDC1 = Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, distance1ColorMap, distance1Columns, true, true, beta, alpha));
      FDC1->setDifferenceMethod(dtype);

      iJac = FDC1;
      J = FDC1;
    }
#endif

    // Broyden. A strange idea that starts with a real Jacobian and
    // does modify it along the steps of the nonlinear solve.
    else if (jacobian=="Broyden")
    {
#if 0
      // we need a Jacobian to start with.
      if (is_null(FDC))
      {
        FDC = Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, colorMap, columns, true));
      }
      FDC->computeJacobian(*soln());

      Teuchos::RefCountPtr<Epetra_CrsMatrix> mat = Teuchos::rcp(new Epetra_CrsMatrix(FDC->getUnderlyingMatrix()));
      //B = Teuchos::rcp(new NOX::Epetra::BroydenOperator(nlParams, *soln(), mat));
      B = Teuchos::rcp(new NOX::Epetra::BroydenOperator(nlParams, utils, *soln(), mat));

      iJac = B;
      J = B;
#else
      dserror("Broyden Operator still unfinished");
#endif
    }
    else
    {
      dserror("unsupported Jacobian '%s'",jacobian.c_str());
    }

    // ==================================================================

    // No preconditioning at all. This might work. But on large
    // systems it probably won't.
    if (preconditioner=="None")
    {
      if (lsParams.get("Preconditioner", "None")!="None")
      {
        if (comm_.MyPID()==0)
          utils->out() << RED "Warning: Preconditioner turned on in linear solver settings.\n"
                       << "Jacobian operator will be used for preconditioning as well." END_COLOR "\n";
      }

      if (Teuchos::is_null(iJac))
      {
        // if no Jacobian has been set this better be the fix point
        // method.
        if (dirParams.get("Method","Newton")!="User Defined")
        {
          if (comm_.MyPID()==0)
            utils->out() << RED "Warning: No Jacobian for solver " << dirParams.get("Method","Newton") << END_COLOR "\n";
        }
        linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iReq, noxSoln));
      }
      else
      {
        linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iReq, iJac, J, noxSoln));
      }
    }

#if 0
    // Finite Difference with coloring. The best (and most expensive)
    // Jacobian we can get. It might do good on really hard problems.
    else if (preconditioner=="Finite Difference")
    {
      if (lsParams.get("Preconditioner", "None")=="None")
      {
        if (comm_.MyPID()==0)
          utils->out() << RED "Warning: Preconditioner turned off in linear solver settings." END_COLOR "\n";
      }

      // A real (approximated) matrix for preconditioning
#if 0
      if (is_null(FDC))
      {
        // FiniteDifferenceColoring might be a good preconditioner for
        // really hard problems. But you should construct it only once
        // a time step. Probably.
        FDC = Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, colorMap, columns, true));
      }
      iPrec = FDC;
      M = FDC;
#else
      Teuchos::RefCountPtr<NOX::Epetra::FiniteDifferenceColoring> precFDC =
        Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, colorMap, columns, true));
      iPrec = precFDC;
      M = precFDC;
#endif

      linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iJac, J, iPrec, M, noxSoln));
    }
#endif

#if 0
    // Finite Difference with distance 1 coloring. Supposed to be a
    // suitable (and cheap) preconditioner
    else if (preconditioner=="Finite Difference 1")
    {
      if (lsParams.get("Preconditioner", "None")=="None")
      {
        if (comm_.MyPID()==0)
          utils->out() << RED "Warning: Preconditioner truned off in linear solver settings." END_COLOR "\n";
      }
#if 0
      if (is_null(FDC1))
      {
        FDC1 = Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, distance1ColorMap, distance1Columns, true, true));
      }
      iPrec = FDC1;
      M = FDC1;
#else
      Teuchos::RefCountPtr<NOX::Epetra::FiniteDifferenceColoring> precFDC =
        Teuchos::rcp(new NOX::Epetra::FiniteDifferenceColoring(printParams, interface, noxSoln, rawGraph_, distance1ColorMap, distance1Columns, true, true));
      iPrec = precFDC;
      M = precFDC;
#endif

      linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iJac, J, iPrec, M, noxSoln));
    }
#endif

    else
    {
      dserror("unsupported preconditioner '%s'",preconditioner.c_str());
    }

    // ==================================================================
    // Convergence Tests

    // Create the Group
    Teuchos::RefCountPtr<NOX::Epetra::Group> grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, iReq, noxSoln, linSys));

#if 0
    // debug FD Jacobian
    grp->computeJacobian();
    if (!Teuchos::is_null(FD))
      EpetraExt::RowMatrixToMatlabFile("UnderlyingMatrix.fd",FD->getUnderlyingMatrix());
    else
      EpetraExt::RowMatrixToMatlabFile("UnderlyingMatrix.fdc",FDC->getUnderlyingMatrix());
    exit(1);
#endif

    // Create the convergence tests
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo       = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> converged   = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

    Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 20)));
    Teuchos::RefCountPtr<NOX::StatusTest::FiniteValue> fv    = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

    combo->addStatusTest(fv);
    combo->addStatusTest(converged);
    combo->addStatusTest(maxiters);

    Teuchos::RefCountPtr<NOX::StatusTest::NormF> absresid = Teuchos::rcp(new NOX::StatusTest::NormF(nlParams.get("Norm abs F", 1.0e-6)));
    converged->addStatusTest(absresid);

    if (nlParams.isParameter("Norm Update"))
    {
      Teuchos::RefCountPtr<NOX::StatusTest::NormUpdate> update = Teuchos::rcp(new NOX::StatusTest::NormUpdate(nlParams.get("Norm Update", 1.0e-5)));
      converged->addStatusTest(update);
    }

    if (nlParams.isParameter("Norm rel F"))
    {
      Teuchos::RefCountPtr<NOX::StatusTest::NormF> relresid = Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), nlParams.get("Norm rel F", 1.0e-2)));
      converged->addStatusTest(relresid);
    }

    //Teuchos::RefCountPtr<NOX::StatusTest::NormWRMS> wrms     = Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
    //converged->addStatusTest(wrms);

    //////////////////////////////////////////////////////////////////

    // Create the method
    solver_ = Teuchos::rcp(new NOX::Solver::Manager(grp, combo, globalparameterlist));
    NOX::StatusTest::StatusType status = solver_->solve();

    if (status != NOX::StatusTest::Converged)
      if (comm_.MyPID()==0)
        utils->out() << "Nonlinear solver failed to converge!" << endl;

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver_->getSolutionGroup());
    const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
    //const Epetra_Vector& finalF        = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getF())).getEpetraVector();

    idispn_->Update(1.0, finalSolution, 0.0);

    // End Nonlinear Solver **************************************

    // Output the parameter list
    if (utils->isPrintType(NOX::Utils::Parameters))
      if (comm_.MyPID()==0)
      {
        utils->out() << endl
                     << "Final Parameters" << endl
                     << "****************" << endl;
        solver_->getList().print(utils->out());
        utils->out() << endl;
      }

    // ==================================================================
    // return results

    if (displacementcoupling_)
    {
    }
    else
    {
      // do we need to distribute the new interface forces?
    }


    // ==================================================================

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
      for (std::vector<int>::size_type i=0; i<counter_.size(); ++i)
      {
        (*log) << " " << counter_[i];
      }
      (*log) << std::endl;
      log->flush();
    }

    // ==================================================================

    Update();

    /* write current solution */
    Output();
  }
}


Teuchos::RefCountPtr<Epetra_Vector> FSI::DirichletNeumannCoupling::InterfaceDisp()
{
  // extract displacements
  return structure_->ExtractInterfaceDisplacement();
}


bool FSI::DirichletNeumannCoupling::computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  char* flags[] = { "Residual", "Jac", "Prec", "FD_Res", "MF_Res", "MF_Jac", "User", NULL };

  if (comm_.MyPID()==0)
    cout << "\n==================================================================================================\n"
         << "FSI::DirichletNeumannCoupling::computeF: fillFlag = " RED << flags[fillFlag] << END_COLOR "\n\n";

  // we count the number of times the residuum is build
  counter_[fillFlag] += 1;

  if (!x.Map().UniqueGIDs())
    dserror("source map not unique");

  Teuchos::RefCountPtr<Epetra_Vector> idispn = rcp(new Epetra_Vector(x));

  Teuchos::RefCountPtr<Epetra_Vector> iforce = FluidOp(idispn, fillFlag);
  Teuchos::RefCountPtr<Epetra_Vector> idispnp = StructOp(iforce, fillFlag);

  F.Update(1.0, *idispnp, -1.0, *idispn, 0.0);

  return true;
}


Teuchos::RefCountPtr<Epetra_Vector>
FSI::DirichletNeumannCoupling::FluidOp(Teuchos::RefCountPtr<Epetra_Vector> idisp,
                                       const FillType fillFlag)
{
  if (fillFlag==User)
  {
    // SD relaxation calculation

    // Do we need to solve the ale here? Would the approximation
    // suffer otherwise?
    // Lets start with an unperturbed mesh
    //ale_->ApplyInterfaceDisplacements(StructToAle(idisp));
    //ale_->Solve();

    // the displacement -> velocity conversion at the interface
    Teuchos::RefCountPtr<Epetra_Vector> ivel = rcp(new Epetra_Vector(*idisp));
    ivel->Scale(1./dt_);

    //Teuchos::RefCountPtr<Epetra_Vector> fluiddisp = AleToFluid(ale_->ExtractDisplacement());

    return FluidToStruct(fluid_->RelaxationSolve(StructToFluid(ivel)));
  }
  else
  {
    // normal fluid solve

    ale_->ApplyInterfaceDisplacements(StructToAle(idisp));
    ale_->Solve();

    // the displacement -> velocity conversion at the interface
    Teuchos::RefCountPtr<Epetra_Vector> ivel = InterfaceVelocity(idisp);
    Teuchos::RefCountPtr<Epetra_Vector> fluiddisp = AleToFluid(ale_->ExtractDisplacement());

    fluid_->ApplyInterfaceVelocities(StructToFluid(ivel));
    fluid_->ApplyMeshDisplacement(fluiddisp);
    fluid_->NonlinearSolve();
    return FluidToStruct(fluid_->ExtractInterfaceForces());
  }
}


Teuchos::RefCountPtr<Epetra_Vector>
FSI::DirichletNeumannCoupling::StructOp(Teuchos::RefCountPtr<Epetra_Vector> iforce,
                                        const FillType fillFlag)
{
  if (fillFlag==User)
  {
    // SD relaxation calculation
    return structure_->RelaxationSolve(iforce);
  }
  else
  {
    // normal structure solve
    structure_->ApplyInterfaceForces(iforce);
    structure_->FullNewton();
    return structure_->ExtractInterfaceDisplacement();
  }
}


Teuchos::RefCountPtr<Epetra_Vector> FSI::DirichletNeumannCoupling::InterfaceVelocity(Teuchos::RefCountPtr<Epetra_Vector> idispnp)
{
  Teuchos::RefCountPtr<Epetra_Vector> ivel = rcp(new Epetra_Vector(*idispn_));
  ivel->Update(1.0, *idispnp, -1.0);
  ivel->Scale(1./dt_);
  return ivel;
}

Teuchos::RefCountPtr<Epetra_Vector> FSI::DirichletNeumannCoupling::StructToAle(Teuchos::RefCountPtr<Epetra_Vector> iv)
{
  return coupsa_.MasterToSlave(iv);
}


Teuchos::RefCountPtr<Epetra_Vector> FSI::DirichletNeumannCoupling::StructToFluid(Teuchos::RefCountPtr<Epetra_Vector> iv)
{
  return coupsf_.MasterToSlave(iv);
}


Teuchos::RefCountPtr<Epetra_Vector> FSI::DirichletNeumannCoupling::FluidToStruct(Teuchos::RefCountPtr<Epetra_Vector> iv)
{
  return coupsf_.SlaveToMaster(iv);
}


Teuchos::RefCountPtr<Epetra_Vector> FSI::DirichletNeumannCoupling::AleToFluid(Teuchos::RefCountPtr<Epetra_Vector> iv)
{
  return coupfa_.SlaveToMaster(iv);
}


void FSI::DirichletNeumannCoupling::PrepareTimeStep()
{
  step_ += 1;
  time_ += dt_;

  // update the fields
  //structure_->NewStep(step_, time_);
  //fluid_->NewStep(step_, time_);
  //ale_->NewStep(step_, time_);

  structure_->PrepareTimeStep();
  fluid_->    PrepareTimeStep();
  ale_->      PrepareTimeStep();
}

void FSI::DirichletNeumannCoupling::Update()
{
  structure_->UpdateandOutput();
  fluid_    ->TimeUpdate();
  ale_      ->Update();
}


void FSI::DirichletNeumannCoupling::Output()
{
  fluid_->Output();
  ale_  ->Output();
}


#endif
#endif
