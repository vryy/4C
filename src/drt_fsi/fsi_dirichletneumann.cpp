
#ifdef CCADISCRET

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
#include "fsi_nox_linearsystem_gcr.H"
#include "fsi_nox_mpe.H"
#include "fsi_nox_epsilon.H"

#include <string>
#include <Epetra_Time.h>


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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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

  if (globalparameterlist==Teuchos::null)
    globalparameterlist = Teuchos::rcp(new Teuchos::ParameterList);

  displacementcoupling_ = globalparameterlist->get("Displacement Coupling", true);

  SetupStructure();
  SetupFluid();
  SetupAle();

  //cout << structure_->Discretization();
  //cout << fluid_->Discretization();
  //cout << ale_->Discretization();

  std::string method = globalparameterlist->get("Coupling Method", "Matching Nodes");
  if (method=="Matching Nodes")
  {
    matchingnodes_ = true;
    coupsf_.SetupConditionCoupling(structure_->Discretization(),
                                   fluid_->Discretization(),
                                   "FSICoupling");

    coupsa_.SetupConditionCoupling(structure_->Discretization(),
                                   *ale_->Discretization(),
                                   "FSICoupling");

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
  }
  else if (method=="Mortar")
  {
    matchingnodes_ = false;
    coupsfm_.Setup( structure_->Discretization(),
                    fluid_->Discretization(),
                    comm );

    // This is cheating. We setup the coupling of interface dofs between fluid
    // and ale. But we use the variable from the matching version.
    coupsa_.SetupConditionCoupling(fluid_->Discretization(),
                                   *ale_->Discretization(),
                                   "FSICoupling");

    // init transfer from interface to field
    structure_->SetInterfaceMap(coupsfm_.MasterDofMap());
    fluid_    ->SetInterfaceMap(coupsfm_.SlaveDofMap());
    ale_      ->SetInterfaceMap(coupsa_.SlaveDofMap());
  }
  else
  {
    dserror("unsupported coupling method '%s'",method.c_str());
  }

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = fluid_->Discretization().NodeRowMap();
  const Epetra_Map* alenodemap   = ale_->Discretization()->NodeRowMap();

  coupfa_.SetupCoupling(fluid_->Discretization(),
                        *ale_->Discretization(),
                        *fluidnodemap,
                        *alenodemap);

  fluid_->SetMeshMap(coupfa_.MasterDofMap());

#if 0
  // create connection graph of interface elements
  Teuchos::RefCountPtr<Epetra_Map> imap = structure_->InterfaceMap();

  vector<int> rredundant;
  DRT::Utils::AllreduceEMap(rredundant, *imap);

  rawGraph_ = Teuchos::rcp(new Epetra_CrsGraph(Copy,*imap,12));
  for (int i=0; i<imap->NumMyElements(); ++i)
  {
    int err = rawGraph_->InsertGlobalIndices(imap->GID(i),rredundant.size(),&rredundant[0]);
    if (err < 0)
      dserror("Epetra_CrsGraph::InsertGlobalIndices returned %d", err);
  }
  rawGraph_->FillComplete();
#endif
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannCoupling::ReadRestart(int step)
{
  structure_->ReadRestart(step);
  fluid_->ReadRestart(step);
  ale_->ReadRestart(step);

  time_ = fluid_->time();
  step_ = fluid_->step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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
  genalphaparams->set<int>   ("io stress every nstep",sdyn->updevry_stress);

  genalphaparams->set<int>   ("restart",genprob.restart);
  genalphaparams->set<int>   ("write restart every",fsidyn->uprestart);

  genalphaparams->set<bool>  ("print to screen",true);
  genalphaparams->set<bool>  ("print to err",true);
  genalphaparams->set<FILE*> ("err file",allfiles.out_err);

  switch (sdyn->nlnSolvTyp)
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
  switch (sdyn->predtype)
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

  structure_ = rcp(new Structure(genalphaparams,actdis,solver,output));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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

  // ----------------------------------------------- restart and output
  // restart
  fluidtimeparams->set                 ("write restart every"       ,fsidyn->uprestart);
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
  fluid_ = rcp(new Fluid(actdis, solver, fluidtimeparams, output));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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

  FSI_DYNAMIC *fsidyn   = alldyn[3].fsidyn;
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

  // ----------------------------------------------- restart and output
  // restart
  params->set<int>("write restart every", fsidyn->uprestart);

  ale_ = rcp(new AleLinear(actdis, solver, params, output));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannCoupling::Timeloop(const Teuchos::RefCountPtr<NOX::Epetra::Interface::Required>& interface)
{
  bool secondsolver = false;

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = *globalparameterlist;

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method","Newton"));
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  //Teuchos::ParameterList& searchParams = nlParams.sublist("Line Search");
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", comm_.MyPID());

  // set default output flag to no output
  // The field solver will output a lot, anyway.
  printParams.get("Output Information",0);

  // Create printing utilities
  utils_ = Teuchos::rcp(new NOX::Utils(printParams));

  // Set user defined aitken line search object.
  if (nlParams.sublist("Line Search").get("Method","Aitken")=="Aitken")
  {
    // insert user defined aitken relaxation
    Teuchos::ParameterList& linesearch = nlParams.sublist("Line Search");
    Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> aitkenfactory =
      Teuchos::rcp(new NOX::FSI::AitkenFactory());

    // We change the method here.
    linesearch.set("Method","User Defined");
    linesearch.set("User Defined Line Search Factory", aitkenfactory);

    Teuchos::ParameterList& aitkenList = linesearch.sublist("Aitken");
    if (aitkenList.get("start steps only", false))
    {
      // we start with aitken and switch to a different solver
      // afterwards
      secondsolver = true;
    }
  }

  // Set user defined steepest descent line search object.
  else if (nlParams.sublist("Line Search").get("Method","Aitken")=="SD")
  {
    // insert user defined aitken relaxation
    Teuchos::ParameterList& linesearch = nlParams.sublist("Line Search");
    Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> sdfactory =
      Teuchos::rcp(new NOX::FSI::SDFactory());

    // We change the method here.
    linesearch.set("Method","User Defined");
    linesearch.set("User Defined Line Search Factory", sdfactory);
  }

#if 0
  // the very special experimental extrapolation
  else if (nlParams.sublist("Line Search").get("Method","Aitken")=="Extrapolate")
  {
    // insert user defined extrapolate relaxation
    Teuchos::ParameterList& linesearch = nlParams.sublist("Line Search");
    Teuchos::RefCountPtr<NOX::LineSearch::Generic> extrapolate = Teuchos::rcp(new NOX::FSI::Extrapolate(utils_,linesearch));

    // We change the method here.
    linesearch.set("Method","User Defined");
    linesearch.set("User Defined Line Search",extrapolate);
  }
#endif

  // ==================================================================

  // log solver iterations

  Teuchos::RefCountPtr<std::ofstream> log;
  if (comm_.MyPID()==0)
  {
    std::string s = allfiles.outputfile_kenner;
    s.append(".iteration");
    log = Teuchos::rcp(new std::ofstream(s.c_str()));
    (*log) << "# num procs      = " << comm_.NumProc() << "\n"
           << "# Method         = " << nlParams.sublist("Direction").get("Method","Newton") << "\n"
           << "# Jacobian       = " << nlParams.get("Jacobian", "None") << "\n"
           << "# Preconditioner = " << nlParams.get("Preconditioner","None") << "\n"
           << "# Line Search    = " << nlParams.sublist("Line Search").get("Method","Aitken") << "\n"
           << "#\n"
           << "# step  time/step  #nliter  |R|  #liter  Residual  Jac  Prec  FD_Res  MF_Res  MF_Jac  User\n"
      ;
  }

  idispn_ = InterfaceDisp();

  Teuchos::Time timer("time step timer");

  // ==================================================================

  while (step_ < nstep_ and time_ <= maxtime_)
  {
    // Increment all field counters and predict field values whenever
    // appropriate.
    PrepareTimeStep();

    // reset all counters
    std::fill(counter_.begin(),counter_.end(),0);
    lsParams.sublist("Output").set("Total Number of Linear Iterations",0);
    linsolvcount_.resize(0);

    // start time measurement
    Teuchos::RefCountPtr<Teuchos::TimeMonitor> timemonitor = rcp(new Teuchos::TimeMonitor(timer,true));

    /*----------------- CSD - predictor for itnum==0 --------------------*/

    // Begin Nonlinear Solver ************************************

    // Get initial guess.
    Teuchos::RefCountPtr<Epetra_Vector> soln;
    if (displacementcoupling_)
    {
      // Here some predictor could be useful.
      // On the other hand, the structural algorithm has its own
      // predictor and it has been run. We extract the predicted
      // interface displacements. Is there any sense in trying to
      // predict another time?

      // d(n)
      soln = InterfaceDisp();

      // These predictors are available within the old code.

      // d(n)+dt*(1.5*v(n)-0.5*v(n-1))
      // d(n)+dt*v(n)
      // d(n)+dt*v(n)+0.5*dt^2*a(n)

#if 0
      Teuchos::RefCountPtr<Epetra_Vector> veln = structure_->ExtractInterfaceVel();
      Teuchos::RefCountPtr<Epetra_Vector> accn = structure_->ExtractInterfaceAcc();
      soln->Update(dt_,*veln,0.5*dt_*dt_,*accn,1.);
#endif
    }
    else
    {
      soln = InterfaceForce();
    }

    NOX::Epetra::Vector noxSoln(soln, NOX::Epetra::Vector::CreateView);

    // Create the linear system
    Teuchos::RefCountPtr<NOX::Epetra::LinearSystem> linSys =
      CreateLinearSystem(nlParams, interface, noxSoln, utils_);

    // Create the Group
    Teuchos::RefCountPtr<NOX::Epetra::Group> grp =
      Teuchos::rcp(new NOX::Epetra::Group(printParams, interface, noxSoln, linSys));

    // Convergence Tests
    Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo = CreateStatusTest(nlParams, grp);

    // Create the solver
    Teuchos::RefCountPtr<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp,combo,RefCountPtr<ParameterList>(&nlParams,false));

#if 0
    if ((step_ % 10) == 0)
    {
      Teuchos::ParameterList& fdParams = nlParams.sublist("Finite Difference");
      double alpha = fdParams.get("alpha", 1.0e-6);
      double beta  = fdParams.get("beta",  1.0e-4);

      ostringstream filename;
      filename << allfiles.outputfile_kenner << "_1_" << step_ << ".dump";
      FSI::Utils::DumpJacobian(*this, alpha, beta, soln, filename.str());
    }
#endif

    // solve the whole thing
    NOX::StatusTest::StatusType status = solver->solve();

    // sometimes we might want to do it again
    if (status != NOX::StatusTest::Converged and secondsolver)
    {
      if (comm_.MyPID()==0)
        utils_->out() << YELLOW "second solver" END_COLOR << endl;

      // Get the Epetra_Vector with the final solution from the solver
      const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
      const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();

      // Start the second solver from the final solution of the first
      // one. Remember that noxSoln is just a view to soln.
      *soln = finalSolution;

      // Create the linear system
      linSys = CreateLinearSystem(nlParams.sublist("Second"), interface, noxSoln, utils_);

      // Create the Group
      grp = Teuchos::rcp(new NOX::Epetra::Group(printParams, interface, noxSoln, linSys));

      // Convergence Tests
      combo = CreateStatusTest(nlParams.sublist("Second"), grp);

      // Create the solver
      solver = NOX::Solver::buildSolver(grp, combo, RefCountPtr<ParameterList>(&nlParams.sublist("Second"),false));

      // solve the whole thing again
      status = solver->solve();
    }

    if (status != NOX::StatusTest::Converged)
      if (comm_.MyPID()==0)
        utils_->out() << RED "Nonlinear solver failed to converge!" END_COLOR << endl;

    // Get the Epetra_Vector with the final solution from the solver
    const NOX::Epetra::Group& finalGroup = dynamic_cast<const NOX::Epetra::Group&>(solver->getSolutionGroup());
    const Epetra_Vector& finalSolution = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getX())).getEpetraVector();
    //const Epetra_Vector& finalF        = (dynamic_cast<const NOX::Epetra::Vector&>(finalGroup.getF())).getEpetraVector();

#if 0
    if ((step_ % 10) == 0)
    {
      Teuchos::ParameterList& fdParams = nlParams.sublist("Finite Difference");
      double alpha = fdParams.get("alpha", 1.0e-6);
      double beta  = fdParams.get("beta",  1.0e-4);

      ostringstream filename;
      filename << allfiles.outputfile_kenner << "_2_" << step_ << ".dump";
      *soln = finalSolution;
      FSI::Utils::DumpJacobian(*this, alpha, beta, soln, filename.str());
    }
#endif

    if (displacementcoupling_)
      idispn_->Update(1.0, finalSolution, 0.0);
    else
      idispn_ = InterfaceDisp();

    // End Nonlinear Solver **************************************

    // Output the parameter list
    if (utils_->isPrintType(NOX::Utils::Parameters))
      if (step_==1 and comm_.MyPID()==0)
      {
        utils_->out() << endl
                      << "Final Parameters" << endl
                      << "****************" << endl;
        solver->getList().print(utils_->out());
        utils_->out() << endl;
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RefCountPtr<NOX::Epetra::LinearSystem>
FSI::DirichletNeumannCoupling::CreateLinearSystem(ParameterList& nlParams,
                                                  const Teuchos::RefCountPtr<NOX::Epetra::Interface::Required>& interface,
                                                  NOX::Epetra::Vector& noxSoln,
                                                  Teuchos::RefCountPtr<NOX::Utils> utils)
{
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method","Aitken"));
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

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

  Teuchos::RefCountPtr<NOX::Epetra::LinearSystem> linSys;

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
    mfresitemax_ = mfParams.get("itemax", -1);

    bool kelleyPerturbation = mfParams.get("Kelley Perturbation", false);

    // MatrixFree seems to be the most interessting choice. But you
    // must set a rather low tolerance for the linear solver.

    MF = Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface, noxSoln, kelleyPerturbation));
    MF->setLambda(lambda);
    iJac = MF;
    J = MF;
  }

  // No Jacobian at all. Do a fix point iteration. This is a user
  // extension, so we have to modify the parameter list here.
  else if (jacobian=="None")
  {
    dirParams.set("Method","User Defined");
    Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
      Teuchos::rcp(new NOX::FSI::FixPointFactory());
    dirParams.set("User Defined Direction Factory",fixpointfactory);
    lsParams.set("Preconditioner","None");
    preconditioner="None";
  }

  // Minimal Polynomial vector extrapolation
  else if (jacobian=="MPE")
  {
    dirParams.set("Method","User Defined");
    Teuchos::RCP<NOX::Direction::UserDefinedFactory> mpefactory =
      Teuchos::rcp(new NOX::FSI::MinimalPolynomialFactory());
    dirParams.set("User Defined Direction Factory",mpefactory);
    lsParams.set("Preconditioner","None");
    preconditioner="None";
  }

#if 0
  // epsilon vector extrapolation
  else if (jacobian=="Epsilon")
  {
    Teuchos::RefCountPtr<NOX::Direction::Generic> mpe = Teuchos::rcp(new NOX::FSI::EpsilonExtrapolation(utils,nlParams));
    dirParams.set("Method","User Defined");
    dirParams.set("User Defined Direction",mpe);
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
#endif

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
      linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, interface, noxSoln));
    }
    else
    {
      // there are different linear solvers available
#if 0
      linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, interface, iJac, J, noxSoln));
#else
      linSys = Teuchos::rcp(new NOX::FSI::LinearSystemGCR(printParams, lsParams, interface, iJac, J, noxSoln));
#endif
    }
  }

  else if (preconditioner=="Dump Finite Difference")
  {
    if (lsParams.get("Preconditioner", "None")=="None")
    {
      if (comm_.MyPID()==0)
        utils->out() << RED "Warning: Preconditioner turned off in linear solver settings." END_COLOR "\n";
    }

    Teuchos::ParameterList& fdParams = nlParams.sublist("Finite Difference");
    //fdresitemax_ = fdParams.get("itemax", -1);
    double alpha = fdParams.get("alpha", 1.0e-4);
    double beta  = fdParams.get("beta",  1.0e-6);

    Teuchos::RefCountPtr<NOX::Epetra::FiniteDifference> precFD =
      Teuchos::rcp(new NOX::Epetra::FiniteDifference(printParams, interface, noxSoln, rawGraph_, beta, alpha));
    iPrec = precFD;
    M = precFD;

    linSys = Teuchos::rcp(new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iJac, J, iPrec, M, noxSoln));
  }

  else
  {
    dserror("unsupported preconditioner '%s'",preconditioner.c_str());
  }

  return linSys;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RefCountPtr<NOX::StatusTest::Combo>
FSI::DirichletNeumannCoupling::CreateStatusTest(ParameterList& nlParams,
                                                Teuchos::RefCountPtr<NOX::Epetra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> combo       = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RefCountPtr<NOX::StatusTest::Combo> converged   = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  Teuchos::RefCountPtr<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RefCountPtr<NOX::StatusTest::FiniteValue> fv    = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  Teuchos::RefCountPtr<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(nlParams.get("Norm abs F", 1.0e-6)));
  converged->addStatusTest(absresid);

  if (nlParams.isParameter("Norm Update"))
  {
    Teuchos::RefCountPtr<NOX::StatusTest::NormUpdate> update =
      Teuchos::rcp(new NOX::StatusTest::NormUpdate(nlParams.get("Norm Update", 1.0e-5)));
    converged->addStatusTest(update);
  }

  if (nlParams.isParameter("Norm rel F"))
  {
    Teuchos::RefCountPtr<NOX::StatusTest::NormF> relresid =
      Teuchos::rcp(new NOX::StatusTest::NormF(*grp.get(), nlParams.get("Norm rel F", 1.0e-2)));
    converged->addStatusTest(relresid);
  }

  //Teuchos::RefCountPtr<NOX::StatusTest::NormWRMS> wrms     = Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  //converged->addStatusTest(wrms);

  return combo;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RefCountPtr<Epetra_Vector> FSI::DirichletNeumannCoupling::InterfaceDisp()
{
  // extract displacements
  return structure_->ExtractInterfaceDisplacement();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RefCountPtr<Epetra_Vector> FSI::DirichletNeumannCoupling::InterfaceForce()
{
  // extract forces
  return FluidToStruct(fluid_->ExtractInterfaceForces());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::DirichletNeumannCoupling::computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  char* flags[] = { "Residual", "Jac", "Prec", "FD_Res", "MF_Res", "MF_Jac", "User", NULL };

  Epetra_Time timer(x.Comm());
  double startTime = timer.WallTime();

  if (comm_.MyPID()==0)
    utils_->out() << "\n "
                  << YELLOW_LIGHT << "FSI residual calculation" << END_COLOR
                  << ".\n fillFlag = " RED << flags[fillFlag] << END_COLOR "\n";

  // we count the number of times the residuum is build
  counter_[fillFlag] += 1;

  if (!x.Map().UniqueGIDs())
    dserror("source map not unique");

  if (displacementcoupling_)
  {
#if 0
    if (fillFlag!=User)
  if (comm_.NumProc()==1)
  {
    static int in_counter;
    std::ostringstream filename;
    filename << allfiles.outputfile_kenner
             << ".x"
             << "." << in_counter
             << ".plot";

    std::cout << "write '" YELLOW_LIGHT << filename.str() << END_COLOR "'\n";
    std::ofstream out(filename.str().c_str());
    for (int i=0; i<x.MyLength()-1; i+=genprob.ndim)
    {
      out << i;
      for (int j=0; j<genprob.ndim; ++j)
        out << " " << x[i+j];
      out << "\n";
    }
    in_counter += 1;
  }
#endif

    Teuchos::RefCountPtr<Epetra_Vector> idispn = rcp(new Epetra_Vector(x));

    Teuchos::RefCountPtr<Epetra_Vector> iforce = FluidOp(idispn, fillFlag);
    Teuchos::RefCountPtr<Epetra_Vector> idispnp = StructOp(iforce, fillFlag);

#if 0
  if (comm_.NumProc()==1)
  {
    static int f_counter;
    std::ostringstream filename;
    filename << allfiles.outputfile_kenner
             << ".i"
             << "." << f_counter
             << ".plot";

    std::cout << "write '" YELLOW_LIGHT << filename.str() << END_COLOR "'\n";
    std::ofstream out(filename.str().c_str());
    for (int i=0; i<iforce->MyLength()-1; i+=2)
    {
      out << i << " " << (*iforce)[i] << " " << (*iforce)[i+1] << "\n";
    }
    f_counter += 1;
  }
#endif

#if 0
    if (fillFlag!=User)
  if (comm_.NumProc()==1)
  {
    static int f_counter;
    std::ostringstream filename;
    filename << allfiles.outputfile_kenner
             << ".d"
             << "." << f_counter
             << ".plot";

    std::cout << "write '" YELLOW_LIGHT << filename.str() << END_COLOR "'\n";
    std::ofstream out(filename.str().c_str());
    for (int i=0; i<idispnp->MyLength()-1; i+=genprob.ndim)
    {
      out << i;
      for (int j=0; j<genprob.ndim; ++j)
        out << " " << (*idispnp)[i+j];
      out << "\n";
    }
    f_counter += 1;
  }
#endif

  F.Update(1.0, *idispnp, -1.0, *idispn, 0.0);

#if 0
  if (comm_.NumProc()==1)
  {
    static int out_counter;
    std::ostringstream filename;
    filename << allfiles.outputfile_kenner
             << ".f"
             << "." << out_counter
             << ".plot";

    std::cout << "write '" YELLOW_LIGHT << filename.str() << END_COLOR "'\n";
    std::ofstream out(filename.str().c_str());
    for (int i=0; i<F.MyLength()-1; i+=2)
    {
      out << i << " " << F[i] << " " << F[i+1] << "\n";
    }
    out_counter += 1;
  }
#endif

  }
  else
  {
    Teuchos::RefCountPtr<Epetra_Vector> iforcen = rcp(new Epetra_Vector(x));

    Teuchos::RefCountPtr<Epetra_Vector> idisp = StructOp(iforcen, fillFlag);
    Teuchos::RefCountPtr<Epetra_Vector> iforcenp = FluidOp(idisp, fillFlag);

    F.Update(1.0, *iforcenp, -1.0, *iforcen, 0.0);
  }

  double endTime = timer.WallTime();
  if (comm_.MyPID()==0)
    utils_->out() << "\nTime for residual calculation: " << endTime-startTime << "\n\n";
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RefCountPtr<Epetra_Vector>
FSI::DirichletNeumannCoupling::FluidOp(Teuchos::RefCountPtr<Epetra_Vector> idisp,
                                       const FillType fillFlag)
{
  if (comm_.MyPID()==0 and utils_->isPrintType(NOX::Utils::OuterIteration))
    utils_->out() << "\nFluid operator\n";

  if (fillFlag==User)
  {
    // SD relaxation calculation

    // Here we have a mesh position independent of the
    // given trial vector, but still the grid velocity depends on the
    // trial vector only.

    // grid velocity
    ale_->ApplyInterfaceDisplacements(StructToAle(idisp));
    ale_->Solve();
    Teuchos::RefCountPtr<Epetra_Vector> fluiddisp = AleToFluid(ale_->ExtractDisplacement());
    fluiddisp->Scale(1./dt_);

    fluid_->ApplyMeshVelocity(fluiddisp);

    // grid position is done inside RelaxationSolve

    // the displacement -> velocity conversion at the interface
    Teuchos::RefCountPtr<Epetra_Vector> ivel = rcp(new Epetra_Vector(*idisp));
    ivel->Scale(1./dt_);

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

    int itemax = fluid_->Itemax();
    if (fillFlag==MF_Res and mfresitemax_ > 0)
      fluid_->SetItemax(mfresitemax_ + 1);

    fluid_->NonlinearSolve();
    fluid_->SetItemax(itemax);

    return FluidToStruct(fluid_->ExtractInterfaceForces());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RefCountPtr<Epetra_Vector>
FSI::DirichletNeumannCoupling::StructOp(Teuchos::RefCountPtr<Epetra_Vector> iforce,
                                        const FillType fillFlag)
{
  if (comm_.MyPID()==0 and utils_->isPrintType(NOX::Utils::OuterIteration))
    utils_->out() << "\nStructural operator\n";

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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RefCountPtr<Epetra_Vector> FSI::DirichletNeumannCoupling::InterfaceVelocity(Teuchos::RefCountPtr<Epetra_Vector> idispnp)
{
  Teuchos::RefCountPtr<Epetra_Vector> ivel = rcp(new Epetra_Vector(*idispn_));
  ivel->Update(1.0, *idispnp, -1.0);
  ivel->Scale(1./dt_);
  return ivel;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RefCountPtr<Epetra_Vector> FSI::DirichletNeumannCoupling::StructToAle(Teuchos::RefCountPtr<Epetra_Vector> iv)
{
  if (matchingnodes_)
  {
    return coupsa_.MasterToSlave(iv);
  }
  else
  {
    // We cannot go from structure to ale directly. So go via the fluid field.
    Teuchos::RefCountPtr<Epetra_Vector> fdisp = coupsfm_.MasterToSlave(iv);
    return coupsa_.MasterToSlave(fdisp);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RefCountPtr<Epetra_Vector> FSI::DirichletNeumannCoupling::StructToFluid(Teuchos::RefCountPtr<Epetra_Vector> iv)
{
  if (matchingnodes_)
  {
    return coupsf_.MasterToSlave(iv);
  }
  else
  {
    return coupsfm_.MasterToSlave(iv);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RefCountPtr<Epetra_Vector> FSI::DirichletNeumannCoupling::FluidToStruct(Teuchos::RefCountPtr<Epetra_Vector> iv)
{
  if (matchingnodes_)
  {
    return coupsf_.SlaveToMaster(iv);
  }
  else
  {
    // Translate consistent nodal forces to interface loads
    RefCountPtr<Epetra_Vector> ishape = fluid_->IntegrateInterfaceShape();
    RefCountPtr<Epetra_Vector> iforce = rcp(new Epetra_Vector(iv->Map()));

    if ( iforce->ReciprocalMultiply( 1.0, *ishape, *iv, 0.0 ) )
      dserror("ReciprocalMultiply failed");

    return coupsfm_.SlaveToMaster(iforce);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RefCountPtr<Epetra_Vector> FSI::DirichletNeumannCoupling::AleToFluid(Teuchos::RefCountPtr<Epetra_Vector> iv)
{
  return coupfa_.SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannCoupling::Update()
{
  structure_->UpdateandOutput();
  fluid_    ->TimeUpdate();
  ale_      ->Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannCoupling::Output()
{
  fluid_->Output();
  ale_  ->Output();
}


#endif
