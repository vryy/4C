
#ifdef CCADISCRET

#include "fsi_dirichletneumann.H"
#include "fsi_utils.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_validparameters.H"

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
#include <Teuchos_StandardParameterEntryValidators.hpp>


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
FSI::DirichletNeumannCoupling::DirichletNeumannCoupling(Epetra_Comm& comm)
  : Algorithm(comm),
    counter_(7)
{
  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  displacementcoupling_ = fsidyn.get<std::string>("COUPVARIABLE") == "Displacement";

  SetDefaultParameters(fsidyn,noxparameterlist);

  //cout << structure_->Discretization();
  //cout << fluid_->Discretization();
  //cout << ale_->Discretization();

  if (Teuchos::getIntegralValue<int>(fsidyn,"COUPMETHOD"))
  {
    matchingnodes_ = true;
    coupsf_.SetupConditionCoupling(StructureField().Discretization(),
                                   StructureField().Interface(),
                                   FluidField().Discretization(),
                                   FluidField().Interface(),
                                   "FSICoupling");

    coupsa_.SetupConditionCoupling(StructureField().Discretization(),
                                   StructureField().Interface(),
                                   *AleField().Discretization(),
                                   AleField().Interface(),
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
    StructureField().SetInterfaceMap(coupsf_.MasterDofMap());
    FluidField().    SetInterfaceMap(coupsf_.SlaveDofMap());
    AleField().      SetInterfaceMap(coupsa_.SlaveDofMap());
  }
  else
  {
    matchingnodes_ = false;
    coupsfm_.Setup( StructureField().Discretization(),
                    FluidField().Discretization(),
                    comm );

    // This is cheating. We setup the coupling of interface dofs between fluid
    // and ale. But we use the variable from the matching version.
    coupsa_.SetupConditionCoupling(FluidField().Discretization(),
                                   FluidField().Interface(),
                                   *AleField().Discretization(),
                                   AleField().Interface(),
                                   "FSICoupling");

    // init transfer from interface to field
    StructureField().SetInterfaceMap(coupsfm_.MasterDofMap());
    FluidField().    SetInterfaceMap(coupsfm_.SlaveDofMap());
    AleField().      SetInterfaceMap(coupsa_.SlaveDofMap());
  }

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = FluidField().Discretization().NodeRowMap();
  const Epetra_Map* alenodemap   = AleField().Discretization()->NodeRowMap();

  coupfa_.SetupCoupling(FluidField().Discretization(),
                        *AleField().Discretization(),
                        *fluidnodemap,
                        *alenodemap);

  FluidField().SetMeshMap(coupfa_.MasterDofMap());

#if 0
  // create connection graph of interface elements
  Teuchos::RCP<Epetra_Map> imap = StructureField().InterfaceMap();

  vector<int> rredundant;
  DRT::UTILS::AllreduceEMap(rredundant, *imap);

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
void FSI::DirichletNeumannCoupling::SetDefaultParameters(const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list)
{
  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = list;

  nlParams.set("Nonlinear Solver", "Line Search Based");
  nlParams.set("Preconditioner", "None");
  nlParams.set("Norm abs F", fsidyn.get<double>("CONVTOL"));
  nlParams.set("Max Iterations", fsidyn.get<int>("ITEMAX"));

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& lineSearchParams = nlParams.sublist("Line Search");

  //
  // Set parameters for NOX to chose the solver direction and line
  // search step.
  //

  switch (Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO"))
  {
  case fsi_iter_stagg_fixed_rel_param:
  {
    // fixed-point solver with fixed relaxation parameter
    SetMethod("ITERATIVE STAGGERED SCHEME WITH FIXED RELAXATION PARAMETER");

    nlParams.set("Jacobian", "None");

    dirParams.set("Method","User Defined");
    Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
      Teuchos::rcp(new NOX::FSI::FixPointFactory());
    dirParams.set("User Defined Direction Factory",fixpointfactory);

    //Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
    //lsParams.set("Preconditioner","None");

    lineSearchParams.set("Method", "Full Step");
    lineSearchParams.sublist("Full Step").set("Full Step", fsidyn.get<double>("RELAX"));
    break;
  }
  case fsi_iter_stagg_AITKEN_rel_param:
  {
    // fixed-point solver with Aitken relaxation parameter
    SetMethod("ITERATIVE STAGGERED SCHEME WITH RELAXATION PARAMETER VIA AITKEN ITERATION");

    nlParams.set("Jacobian", "None");

    dirParams.set("Method","User Defined");
    Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
      Teuchos::rcp(new NOX::FSI::FixPointFactory());
    dirParams.set("User Defined Direction Factory",fixpointfactory);

    Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> aitkenfactory =
      Teuchos::rcp(new NOX::FSI::AitkenFactory());
    lineSearchParams.set("Method","User Defined");
    lineSearchParams.set("User Defined Line Search Factory", aitkenfactory);

    lineSearchParams.sublist("Aitken").set("max step size", fsidyn.get<double>("MAXOMEGA"));
    break;
  }
  case fsi_iter_stagg_steep_desc:
  {
    // fixed-point solver with steepest descent relaxation parameter
    SetMethod("ITERATIVE STAGGERED SCHEME WITH RELAXATION PARAMETER VIA STEEPEST DESCENT METHOD");

    nlParams.set("Jacobian", "None");

    dirParams.set("Method","User Defined");
    Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
      Teuchos::rcp(new NOX::FSI::FixPointFactory());
    dirParams.set("User Defined Direction Factory",fixpointfactory);

    Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> sdfactory =
      Teuchos::rcp(new NOX::FSI::SDFactory());
    lineSearchParams.set("Method","User Defined");
    lineSearchParams.set("User Defined Line Search Factory", sdfactory);
    break;
  }
  case fsi_iter_stagg_NLCG:
  {
    // nonlinear CG solver (pretty much steepest descent with finite
    // difference Jacobian)
    SetMethod("ITERATIVE STAGGERED SCHEME WITH NONLINEAR CG SOLVER");

    nlParams.set("Jacobian", "None");
    dirParams.set("Method", "NonlinearCG");
    lineSearchParams.set("Method", "NonlinearCG");
    break;
  }
  case fsi_iter_stagg_MFNK_FD:
  {
    // matrix free Newton Krylov with finite difference Jacobian
    SetMethod("MATRIX FREE NEWTON KRYLOV SOLVER BASED ON FINITE DIFFERENCES");

    nlParams.set("Jacobian", "Matrix Free");

    Teuchos::ParameterList& mfParams = nlParams.sublist("Matrix Free");
    mfParams.set("lambda", 1.0e-4);
    mfParams.set("itemax", -1);
    mfParams.set("Kelley Perturbation", false);

    lineSearchParams.set("Method", "Full Step");
    lineSearchParams.sublist("Full Step").set("Full Step", 1.0);
    break;
  }
  case fsi_iter_stagg_MFNK_FSI:
  {
    // matrix free Newton Krylov with FSI specific Jacobian
    SetMethod("MATRIX FREE NEWTON KRYLOV SOLVER BASED ON FSI SPECIFIC JACOBIAN APPROXIMATION");

    nlParams.set("Jacobian", "FSI Matrix Free");

    lineSearchParams.set("Method", "Full Step");
    lineSearchParams.sublist("Full Step").set("Full Step", 1.0);
    break;
  }
  case fsi_iter_stagg_MPE:
  {
    // minimal polynomial extrapolation
    SetMethod("ITERATIVE STAGGERED SCHEME WITH MINIMAL POLYNOMIAL EXTRAPOLATION");

    nlParams.set("Jacobian", "None");
    dirParams.set("Method","User Defined");

    Teuchos::RCP<NOX::Direction::UserDefinedFactory> factory =
      Teuchos::rcp(new NOX::FSI::MinimalPolynomialFactory());
    dirParams.set("User Defined Direction Factory",factory);

    Teuchos::ParameterList& exParams = dirParams.sublist("Extrapolation");
    exParams.set("Tolerance", 1e-01);
    exParams.set("omega", fsidyn.get<double>("RELAX"));
    exParams.set("kmax", 10);
    exParams.set("Method", "MPE");

    //lsParams.set("Preconditioner","None");

    lineSearchParams.set("Method", "Full Step");
    lineSearchParams.sublist("Full Step").set("Full Step", 1.0);
    break;
  }
  case fsi_iter_stagg_RRE:
  {
    // reduced rank extrapolation
    SetMethod("ITERATIVE STAGGERED SCHEME WITH REDUCED RANK EXTRAPOLATION");

    nlParams.set("Jacobian", "None");
    dirParams.set("Method","User Defined");

    Teuchos::RCP<NOX::Direction::UserDefinedFactory> factory =
      Teuchos::rcp(new NOX::FSI::MinimalPolynomialFactory());
    dirParams.set("User Defined Direction Factory",factory);

    Teuchos::ParameterList& exParams = dirParams.sublist("Extrapolation");
    exParams.set("Tolerance", 1e-01);
    exParams.set("omega", fsidyn.get<double>("RELAX"));
    exParams.set("kmax", 10);
    exParams.set("Method", "RRE");

    //lsParams.set("Preconditioner","None");

    lineSearchParams.set("Method", "Full Step");
    lineSearchParams.sublist("Full Step").set("Full Step", 1.0);
    break;
  }
  case fsi_iter_nox:
    dserror("obsolete");
    break;
  case fsi_iter_monolithic:
    dserror("No monolithic coupling with Dirichlet-Neumann partitioning. Panic.");
    break;
  case fsi_basic_sequ_stagg:
  {
    // sequential coupling (no iteration!)
    SetMethod("BASIC SEQUENTIAL STAGGERED SCHEME");

    nlParams.set("Jacobian", "None");
    nlParams.set("Max Iterations", 1.);

    dirParams.set("Method","User Defined");
    Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
      Teuchos::rcp(new NOX::FSI::FixPointFactory());
    dirParams.set("User Defined Direction Factory",fixpointfactory);

    lineSearchParams.set("Method", "Full Step");
    lineSearchParams.sublist("Full Step").set("Full Step", 1.0);
    break;
  }
  case fsi_sequ_stagg_pred:
  case fsi_sequ_stagg_shift:
  default:
    dserror("coupling method type '%s' unsupported", fsidyn.get<string>("COUPALGO").c_str());
  }

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm().MyPID());

  // set default output flag to no output
  // The field solver will output a lot, anyway.
  printParams.get("Output Information",
                  ::NOX::Utils::Warning |
                  ::NOX::Utils::OuterIteration |
                  ::NOX::Utils::OuterIterationStatusTest
    );

}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::DirichletNeumannCoupling::Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
  // not such a smart idea?!
  bool secondsolver = false;

  const Teuchos::ParameterList& fsidyn   = DRT::Problem::Instance()->FSIDynamicParams();

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = noxparameterlist;

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method","Newton"));
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");

  // Create printing utilities
  utils_ = Teuchos::rcp(new NOX::Utils(printParams));

  // ==================================================================

  // log solver iterations

  Teuchos::RCP<std::ofstream> log;
  if (Comm().MyPID()==0)
  {
    std::string s = allfiles.outputfile_kenner;
    s.append(".iteration");
    log = Teuchos::rcp(new std::ofstream(s.c_str()));
    (*log) << "# num procs      = " << Comm().NumProc() << "\n"
           << "# Method         = " << nlParams.sublist("Direction").get("Method","Newton") << "\n"
           << "# Jacobian       = " << nlParams.get("Jacobian", "None") << "\n"
           << "# Preconditioner = " << nlParams.get("Preconditioner","None") << "\n"
           << "# Line Search    = " << nlParams.sublist("Line Search").get("Method","Aitken") << "\n"
           << "# Predictor      = '" << fsidyn.get<std::string>("PREDICTOR") << "'\n"
           << "#\n"
           << "# step  time/step  #nliter  |R|  #liter  Residual  Jac  Prec  FD_Res  MF_Res  MF_Jac  User\n"
      ;
  }

  // get an idea of interface displacement
  idispn_ = StructureField().ExtractInterfaceDispn();

  Teuchos::Time timer("time step timer");

  // ==================================================================

  while (NotFinished())
  {
    // Increment all field counters and predict field values whenever
    // appropriate.
    PrepareTimeStep();

    // reset all counters
    std::fill(counter_.begin(),counter_.end(),0);
    lsParams.sublist("Output").set("Total Number of Linear Iterations",0);
    linsolvcount_.resize(0);

    // start time measurement
    Teuchos::RCP<Teuchos::TimeMonitor> timemonitor = rcp(new Teuchos::TimeMonitor(timer,true));

    /*----------------- CSD - predictor for itnum==0 --------------------*/

    // Begin Nonlinear Solver ************************************

    // Get initial guess.
    Teuchos::RCP<Epetra_Vector> soln;
    if (displacementcoupling_)
    {
      // predict displacement
      soln = StructureField().PredictInterfaceDisplacement();
    }
    else
    {
      if (Teuchos::getIntegralValue<int>(fsidyn,"PREDICTOR")!=1)
      {
        dserror("unknown interface force predictor '%s'",
                fsidyn.get<string>("PREDICTOR").c_str());
      }
      soln = InterfaceForce();
    }

    NOX::Epetra::Vector noxSoln(soln, NOX::Epetra::Vector::CreateView);

    // Create the linear system
    Teuchos::RCP<NOX::Epetra::LinearSystem> linSys =
      CreateLinearSystem(nlParams, interface, noxSoln, utils_);

    // Create the Group
    Teuchos::RCP<NOX::Epetra::Group> grp =
      Teuchos::rcp(new NOX::Epetra::Group(printParams, interface, noxSoln, linSys));

    // Convergence Tests
    Teuchos::RCP<NOX::StatusTest::Combo> combo = CreateStatusTest(nlParams, grp);

    // Create the solver
    Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(grp,combo,RCP<ParameterList>(&nlParams,false));

#if 0
    if ((step_ % 10) == 0)
    {
      Teuchos::ParameterList& fdParams = nlParams.sublist("Finite Difference");
      double alpha = fdParams.get("alpha", 1.0e-6);
      double beta  = fdParams.get("beta",  1.0e-4);

      ostringstream filename;
      filename << allfiles.outputfile_kenner << "_1_" << step_ << ".dump";
      FSI::UTILS::DumpJacobian(*this, alpha, beta, soln, filename.str());
    }
#endif

    // solve the whole thing
    NOX::StatusTest::StatusType status = solver->solve();

    // sometimes we might want to do it again
    if (status != NOX::StatusTest::Converged and secondsolver)
    {
      if (Comm().MyPID()==0)
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
      solver = NOX::Solver::buildSolver(grp, combo, RCP<ParameterList>(&nlParams.sublist("Second"),false));

      // solve the whole thing again
      status = solver->solve();
    }

    if (status != NOX::StatusTest::Converged)
      if (Comm().MyPID()==0)
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
      FSI::UTILS::DumpJacobian(*this, alpha, beta, soln, filename.str());
    }
#endif

    if (displacementcoupling_)
      idispn_->Update(1.0, finalSolution, 0.0);
    else
      idispn_ = InterfaceDisp();

    // End Nonlinear Solver **************************************

    // Output the parameter list
    if (utils_->isPrintType(NOX::Utils::Parameters))
      if (Step()==1 and Comm().MyPID()==0)
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

    if (Comm().MyPID()==0)
    {
      (*log) << Step()
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
Teuchos::RCP<NOX::Epetra::LinearSystem>
FSI::DirichletNeumannCoupling::CreateLinearSystem(ParameterList& nlParams,
                                                  const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface,
                                                  NOX::Epetra::Vector& noxSoln,
                                                  Teuchos::RCP<NOX::Utils> utils)
{
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method","Aitken"));
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  Teuchos::RCP<NOX::FSI::FSIMatrixFree> FSIMF;
  Teuchos::RCP<NOX::Epetra::MatrixFree> MF;
  Teuchos::RCP<NOX::Epetra::FiniteDifference> FD;
  Teuchos::RCP<NOX::Epetra::FiniteDifferenceColoring> FDC;
  Teuchos::RCP<NOX::Epetra::FiniteDifferenceColoring> FDC1;
  Teuchos::RCP<NOX::Epetra::BroydenOperator> B;

  Teuchos::RCP<NOX::Epetra::Interface::Jacobian> iJac;
  Teuchos::RCP<NOX::Epetra::Interface::Preconditioner> iPrec;

  Teuchos::RCP<Epetra_Operator> J;
  Teuchos::RCP<Epetra_Operator> M;

  Teuchos::RCP<NOX::Epetra::LinearSystem> linSys;

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

  // Matrix Free Newton Krylov.
  else if (jacobian=="Matrix Free")
  {
    Teuchos::ParameterList& mfParams = nlParams.sublist("Matrix Free");
    double lambda = mfParams.get("lambda", 1.0e-4);
    mfresitemax_ = mfParams.get("itemax", -1);

    bool kelleyPerturbation = mfParams.get("Kelley Perturbation", false);

    // MatrixFree seems to be the most interessting choice. But you
    // must set a rather low tolerance for the linear solver.

    MF = Teuchos::rcp(new NOX::Epetra::MatrixFree(printParams, interface, noxSoln, kelleyPerturbation));
    MF->setLambda(lambda);
    iJac = MF;
    J = MF;
  }

  // No Jacobian at all. Do a fix point iteration.
  else if (jacobian=="None")
  {
    preconditioner="None";
  }

  // This is pretty much debug code. Or rather research code.
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

  // No preconditioning at all.
  if (preconditioner=="None")
  {
    if (Teuchos::is_null(iJac))
    {
      // if no Jacobian has been set this better be the fix point
      // method.
      if (dirParams.get("Method","Newton")!="User Defined")
      {
        if (Comm().MyPID()==0)
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
      if (Comm().MyPID()==0)
        utils->out() << RED "Warning: Preconditioner turned off in linear solver settings." END_COLOR "\n";
    }

    Teuchos::ParameterList& fdParams = nlParams.sublist("Finite Difference");
    //fdresitemax_ = fdParams.get("itemax", -1);
    double alpha = fdParams.get("alpha", 1.0e-4);
    double beta  = fdParams.get("beta",  1.0e-6);

    Teuchos::RCP<NOX::Epetra::FiniteDifference> precFD =
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
Teuchos::RCP<NOX::StatusTest::Combo>
FSI::DirichletNeumannCoupling::CreateStatusTest(ParameterList& nlParams,
                                                Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo       = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::StatusTest::Combo> converged   = Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters = Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv    = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  Teuchos::RCP<NOX::StatusTest::NormF> absresid =
    Teuchos::rcp(new NOX::StatusTest::NormF(nlParams.get("Norm abs F", 1.0e-6)));
  converged->addStatusTest(absresid);

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

  //Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms     = Teuchos::rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
  //converged->addStatusTest(wrms);

  return combo;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannCoupling::InterfaceDisp()
{
  // extract displacements
  return StructureField().ExtractInterfaceDisplacement();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannCoupling::InterfaceForce()
{
  // extract forces
  return FluidToStruct(FluidField().ExtractInterfaceForces());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::DirichletNeumannCoupling::computeF(const Epetra_Vector &x, Epetra_Vector &F, const FillType fillFlag)
{
  char* flags[] = { "Residual", "Jac", "Prec", "FD_Res", "MF_Res", "MF_Jac", "User", NULL };

  Epetra_Time timer(x.Comm());
  double startTime = timer.WallTime();

  if (Comm().MyPID()==0)
  {
    utils_->out() << "\n "
                  << YELLOW_LIGHT << "FSI residual calculation" << END_COLOR
                  << ".\n";
    if (fillFlag!=Residual)
      utils_->out() << " fillFlag = " RED << flags[fillFlag] << END_COLOR "\n";
  }

  // we count the number of times the residuum is build
  counter_[fillFlag] += 1;

  if (!x.Map().UniqueGIDs())
    dserror("source map not unique");

  if (displacementcoupling_)
  {
#if 0
    if (fillFlag!=User)
  if (Comm().NumProc()==1)
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

    Teuchos::RCP<Epetra_Vector> idispn = rcp(new Epetra_Vector(x));

    Teuchos::RCP<Epetra_Vector> iforce = FluidOp(idispn, fillFlag);
    Teuchos::RCP<Epetra_Vector> idispnp = StructOp(iforce, fillFlag);

#if 0
  if (Comm().NumProc()==1)
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
  if (Comm().NumProc()==1)
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
  if (Comm().NumProc()==1)
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
    Teuchos::RCP<Epetra_Vector> iforcen = rcp(new Epetra_Vector(x));

    Teuchos::RCP<Epetra_Vector> idisp = StructOp(iforcen, fillFlag);
    Teuchos::RCP<Epetra_Vector> iforcenp = FluidOp(idisp, fillFlag);

    F.Update(1.0, *iforcenp, -1.0, *iforcen, 0.0);
  }

  double endTime = timer.WallTime();
  if (Comm().MyPID()==0)
    utils_->out() << "\nTime for residual calculation: " << endTime-startTime << "\n\n";
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
FSI::DirichletNeumannCoupling::FluidOp(Teuchos::RCP<Epetra_Vector> idisp,
                                       const FillType fillFlag)
{
  if (Comm().MyPID()==0 and utils_->isPrintType(NOX::Utils::OuterIteration))
    utils_->out() << "\nFluid operator\n";

  if (fillFlag==User)
  {
    // SD relaxation calculation

    // Here we have a mesh position independent of the
    // given trial vector, but still the grid velocity depends on the
    // trial vector only.

    // grid velocity
    AleField().ApplyInterfaceDisplacements(StructToAle(idisp));
    AleField().Solve();
    Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluid(AleField().ExtractDisplacement());
    fluiddisp->Scale(1./Dt());

    FluidField().ApplyMeshVelocity(fluiddisp);

    // grid position is done inside RelaxationSolve

    // the displacement -> velocity conversion at the interface
    Teuchos::RCP<Epetra_Vector> ivel = rcp(new Epetra_Vector(*idisp));
    ivel->Scale(1./Dt());

    return FluidToStruct(FluidField().RelaxationSolve(StructToFluid(ivel)));
  }
  else
  {
    // normal fluid solve

    AleField().ApplyInterfaceDisplacements(StructToAle(idisp));
    AleField().Solve();

    // the displacement -> velocity conversion at the interface
    Teuchos::RCP<Epetra_Vector> ivel = InterfaceVelocity(idisp);
    Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluid(AleField().ExtractDisplacement());

    FluidField().ApplyInterfaceVelocities(StructToFluid(ivel));
    FluidField().ApplyMeshDisplacement(fluiddisp);

    int itemax = FluidField().Itemax();
    if (fillFlag==MF_Res and mfresitemax_ > 0)
      FluidField().SetItemax(mfresitemax_ + 1);

    FluidField().NonlinearSolve();
    FluidField().SetItemax(itemax);

    return FluidToStruct(FluidField().ExtractInterfaceForces());
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector>
FSI::DirichletNeumannCoupling::StructOp(Teuchos::RCP<Epetra_Vector> iforce,
                                        const FillType fillFlag)
{
  if (Comm().MyPID()==0 and utils_->isPrintType(NOX::Utils::OuterIteration))
    utils_->out() << "\nStructural operator\n";

  if (fillFlag==User)
  {
    // SD relaxation calculation
    return StructureField().RelaxationSolve(iforce);
  }
  else
  {
    // normal structure solve
    StructureField().ApplyInterfaceForces(iforce);
    StructureField().Solve();
    return StructureField().ExtractInterfaceDisplacement();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannCoupling::InterfaceVelocity(Teuchos::RCP<Epetra_Vector> idispnp)
{
  Teuchos::RCP<Epetra_Vector> ivel = rcp(new Epetra_Vector(*idispn_));
  ivel->Update(1.0, *idispnp, -1.0);
  ivel->Scale(1./Dt());
  return ivel;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannCoupling::StructToAle(Teuchos::RCP<Epetra_Vector> iv)
{
  if (matchingnodes_)
  {
    return coupsa_.MasterToSlave(iv);
  }
  else
  {
    // We cannot go from structure to ale directly. So go via the fluid field.
    Teuchos::RCP<Epetra_Vector> fdisp = coupsfm_.MasterToSlave(iv);
    return coupsa_.MasterToSlave(fdisp);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannCoupling::StructToFluid(Teuchos::RCP<Epetra_Vector> iv)
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
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannCoupling::FluidToStruct(Teuchos::RCP<Epetra_Vector> iv)
{
  if (matchingnodes_)
  {
    return coupsf_.SlaveToMaster(iv);
  }
  else
  {
    // Translate consistent nodal forces to interface loads
    Teuchos::RCP<Epetra_Vector> ishape = FluidField().IntegrateInterfaceShape();
    Teuchos::RCP<Epetra_Vector> iforce = rcp(new Epetra_Vector(iv->Map()));

    if ( iforce->ReciprocalMultiply( 1.0, *ishape, *iv, 0.0 ) )
      dserror("ReciprocalMultiply failed");

    return coupsfm_.SlaveToMaster(iforce);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::DirichletNeumannCoupling::AleToFluid(Teuchos::RCP<Epetra_Vector> iv)
{
  return coupfa_.SlaveToMaster(iv);
}


#endif
