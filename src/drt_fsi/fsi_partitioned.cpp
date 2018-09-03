/*----------------------------------------------------------------------*/
/*!
 * \file fsi_partitioned.cpp
\brief Partitioned FSI base

\level 1

\maintainer Matthias Mayr
            mayr@mhpc.mw.tum.de
            089 - 289-10362

*/
/*----------------------------------------------------------------------*/



#include "fsi_partitioned.H"
#include "fsi_utils.H"

#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_adapter/ad_fld_fluid_xfem.H"
#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_adapter/ad_fld_fluid_xfsi.H"
#include "../drt_adapter/ad_ale_fluid.H"
#include "../drt_adapter/ad_fld_fluid_immersed.H"

#include "../drt_fluid_xfluid/xfluid.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_fsi.H"
#include "../drt_lib/drt_colors.H"

#include "../drt_io/io_control.H"
#include "../drt_structure/stru_aux.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "fsi_debugwriter.H"

#include "fsi_nox_aitken.H"
#include "fsi_nox_fixpoint.H"
#include "fsi_nox_jacobian.H"
#include "fsi_nox_sd.H"
#include "fsi_nox_linearsystem_gcr.H"
#include "fsi_nox_mpe.H"

#include <string>
#include <Epetra_Time.h>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::Partitioned::Partitioned(const Epetra_Comm& comm)
    : Algorithm(comm),
      idispn_(Teuchos::null),
      iveln_(Teuchos::null),
      rawGraph_(Teuchos::null),
      counter_(7),
      mfresitemax_(0),
      coupsfm_(Teuchos::null),
      matchingnodes_(false),
      debugwriter_(Teuchos::null)
{
  // empty constructor
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Partitioned::Setup()
{
  // call setup of base class
  FSI::Algorithm::Setup();

  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  SetDefaultParameters(fsidyn, noxparameterlist_);
  SetupCoupling(fsidyn, Comm());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Partitioned::SetupCoupling(const Teuchos::ParameterList& fsidyn, const Epetra_Comm& comm)
{
  if (Comm().MyPID() == 0) std::cout << "\n SetupCoupling in FSI::Partitioned ..." << std::endl;

  ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  coupsfm_ = Teuchos::rcp(new ADAPTER::CouplingMortar());


  if ((DRT::INPUT::IntegralValue<int>(fsidyn.sublist("PARTITIONED SOLVER"), "COUPMETHOD") ==
          1)  // matching meshes
      and (DRT::Problem::Instance()->ProblemType() != prb_fsi_xfem))
  {
    matchingnodes_ = true;
    const int ndim = DRT::Problem::Instance()->NDim();
    coupsf.SetupConditionCoupling(*StructureField()->Discretization(),
        StructureField()->Interface()->FSICondMap(), *MBFluidField()->Discretization(),
        MBFluidField()->Interface()->FSICondMap(), "FSICoupling", ndim);

    if (coupsf.MasterDofMap()->NumGlobalElements() == 0)
      dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");
  }
  else if ((DRT::INPUT::IntegralValue<int>(fsidyn.sublist("PARTITIONED SOLVER"), "COUPMETHOD") ==
               1)  // matching meshes coupled via XFEM
           and (DRT::Problem::Instance()->ProblemType() == prb_fsi_xfem))
  {
    matchingnodes_ = true;  // matching between structure and boundary dis! non-matching between
                            // boundary dis and fluid is handled bei XFluid itself
    const int ndim = DRT::Problem::Instance()->NDim();

    Teuchos::RCP<ADAPTER::FluidXFEM> x_movingboundary =
        Teuchos::rcp_dynamic_cast<ADAPTER::FluidXFEM>(MBFluidField());
    coupsf.SetupConditionCoupling(*StructureField()->Discretization(),
        StructureField()->Interface()->FSICondMap(),
        *x_movingboundary->BoundaryDiscretization(),  // use the matching boundary discretization
        x_movingboundary->StructInterface()->FSICondMap(), "FSICoupling", ndim);

    if (coupsf.MasterDofMap()->NumGlobalElements() == 0)
      dserror("No nodes in matching FSI interface. Empty FSI coupling condition?");
  }
  else if (DRT::INPUT::IntegralValue<int>(fsidyn.sublist("PARTITIONED SOLVER"), "COUPMETHOD") ==
               0  // mortar coupling
           and (DRT::Problem::Instance()->ProblemType() != prb_fsi_xfem))
  {
    // coupling condition at the fsi interface: displacements (=number of spatial dimensions) are
    // coupled e.g.: 3D: coupleddof = [1, 1, 1]
    std::vector<int> coupleddof(DRT::Problem::Instance()->NDim(), 1);

    matchingnodes_ = false;
    coupsfm_->Setup(StructureField()->Discretization(), MBFluidField()->Discretization(),
        (Teuchos::rcp_dynamic_cast<ADAPTER::FluidAle>(MBFluidField()))
            ->AleField()
            ->WriteAccessDiscretization(),
        coupleddof, "FSICoupling", comm, true);
  }
  else
  {
    dserror("You should not arrive here");
  }

  // enable debugging
  if (DRT::INPUT::IntegralValue<int>(fsidyn, "DEBUGOUTPUT"))
    debugwriter_ = Teuchos::rcp(new UTILS::DebugWriter(StructureField()->Discretization()));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Partitioned::SetDefaultParameters(
    const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list)
{
  if (Comm().MyPID() == 0)
    std::cout << "\n SetDefaultParameters in FSI::Partitioned ..." << std::endl;

  // extract sublist with settings for partitioned solver
  const Teuchos::ParameterList& fsipart = fsidyn.sublist("PARTITIONED SOLVER");

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = list;

  nlParams.set("Nonlinear Solver", "Line Search Based");
  nlParams.set("Preconditioner", "None");
  nlParams.set("Norm abs F", fsipart.get<double>("CONVTOL"));
  nlParams.set("Max Iterations", fsipart.get<int>("ITEMAX"));

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& lineSearchParams = nlParams.sublist("Line Search");

  //
  // Set parameters for NOX to chose the solver direction and line
  // search step.
  //

  switch (DRT::INPUT::IntegralValue<int>(fsidyn, "COUPALGO"))
  {
    case fsi_iter_stagg_fixed_rel_param:
    {
      // fixed-point solver with fixed relaxation parameter
      SetMethod("ITERATIVE STAGGERED SCHEME WITH FIXED RELAXATION PARAMETER");

      nlParams.set("Jacobian", "None");

      dirParams.set("Method", "User Defined");
      Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
          Teuchos::rcp(new NOX::FSI::FixPointFactory());
      dirParams.set("User Defined Direction Factory", fixpointfactory);

      // Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
      // lsParams.set("Preconditioner","None");

      lineSearchParams.set("Method", "Full Step");
      lineSearchParams.sublist("Full Step").set("Full Step", fsipart.get<double>("RELAX"));
      break;
    }
    case fsi_iter_stagg_AITKEN_rel_param:
    {
      // fixed-point solver with Aitken relaxation parameter
      SetMethod("ITERATIVE STAGGERED SCHEME WITH RELAXATION PARAMETER VIA AITKEN ITERATION");

      nlParams.set("Jacobian", "None");

      dirParams.set("Method", "User Defined");
      Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
          Teuchos::rcp(new NOX::FSI::FixPointFactory());
      dirParams.set("User Defined Direction Factory", fixpointfactory);

      Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> linesearchfactory =
          Teuchos::rcp(new NOX::FSI::AitkenFactory());
      lineSearchParams.set("Method", "User Defined");
      lineSearchParams.set("User Defined Line Search Factory", linesearchfactory);

      lineSearchParams.sublist("Aitken").set("max step size", fsipart.get<double>("MAXOMEGA"));
      lineSearchParams.sublist("Aitken").set("min step size", fsipart.get<double>("MINOMEGA"));
      break;
    }
    case fsi_iter_stagg_steep_desc:
    {
      // fixed-point solver with steepest descent relaxation parameter
      SetMethod("ITERATIVE STAGGERED SCHEME WITH RELAXATION PARAMETER VIA STEEPEST DESCENT METHOD");

      nlParams.set("Jacobian", "None");

      dirParams.set("Method", "User Defined");
      Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
          Teuchos::rcp(new NOX::FSI::FixPointFactory());
      dirParams.set("User Defined Direction Factory", fixpointfactory);

      Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> linesearchfactory =
          Teuchos::rcp(new NOX::FSI::SDFactory());
      lineSearchParams.set("Method", "User Defined");
      lineSearchParams.set("User Defined Line Search Factory", linesearchfactory);
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
      mfParams.set("itemax", 1);
      mfParams.set("Kelley Perturbation", false);

      lineSearchParams.set("Method", "Full Step");
      lineSearchParams.sublist("Full Step").set("Full Step", 1.0);

      Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
      Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method", "Newton"));
      Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

      lsParams.set("Tolerance", fsipart.get<double>("BASETOL"));

      break;
    }
    case fsi_iter_stagg_MFNK_FSI:
    {
      // matrix free Newton Krylov with FSI specific Jacobian
      SetMethod("MATRIX FREE NEWTON KRYLOV SOLVER BASED ON FSI SPECIFIC JACOBIAN APPROXIMATION");

      nlParams.set("Jacobian", "FSI Matrix Free");

      lineSearchParams.set("Method", "Full Step");
      lineSearchParams.sublist("Full Step").set("Full Step", 1.0);

      Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
      Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method", "Newton"));
      Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

      lsParams.set("Tolerance", fsipart.get<double>("BASETOL"));

      break;
    }
    case fsi_iter_stagg_MPE:
    {
      // minimal polynomial extrapolation
      SetMethod("ITERATIVE STAGGERED SCHEME WITH MINIMAL POLYNOMIAL EXTRAPOLATION");

      nlParams.set("Jacobian", "None");
      dirParams.set("Method", "User Defined");

      Teuchos::RCP<NOX::Direction::UserDefinedFactory> factory =
          Teuchos::rcp(new NOX::FSI::MinimalPolynomialFactory());
      dirParams.set("User Defined Direction Factory", factory);

      Teuchos::ParameterList& exParams = dirParams.sublist("Extrapolation");
      exParams.set("Tolerance", fsipart.get<double>("BASETOL"));
      exParams.set("omega", fsipart.get<double>("RELAX"));
      exParams.set("kmax", 25);
      exParams.set("Method", "MPE");

      // lsParams.set("Preconditioner","None");

      lineSearchParams.set("Method", "Full Step");
      lineSearchParams.sublist("Full Step").set("Full Step", 1.0);
      break;
    }
    case fsi_iter_stagg_RRE:
    {
      // reduced rank extrapolation
      SetMethod("ITERATIVE STAGGERED SCHEME WITH REDUCED RANK EXTRAPOLATION");

      nlParams.set("Jacobian", "None");
      dirParams.set("Method", "User Defined");

      Teuchos::RCP<NOX::Direction::UserDefinedFactory> factory =
          Teuchos::rcp(new NOX::FSI::MinimalPolynomialFactory());
      dirParams.set("User Defined Direction Factory", factory);

      Teuchos::ParameterList& exParams = dirParams.sublist("Extrapolation");
      exParams.set("Tolerance", fsipart.get<double>("BASETOL"));
      exParams.set("omega", fsipart.get<double>("RELAX"));
      exParams.set("kmax", 25);
      exParams.set("Method", "RRE");

      // lsParams.set("Preconditioner","None");

      lineSearchParams.set("Method", "Full Step");
      lineSearchParams.sublist("Full Step").set("Full Step", 1.0);
      break;
    }
    case fsi_basic_sequ_stagg:
    {
      // sequential coupling (no iteration!)
      SetMethod("BASIC SEQUENTIAL STAGGERED SCHEME");

      nlParams.set("Jacobian", "None");
      nlParams.set("Max Iterations", 1);

      dirParams.set("Method", "User Defined");
      Teuchos::RCP<NOX::Direction::UserDefinedFactory> fixpointfactory =
          Teuchos::rcp(new NOX::FSI::FixPointFactory());
      dirParams.set("User Defined Direction Factory", fixpointfactory);

      lineSearchParams.set("Method", "Full Step");
      lineSearchParams.sublist("Full Step").set("Full Step", 1.0);
      break;
    }
    default:
    {
      dserror("coupling method type '%s' unsupported", fsidyn.get<std::string>("COUPALGO").c_str());
      break;
    }
  }

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm().MyPID());

  // set default output flag to no output
  // The field solver will output a lot, anyway.
  printParams.get("Output Information",
      ::NOX::Utils::Warning | ::NOX::Utils::OuterIteration | ::NOX::Utils::OuterIterationStatusTest
      // ::NOX::Utils::Parameters
  );

  Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
  solverOptions.set<std::string>("Status Test Check Type", "Complete");
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Partitioned::Timeloop(const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface)
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = noxparameterlist_;

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method", "Newton"));
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");

  // Create printing utilities
  utils_ = Teuchos::rcp(new NOX::Utils(printParams));

  // ==================================================================

  // log solver iterations

  Teuchos::RCP<std::ofstream> log;
  if (Comm().MyPID() == 0)
  {
    std::string s = DRT::Problem::Instance()->OutputControlFile()->FileName();
    s.append(".iteration");
    log = Teuchos::rcp(new std::ofstream(s.c_str()));
    (*log) << "# num procs      = " << Comm().NumProc() << "\n"
           << "# Method         = " << nlParams.sublist("Direction").get("Method", "Newton") << "\n"
           << "# Jacobian       = " << nlParams.get("Jacobian", "None") << "\n"
           << "# Preconditioner = " << nlParams.get("Preconditioner", "None") << "\n"
           << "# Line Search    = " << nlParams.sublist("Line Search").get("Method", "Aitken")
           << "\n"
           << "# Predictor      = '"
           << fsidyn.sublist("PARTITIONED SOLVER").get<std::string>("PREDICTOR") << "'\n"
           << "#\n"
           << "# step | time | time/step | #nliter  |R|  #liter  Residual  Jac  Prec  FD_Res  "
              "MF_Res  MF_Jac  User\n";
  }

  // get an idea of interface displacement
  ExtractPreviousInterfaceSolution();

  Teuchos::Time timer("time step timer");

  // ==================================================================

  while (NotFinished())
  {
    // Increment all field counters and predict field values whenever
    // appropriate.
    PrepareTimeStep();

    if (debugwriter_ != Teuchos::null) debugwriter_->NewTimeStep(Step());

    // reset all counters
    std::fill(counter_.begin(), counter_.end(), 0);
    lsParams.sublist("Output").set("Total Number of Linear Iterations", 0);
    linsolvcount_.resize(0);

    // start time measurement
    Teuchos::RCP<Teuchos::TimeMonitor> timemonitor =
        Teuchos::rcp(new Teuchos::TimeMonitor(timer, true));

    /*----------------- CSD - predictor for itnum==0 --------------------*/

    // Begin Nonlinear Solver ************************************

    // Get initial guess
    Teuchos::RCP<Epetra_Vector> soln = InitialGuess();

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
    Teuchos::RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(
        grp, combo, Teuchos::RCP<Teuchos::ParameterList>(&nlParams, false));

    // solve the whole thing
    NOX::StatusTest::StatusType status = solver->solve();

    if (status != NOX::StatusTest::Converged) dserror("Nonlinear solver failed to converge!");

    // End Nonlinear Solver **************************************

    // Output the parameter list
    if (utils_->isPrintType(NOX::Utils::Parameters))
      if (Step() == 1 and Comm().MyPID() == 0)
      {
        utils_->out() << std::endl
                      << "Final Parameters" << std::endl
                      << "****************" << std::endl;
        solver->getList().print(utils_->out());
        utils_->out() << std::endl;
      }

    // ==================================================================

    // stop time measurement
    timemonitor = Teuchos::null;

    if (Comm().MyPID() == 0)
    {
      (*log) << Step() << "\t" << Time() << "\t" << timer.totalElapsedTime() << "\t"
             << nlParams.sublist("Output").get("Nonlinear Iterations", 0) << "\t"
             << nlParams.sublist("Output").get("2-Norm of Residual", 0.) << "\t"
             << lsParams.sublist("Output").get("Total Number of Linear Iterations", 0);
      for (std::vector<int>::size_type i = 0; i < counter_.size(); ++i)
      {
        (*log) << " " << counter_[i];
      }
      (*log) << std::endl;
      log->flush();
    }

    // ==================================================================

    // In case of sliding ALE interfaces, 'remesh' fluid field
    INPAR::FSI::PartitionedCouplingMethod usedmethod =
        DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(
            fsidyn.sublist("PARTITIONED SOLVER"), "PARTITIONED");

    if (usedmethod == INPAR::FSI::DirichletNeumannSlideale)
    {
      Remeshing();
    }

    // calculate stresses, strains, energies
    PrepareOutput();

    // prepare field variables for new time step
    Update();

    // extract final displacement and velocity
    // since we did update, this is very easy to extract
    ExtractPreviousInterfaceSolution();

    // write current solution
    Output();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::Epetra::LinearSystem> FSI::Partitioned::CreateLinearSystem(
    Teuchos::ParameterList& nlParams,
    const Teuchos::RCP<NOX::Epetra::Interface::Required>& interface, NOX::Epetra::Vector& noxSoln,
    Teuchos::RCP<NOX::Utils> utils)
{
  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist(dirParams.get("Method", "Aitken"));
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

  const std::string jacobian = nlParams.get("Jacobian", "None");
  std::string preconditioner = nlParams.get("Preconditioner", "None");

  // Special FSI based matrix free method
  if (jacobian == "FSI Matrix Free")
  {
    // Teuchos::ParameterList& mfParams = nlParams.sublist("FSI Matrix Free");

    // MatrixFree seems to be the most interessting choice. This
    // version builds on our steepest descent relaxation
    // implementation to approximate the Jacobian times x.

    // This is the default method.

    FSIMF = Teuchos::rcp(new NOX::FSI::FSIMatrixFree(printParams, interface, noxSoln));
    iJac = FSIMF;
    J = FSIMF;
  }

  // Matrix Free Newton Krylov.
  else if (jacobian == "Matrix Free")
  {
    Teuchos::ParameterList& mfParams = nlParams.sublist("Matrix Free");
    double lambda = mfParams.get("lambda", 1.0e-4);
    mfresitemax_ = mfParams.get("itemax", -1);

    bool kelleyPerturbation = mfParams.get("Kelley Perturbation", false);

    // MatrixFree seems to be the most interessting choice. But you
    // must set a rather low tolerance for the linear solver.

    MF = Teuchos::rcp(
        new NOX::Epetra::MatrixFree(printParams, interface, noxSoln, kelleyPerturbation));
    MF->setLambda(lambda);
    iJac = MF;
    J = MF;
  }

  // No Jacobian at all. Do a fix point iteration.
  else if (jacobian == "None")
  {
    preconditioner = "None";
  }

  // This is pretty much debug code. Or rather research code.
  else if (jacobian == "Dumb Finite Difference")
  {
    Teuchos::ParameterList& fdParams = nlParams.sublist("Finite Difference");
    // fdresitemax_ = fdParams.get("itemax", -1);
    double alpha = fdParams.get("alpha", 1.0e-4);
    double beta = fdParams.get("beta", 1.0e-6);
    std::string dt = fdParams.get("Difference Type", "Forward");
    NOX::Epetra::FiniteDifference::DifferenceType dtype = NOX::Epetra::FiniteDifference::Forward;
    if (dt == "Forward")
      dtype = NOX::Epetra::FiniteDifference::Forward;
    else if (dt == "Backward")
      dtype = NOX::Epetra::FiniteDifference::Backward;
    else if (dt == "Centered")
      dtype = NOX::Epetra::FiniteDifference::Centered;
    else
      dserror("unsupported difference type '%s'", dt.c_str());

    FD = Teuchos::rcp(
        new NOX::Epetra::FiniteDifference(printParams, interface, noxSoln, rawGraph_, beta, alpha));
    FD->setDifferenceMethod(dtype);

    iJac = FD;
    J = FD;
  }

  else
  {
    dserror("unsupported Jacobian '%s'", jacobian.c_str());
  }

  // ==================================================================

  // No preconditioning at all.
  if (preconditioner == "None")
  {
    if (Teuchos::is_null(iJac))
    {
      // if no Jacobian has been set this better be the fix point
      // method.
      if (dirParams.get("Method", "Newton") != "User Defined")
      {
        if (Comm().MyPID() == 0)
          utils->out() << "Warning: No Jacobian for solver " << dirParams.get("Method", "Newton")
                       << "\n";
      }
      linSys = Teuchos::rcp(
          new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, interface, noxSoln));
    }
    else
    {
      linSys = Teuchos::rcp(
          new NOX::FSI::LinearSystemGCR(printParams, lsParams, interface, iJac, J, noxSoln));
    }
  }

  else if (preconditioner == "Dump Finite Difference")
  {
    if (lsParams.get("Preconditioner", "None") == "None")
    {
      if (Comm().MyPID() == 0)
        utils->out() << "Warning: Preconditioner turned off in linear solver settings.\n";
    }

    Teuchos::ParameterList& fdParams = nlParams.sublist("Finite Difference");
    // fdresitemax_ = fdParams.get("itemax", -1);
    double alpha = fdParams.get("alpha", 1.0e-4);
    double beta = fdParams.get("beta", 1.0e-6);

    Teuchos::RCP<NOX::Epetra::FiniteDifference> precFD = Teuchos::rcp(
        new NOX::Epetra::FiniteDifference(printParams, interface, noxSoln, rawGraph_, beta, alpha));
    iPrec = precFD;
    M = precFD;

    linSys = Teuchos::rcp(
        new NOX::Epetra::LinearSystemAztecOO(printParams, lsParams, iJac, J, iPrec, M, noxSoln));
  }

  else
  {
    dserror("unsupported preconditioner '%s'", preconditioner.c_str());
  }

  return linSys;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<NOX::StatusTest::Combo> FSI::Partitioned::CreateStatusTest(
    Teuchos::ParameterList& nlParams, Teuchos::RCP<NOX::Epetra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
  Teuchos::RCP<NOX::StatusTest::Combo> converged =
      Teuchos::rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));

  Teuchos::RCP<NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<NOX::StatusTest::FiniteValue> fv = Teuchos::rcp(new NOX::StatusTest::FiniteValue);

  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // setup the real tests
  CreateStatusTest(nlParams, grp, converged);

  return combo;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Partitioned::CreateStatusTest(Teuchos::ParameterList& nlParams,
    Teuchos::RCP<NOX::Epetra::Group> grp, Teuchos::RCP<NOX::StatusTest::Combo> converged)
{
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

  // Teuchos::RCP<NOX::StatusTest::NormWRMS> wrms     = Teuchos::rcp(new
  // NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8)); converged->addStatusTest(wrms);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Partitioned::InitialGuess()
{
  return StructureField()->PredictInterfaceDispnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Partitioned::InterfaceDisp()
{
  // extract displacements
  return StructureField()->ExtractInterfaceDispnp();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Partitioned::InterfaceForce()
{
  // extract forces
  return FluidToStruct(MBFluidField()->ExtractInterfaceForces());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::Partitioned::computeF(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  const char* flags[] = {"Residual", "Jac", "Prec", "FD_Res", "MF_Res", "MF_Jac", "User", NULL};

  Epetra_Time timer(x.Comm());
  const double startTime = timer.WallTime();

  if (Comm().MyPID() == 0)
  {
    utils_->out() << "\n "
                  << "FSI residual calculation"
                  << ".\n";
    if (fillFlag != Residual) utils_->out() << " fillFlag = " << flags[fillFlag] << "\n";
  }

  // we count the number of times the residuum is build
  counter_[fillFlag] += 1;

  if (!x.Map().UniqueGIDs()) dserror("source map not unique");

  if (debugwriter_ != Teuchos::null) debugwriter_->NewIteration();

  // Do the FSI step. The real work is in here.
  FSIOp(x, F, fillFlag);

  if (debugwriter_ != Teuchos::null) debugwriter_->WriteVector("F", F);

  const double endTime = timer.WallTime();
  if (Comm().MyPID() == 0)
    utils_->out() << "\nTime for residual calculation: " << endTime - startTime << " secs\n\n";
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Partitioned::Remeshing() {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Partitioned::FSIOp(const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Partitioned::FluidOp(
    Teuchos::RCP<Epetra_Vector> idisp, const FillType fillFlag)
{
  if (Comm().MyPID() == 0 and utils_->isPrintType(NOX::Utils::OuterIteration))
    utils_->out() << std::endl << "Fluid operator" << std::endl;
  return Teuchos::null;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Partitioned::StructOp(
    Teuchos::RCP<Epetra_Vector> iforce, const FillType fillFlag)
{
  if (Comm().MyPID() == 0 and utils_->isPrintType(NOX::Utils::OuterIteration))
    utils_->out() << std::endl << "Structural operator" << std::endl;
  return Teuchos::null;
}


/*----------------------------------------------------------------------*
 | Calculate interface velocity based on given interface displacement   |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Partitioned::InterfaceVelocity(
    Teuchos::RCP<const Epetra_Vector> idispnp) const
{
  const Teuchos::ParameterList& fsidyn = DRT::Problem::Instance()->FSIDynamicParams();
  Teuchos::RCP<Epetra_Vector> ivel = Teuchos::null;

  if (DRT::INPUT::IntegralValue<int>(fsidyn, "SECONDORDER") == 1)
  {
    ivel = Teuchos::rcp(new Epetra_Vector(*iveln_));
    ivel->Update(2. / Dt(), *idispnp, -2. / Dt(), *idispn_, -1.);
  }
  else
  {
    ivel = Teuchos::rcp(new Epetra_Vector(*idispn_));
    ivel->Update(1. / Dt(), *idispnp, -1. / Dt());
  }
  return ivel;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Partitioned::StructToFluid(Teuchos::RCP<Epetra_Vector> iv)
{
  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  if (matchingnodes_)
  {
    return coupsf.MasterToSlave(iv);
  }
  else
  {
    return coupsfm_->MasterToSlave(iv);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::Partitioned::FluidToStruct(Teuchos::RCP<Epetra_Vector> iv)
{
  const ADAPTER::Coupling& coupsf = StructureFluidCoupling();
  if (matchingnodes_)
  {
    return coupsf.SlaveToMaster(iv);
  }
  else
  {
    // Translate consistent nodal forces to interface loads
    const Teuchos::RCP<Epetra_Vector> ishape = MBFluidField()->IntegrateInterfaceShape();
    const Teuchos::RCP<Epetra_Vector> iforce = Teuchos::rcp(new Epetra_Vector(iv->Map()));

    if (iforce->ReciprocalMultiply(1.0, *ishape, *iv, 0.0)) dserror("ReciprocalMultiply failed");

    return coupsfm_->SlaveToMaster(iforce);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ADAPTER::CouplingMortar& FSI::Partitioned::StructureFluidCouplingMortar() { return *coupsfm_; }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const ADAPTER::CouplingMortar& FSI::Partitioned::StructureFluidCouplingMortar() const
{
  return *coupsfm_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Partitioned::ExtractPreviousInterfaceSolution()
{
  idispn_ = StructureField()->ExtractInterfaceDispn();
  iveln_ = FluidToStruct(MBFluidField()->ExtractInterfaceVeln());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Partitioned::Output()
{
  // call base class version
  FSI::Algorithm::Output();

  switch (DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->FSIDynamicParams(), "COUPALGO"))
  {
    case fsi_iter_stagg_AITKEN_rel_param:
    {
      Teuchos::RCP<NOX::LineSearch::UserDefinedFactory> linesearchfactory =
          noxparameterlist_.sublist("Line Search")
              .get<Teuchos::RCP<NOX::LineSearch::UserDefinedFactory>>(
                  "User Defined Line Search Factory");
      if (linesearchfactory == Teuchos::null)
        dserror("Could not get UserDefinedFactory from noxparameterlist_");
      Teuchos::RCP<NOX::FSI::AitkenFactory> aitkenfactory =
          Teuchos::rcp_dynamic_cast<NOX::FSI::AitkenFactory>(linesearchfactory, true);

      // write aitken relaxation parameter
      MBFluidField()->FluidField()->DiscWriter()->WriteDouble(
          "omega", aitkenfactory->GetAitken()->GetOmega());

      break;
    }  // AITKEN case

  }  // switch-case to be extended for other solver variants if necessary
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::Partitioned::ReadRestart(int step)
{
  // call base class version
  FSI::Algorithm::ReadRestart(step);

  switch (DRT::INPUT::IntegralValue<int>(DRT::Problem::Instance()->FSIDynamicParams(), "COUPALGO"))
  {
    case fsi_iter_stagg_AITKEN_rel_param:
    {
      double omega = -1234.0;

      if (Teuchos::rcp_dynamic_cast<ADAPTER::FluidImmersed>(MBFluidField()) != Teuchos::null)
      {
        IO::DiscretizationReader reader(MBFluidField()->FluidField()->Discretization(), step);
        omega = reader.ReadDouble("omega");
      }
      else if (Teuchos::rcp_dynamic_cast<ADAPTER::FluidAle>(MBFluidField()) != Teuchos::null)
      {
        Teuchos::RCP<ADAPTER::FluidAle> fluidale =
            Teuchos::rcp_dynamic_cast<ADAPTER::FluidAle>(MBFluidField());
        IO::DiscretizationReader reader(
            Teuchos::rcp_const_cast<DRT::Discretization>(fluidale->AleField()->Discretization()),
            step);
        omega = reader.ReadDouble("omega");
      }
      else
        dserror(
            "You want to restart a partitioned FSI scheme with AITKEN relaxation.\n"
            "This is only tested for standard ALE FSI and Immersed FSI.\n"
            "Check the implementation of FSI::Partitioned::ReadRestart.");

      noxparameterlist_.sublist("Line Search").sublist("Aitken").set<int>("restart", step);
      noxparameterlist_.sublist("Line Search")
          .sublist("Aitken")
          .set<double>("restart_omega", omega);

      break;
    }  // AITKEN case

  }  // switch-case to be extended for other solver variants if necessary
}
