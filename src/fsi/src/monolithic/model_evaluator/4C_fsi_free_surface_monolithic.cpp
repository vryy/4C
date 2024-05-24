/*----------------------------------------------------------------------------*/
/*! \file


\level 2

\brief Monolithic solver for free-surface flow
*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
#include "4C_fsi_free_surface_monolithic.hpp"

#include "4C_adapter_ale_fluid.hpp"
#include "4C_ale_utils_mapextractor.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_converter.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fsi_monolithic_linearsystem.hpp"
#include "4C_fsi_nox_newton.hpp"
#include "4C_fsi_overlapprec_fsiamg.hpp"
#include "4C_fsi_statustest.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_control.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_matrixtransform.hpp"
#include "4C_linear_solver_preconditioner_linalg.hpp"

#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicBaseFS::MonolithicBaseFS(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : AlgorithmBase(comm, timeparams)
{
  // ask base algorithm for the fluid time integrator
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluid = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(
      timeparams, GLOBAL::Problem::Instance()->FluidDynamicParams(), "fluid", true));
  fluid_ = fluid->fluid_field();

  // ask base algorithm for the ale time integrator
  Teuchos::RCP<ADAPTER::AleBaseAlgorithm> ale = Teuchos::rcp(
      new ADAPTER::AleBaseAlgorithm(timeparams, GLOBAL::Problem::Instance()->GetDis("ale")));
  ale_ = Teuchos::rcp_dynamic_cast<ADAPTER::AleFluidWrapper>(ale->ale_field());
  if (ale_ == Teuchos::null)
    FOUR_C_THROW("cast from ADAPTER::Ale to ADAPTER::AleFluidWrapper failed");

  coupfa_ = Teuchos::rcp(new CORE::ADAPTER::Coupling());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
CORE::ADAPTER::Coupling& FSI::MonolithicBaseFS::FluidAleCoupling() { return *coupfa_; }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
const CORE::ADAPTER::Coupling& FSI::MonolithicBaseFS::FluidAleCoupling() const { return *coupfa_; }

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBaseFS::read_restart(int step)
{
  fluid_field()->read_restart(step);
  ale_field()->read_restart(step);

  SetTimeStep(fluid_field()->Time(), fluid_field()->Step());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBaseFS::prepare_time_step()
{
  increment_time_and_step();

  PrintHeader();

  fluid_field()->prepare_time_step();
  ale_field()->prepare_time_step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBaseFS::Update()
{
  fluid_field()->Update();
  ale_field()->Update();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicBaseFS::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the Discretizations, which in turn defines the dof number ordering of the
  // Discretizations.
  fluid_field()->Output();
  ale_field()->Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBaseFS::AleToFluid(Teuchos::RCP<Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FSI::MonolithicBaseFS::AleToFluid(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicMainFS::MonolithicMainFS(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : MonolithicBaseFS(comm, timeparams)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicMainFS::Timeloop(
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& interface)
{
  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = NOXParameterList();

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  printParams.set("MyPID", Comm().MyPID());

  // turn on output
  printParams.set("Output Information", 0xffff);

  // Create printing utilities
  utils_ = Teuchos::rcp(new ::NOX::Utils(printParams));

  Teuchos::RCP<std::ofstream> log;
  if (Comm().MyPID() == 0)
  {
    std::string s = GLOBAL::Problem::Instance()->OutputControlFile()->FileName();
    s.append(".iteration");
    log = Teuchos::rcp(new std::ofstream(s.c_str()));
    (*log) << "# num procs      = " << Comm().NumProc() << "\n"
           << "# Method         = " << nlParams.sublist("Direction").get("Method", "Newton") << "\n"
           << "# step | time | time/step | #nliter | res-norm | #liter\n"
           << "#\n";
  }

  Teuchos::Time timer("time step timer");

  while (NotFinished())
  {
    prepare_time_step();

    // start time measurement
    Teuchos::RCP<Teuchos::TimeMonitor> timemonitor =
        Teuchos::rcp(new Teuchos::TimeMonitor(timer, true));

    // calculate initial linear system at current position
    // (no increment)
    // This initializes the field algorithms and creates the first linear
    // systems. And this is the reason we know the initial linear system is
    // there when we create the NOX::Group.
    Evaluate(Teuchos::null);

    // Get initial guess.
    // The initial system is there, so we can happily extract the
    // initial guess. (The Dirichlet conditions are already build in!)
    Teuchos::RCP<Epetra_Vector> initial_guess_v = Teuchos::rcp(new Epetra_Vector(*dof_row_map()));
    initial_guess(initial_guess_v);

    ::NOX::Epetra::Vector noxSoln(initial_guess_v, ::NOX::Epetra::Vector::CreateView);

    // Create the linear system
    Teuchos::RCP<::NOX::Epetra::LinearSystem> linSys =
        create_linear_system(nlParams, noxSoln, utils_);

    // Create the Group
    Teuchos::RCP<NOX::FSI::GroupFS> grp =
        Teuchos::rcp(new NOX::FSI::GroupFS(*this, printParams, interface, noxSoln, linSys));

    // Convergence Tests
    Teuchos::RCP<::NOX::StatusTest::Combo> combo = create_status_test(nlParams, grp);

    // Create the solver
    Teuchos::RCP<::NOX::Solver::Generic> solver = ::NOX::Solver::buildSolver(
        grp, combo, Teuchos::RCP<Teuchos::ParameterList>(&nlParams, false));

    // we know we already have the first linear system calculated
    grp->CaptureSystemState();

    // solve the whole thing
    ::NOX::StatusTest::StatusType status = solver->solve();

    if (status != ::NOX::StatusTest::Converged)
      FOUR_C_THROW("Nonlinear solver failed to converge!");

    // cleanup
    // mat_->Zero();

    // stop time measurement
    timemonitor = Teuchos::null;

    if (Comm().MyPID() == 0)
    {
      (*log) << Step() << "\t" << Time() << " " << timer.totalElapsedTime(true) << " "
             << nlParams.sublist("Output").get("Nonlinear Iterations", 0) << " "
             << nlParams.sublist("Output").get("2-Norm of Residual", 0.) << " "
             << lsParams.sublist("Output").get("Total Number of Linear Iterations", 0);
      (*log) << std::endl;
      lsParams.sublist("Output").set("Total Number of Linear Iterations", 0);
    }

    Update();
    Output();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicMainFS::Evaluate(Teuchos::RCP<const Epetra_Vector> step_increment)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicMainFS::Evaluate");

  Teuchos::RCP<const Epetra_Vector> fx;
  Teuchos::RCP<const Epetra_Vector> ax;

  if (step_increment != Teuchos::null)
  {
    extract_field_vectors(step_increment, fx, ax);
  }

  // Call all elements and assemble rhs and matrices
  // We only need the rhs here because NOX will ask for the rhs
  // only. But the Jacobian is stored internally and will be returned
  // later on without looking at x again!

  Utils()->out() << "\nEvaluate elements\n";

  {
    Teuchos::Time ta("ale", true);
    ale_field()->Evaluate(ax);
    Utils()->out() << "ale      : " << ta.totalElapsedTime(true) << " sec\n";
  }

  // transfer the current ale mesh positions to the fluid field
  Teuchos::RCP<Epetra_Vector> fluiddisp = AleToFluid(ale_field()->Dispnp());
  fluid_field()->apply_mesh_displacement(fluiddisp);

  {
    Teuchos::Time tf("fluid", true);
    fluid_field()->Evaluate(fx);
    Utils()->out() << "fluid    : " << tf.totalElapsedTime(true) << " sec\n";
  }

  Utils()->out() << "\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicMainFS::set_dof_row_maps(
    const std::vector<Teuchos::RCP<const Epetra_Map>>& maps)
{
  Teuchos::RCP<Epetra_Map> fullmap = CORE::LINALG::MultiMapExtractor::MergeMaps(maps);
  blockrowdofmap_.Setup(*fullmap, maps);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicMainFS::set_default_parameters(
    const Teuchos::ParameterList& fsidyn, Teuchos::ParameterList& list)
{
  // monolithic solver settings
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

  // Get the top level parameter list
  Teuchos::ParameterList& nlParams = list;

  nlParams.set<std::string>("Nonlinear Solver", "Line Search Based");
  // nlParams.set("Preconditioner", "None");
  // nlParams.set("Norm abs F", fsimono.get<double>("CONVTOL"));
  nlParams.set("Max Iterations", fsimono.get<int>("ITEMAX"));

  nlParams.set("Norm abs pres", fsimono.get<double>("CONVTOL"));
  nlParams.set("Norm abs vel", fsimono.get<double>("CONVTOL"));
  nlParams.set("Norm abs disp", fsimono.get<double>("CONVTOL"));

  // sublists

  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& solverOptions = nlParams.sublist("Solver Options");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  //   Teuchos::ParameterList& lineSearchParams = nlParams.sublist("Line Search");



  Teuchos::ParameterList& lsParams = newtonParams.sublist("Linear Solver");
  dirParams.set<std::string>("Method", "User Defined");
  Teuchos::RCP<::NOX::Direction::UserDefinedFactory> newtonfactory = Teuchos::rcp(this, false);
  dirParams.set("User Defined Direction Factory", newtonfactory);

  // status tests are expensive, but instructive
  // solverOptions.set<std::string>("Status Test Check Type","Minimal");
  solverOptions.set<std::string>("Status Test Check Type", "Complete");

  // be explicit about linear solver parameters
  lsParams.set<std::string>("Aztec Solver", "GMRES");
  lsParams.set<int>("Size of Krylov Subspace", fsimono.get<int>("KRYLOV_SIZE"));
  lsParams.set<int>("Max Iterations", fsimono.get<int>("KRYLOV_ITEMAX"));
  lsParams.set<std::string>("Preconditioner", "User Defined");
  lsParams.set<int>("Output Frequency", 10);
  lsParams.set<bool>("Output Solver Details", true);

  // adaptive tolerance settings for linear solver
  lsParams.set<double>("base tolerance", fsimono.get<double>("BASETOL"));  // relative tolerance
  lsParams.set<double>(
      "adaptive distance", fsimono.get<double>("ADAPTIVEDIST"));  // adaptive distance
  lsParams.set<INPAR::FSI::Verbosity>(
      "verbosity", CORE::UTILS::IntegralValue<INPAR::FSI::Verbosity>(
                       fsidyn, "VERBOSITY"));  // verbosity level of FSI algorithm
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Direction::Generic> FSI::MonolithicMainFS::buildDirection(
    const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params) const
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
bool FSI::MonolithicMainFS::computeF(
    const Epetra_Vector& x, Epetra_Vector& F, const FillType fillFlag)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicMainFS::computeF");
  Evaluate(Teuchos::rcp(&x, false));
  setup_rhs(F);
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicMainFS::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::MonolithicMainFS::computePreconditioner(
    const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
  return true;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::BlockMonolithicFS::BlockMonolithicFS(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : MonolithicMainFS(comm, timeparams), precondreusecount_(0)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::BlockMonolithicFS::computeJacobian(const Epetra_Vector& x, Epetra_Operator& Jac)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::BlockMonolithicFS::computeJacobian");
  Evaluate(Teuchos::rcp(&x, false));
  CORE::LINALG::BlockSparseMatrixBase& mat =
      Teuchos::dyn_cast<CORE::LINALG::BlockSparseMatrixBase>(Jac);
  setup_system_matrix(mat);
  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FSI::BlockMonolithicFS::computePreconditioner(
    const Epetra_Vector& x, Epetra_Operator& M, Teuchos::ParameterList* precParams)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::BlockMonolithicFS::computePreconditioner");

  if (precondreusecount_ <= 0)
  {
    // Create preconditioner operator. The blocks are already there. This is
    // the perfect place to initialize the block preconditioners.
    SystemMatrix()->SetupPreconditioner();

    const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();
    const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
    precondreusecount_ = fsimono.get<int>("PRECONDREUSE");
  }

  precondreusecount_ -= 1;

  return true;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::BlockMonolithicFS::prepare_time_step()
{
  FSI::MonolithicMainFS::prepare_time_step();

  // new time step, rebuild preconditioner
  precondreusecount_ = 0;
}



/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FSI::MonolithicFS::MonolithicFS(const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : BlockMonolithicFS(comm, timeparams)
{
  aigtransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fmiitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);
  fmgitransform_ = Teuchos::rcp(new CORE::LINALG::MatrixColTransform);

  const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  linearsolverstrategy_ =
      CORE::UTILS::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

  set_default_parameters(fsidyn, NOXParameterList());

  CORE::ADAPTER::Coupling& coupfa = FluidAleCoupling();

  // fluid to ale at the free surface

  icoupfa_ = Teuchos::rcp(new CORE::ADAPTER::Coupling());
  const int ndim = GLOBAL::Problem::Instance()->NDim();
  icoupfa_->setup_condition_coupling(*fluid_field()->Discretization(),
      fluid_field()->Interface()->FSCondMap(), *ale_field()->Discretization(),
      ale_field()->Interface()->FSCondMap(), "FREESURFCoupling", ndim);

  // the fluid-ale coupling always matches
  const Epetra_Map* fluidnodemap = fluid_field()->Discretization()->NodeRowMap();
  const Epetra_Map* alenodemap = ale_field()->Discretization()->NodeRowMap();

  coupfa.setup_coupling(*fluid_field()->Discretization(), *ale_field()->Discretization(),
      *fluidnodemap, *alenodemap, ndim);

  fluid_field()->SetMeshMap(coupfa.MasterDofMap());

  // create combined map

  std::vector<Teuchos::RCP<const Epetra_Map>> vecSpaces;
  vecSpaces.push_back(fluid_field()->dof_row_map());
  vecSpaces.push_back(ale_field()->Interface()->OtherMap());

  set_dof_row_maps(vecSpaces);

  // Use normal matrix for fluid equations but build (splitted) mesh movement
  // linearization (if requested in the input file)
  fluid_field()->UseBlockMatrix(false);

  // build ale system matrix in splitted (at the free surface) system
  ale_field()->CreateSystemMatrix(ale_field()->Interface());

  std::vector<int> pciter;
  std::vector<double> pcomega;
  std::vector<int> spciter;
  std::vector<double> spcomega;
  std::vector<int> fpciter;
  std::vector<double> fpcomega;
  std::vector<int> apciter;
  std::vector<double> apcomega;
  {
    std::string word1;
    std::string word2;
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "PCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsimono, "PCOMEGA"));
      while (pciterstream >> word1) pciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) pcomega.push_back(std::atof(word2.c_str()));
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "STRUCTPCITER"));
      std::istringstream pcomegastream(
          Teuchos::getNumericStringParameter(fsimono, "STRUCTPCOMEGA"));
      while (pciterstream >> word1) spciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) spcomega.push_back(std::atof(word2.c_str()));
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "FLUIDPCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsimono, "FLUIDPCOMEGA"));
      while (pciterstream >> word1) fpciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) fpcomega.push_back(std::atof(word2.c_str()));
    }
    {
      std::istringstream pciterstream(Teuchos::getNumericStringParameter(fsimono, "ALEPCITER"));
      std::istringstream pcomegastream(Teuchos::getNumericStringParameter(fsimono, "ALEPCOMEGA"));
      while (pciterstream >> word1) apciter.push_back(std::atoi(word1.c_str()));
      while (pcomegastream >> word2) apcomega.push_back(std::atof(word2.c_str()));
    }
  }

  // create block system matrix
  switch (linearsolverstrategy_)
  {
    case INPAR::FSI::PreconditionedKrylov:
      systemmatrix_ = Teuchos::rcp(new OverlappingBlockMatrixFS(Extractor(), *fluid_field(),
          *ale_field(), true, CORE::UTILS::IntegralValue<int>(fsimono, "SYMMETRICPRECOND"),
          pcomega[0], pciter[0], fpcomega[0], fpciter[0]));
      break;
    default:
      FOUR_C_THROW("Unsupported type of monolithic free surface solver");
      break;
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFS::setup_rhs(Epetra_Vector& f, bool firstcall)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFS::setup_rhs");

  // see Kue eq. (4.2)
  setup_vector(f, fluid_field()->RHS(), ale_field()->RHS());


  // see Kue eq. (4.21)
  if (firstcall)
  {
    // additional rhs term for ALE equations
    // -dt Aig u(n)
    //
    //    1/dt Delta d(n+1) = theta Delta u(n+1) + u(n)
    //
    // And we are concerned with the u(n) part here.

    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> a = ale_field()->BlockSystemMatrix();
    if (a == Teuchos::null) FOUR_C_THROW("expect ale block matrix");

    // here we extract the free surface submatrices from position 2
    CORE::LINALG::SparseMatrix& aig = a->Matrix(0, 2);

    // extract fluid free surface velocities.
    Teuchos::RCP<Epetra_Vector> fveln = fluid_field()->extract_free_surface_veln();

    Teuchos::RCP<Epetra_Vector> aveln = icoupfa_->MasterToSlave(fveln);

    Teuchos::RCP<Epetra_Vector> rhs = Teuchos::rcp(new Epetra_Vector(aig.RowMap()));
    aig.Apply(*aveln, *rhs);

    rhs->Scale(-1. * Dt());

    Extractor().AddVector(*rhs, 1, f);

    // shape derivatives
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> mmm = fluid_field()->ShapeDerivatives();
    if (mmm != Teuchos::null)
    {
      // here we extract the free surface submatrices from position 2
      CORE::LINALG::SparseMatrix& fmig = mmm->Matrix(0, 2);
      CORE::LINALG::SparseMatrix& fmgg = mmm->Matrix(2, 2);

      rhs = Teuchos::rcp(new Epetra_Vector(fmig.RowMap()));
      fmig.Apply(*fveln, *rhs);
      Teuchos::RCP<Epetra_Vector> veln = fluid_field()->Interface()->InsertOtherVector(rhs);

      rhs = Teuchos::rcp(new Epetra_Vector(fmgg.RowMap()));
      fmgg.Apply(*fveln, *rhs);
      fluid_field()->Interface()->InsertFSCondVector(rhs, veln);

      veln->Scale(-1. * Dt());

      Extractor().AddVector(*veln, 0, f);
    }
  }

  // NOX expects a different sign here.
  f.Scale(-1.);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFS::setup_system_matrix(CORE::LINALG::BlockSparseMatrixBase& mat)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFS::setup_system_matrix");

  // extract Jacobian matrices and put them into composite system
  // matrix W

  // split fluid matrix

  Teuchos::RCP<CORE::LINALG::SparseMatrix> f = fluid_field()->SystemMatrix();

  /*----------------------------------------------------------------------*/

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> a = ale_field()->BlockSystemMatrix();

  if (a == Teuchos::null) FOUR_C_THROW("expect ale block matrix");

  // here we extract the free surface submatrices from position 2
  CORE::LINALG::SparseMatrix& aii = a->Matrix(0, 0);
  CORE::LINALG::SparseMatrix& aig = a->Matrix(0, 2);

  /*----------------------------------------------------------------------*/

  //   double scale     = fluid_field()->ResidualScaling();
  double timescale = fluid_field()->TimeScaling();

  // build block matrix
  // The maps of the block matrix have to match the maps of the blocks we
  // insert here.

  // Uncomplete fluid matrix to be able to deal with slightly defective
  // interface meshes.
  f->UnComplete();

  mat.Assign(0, 0, CORE::LINALG::View, *f);

  (*aigtransform_)(a->FullRowMap(), a->FullColMap(), aig, 1. / timescale,
      CORE::ADAPTER::CouplingSlaveConverter(*icoupfa_), mat.Matrix(1, 0));
  mat.Assign(1, 1, CORE::LINALG::View, aii);

  /*----------------------------------------------------------------------*/
  // add optional fluid linearization with respect to mesh motion block

  Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> mmm = fluid_field()->ShapeDerivatives();
  if (mmm != Teuchos::null)
  {
    // here we extract the free surface submatrices from position 2
    CORE::LINALG::SparseMatrix& fmii = mmm->Matrix(0, 0);
    CORE::LINALG::SparseMatrix& fmig = mmm->Matrix(0, 2);
    CORE::LINALG::SparseMatrix& fmgi = mmm->Matrix(2, 0);
    CORE::LINALG::SparseMatrix& fmgg = mmm->Matrix(2, 2);

    mat.Matrix(0, 0).Add(fmgg, false, 1. / timescale, 1.0);
    mat.Matrix(0, 0).Add(fmig, false, 1. / timescale, 1.0);

    const CORE::ADAPTER::Coupling& coupfa = FluidAleCoupling();

    (*fmgitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmgi, 1.,
        CORE::ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(0, 1), false, false);

    (*fmiitransform_)(mmm->FullRowMap(), mmm->FullColMap(), fmii, 1.,
        CORE::ADAPTER::CouplingMasterConverter(coupfa), mat.Matrix(0, 1), false, true);
  }

  // done. make sure all blocks are filled.
  mat.Complete();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFS::initial_guess(Teuchos::RCP<Epetra_Vector> ig)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFS::initial_guess");

  setup_vector(*ig, fluid_field()->initial_guess(), ale_field()->initial_guess());
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFS::scale_system(CORE::LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& b)
{
  // should we scale the system?
  const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)CORE::UTILS::IntegralValue<int>(fsimono, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    // The matrices are modified here. Do we have to change them back later on?

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(1, 1).EpetraMatrix();

    arowsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    acolsum_ = Teuchos::rcp(new Epetra_Vector(A->RowMap(), false));
    A->InvRowSums(*arowsum_);
    A->InvColSums(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_) or
        mat.Matrix(1, 0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0, 1).EpetraMatrix()->RightScale(*acolsum_))
      FOUR_C_THROW("ale scaling failed");

    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b, 1);

    if (ax->Multiply(1.0, *arowsum_, *ax, 0.0)) FOUR_C_THROW("ale scaling failed");

    Extractor().InsertVector(*ax, 1, b);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFS::unscale_solution(
    CORE::LINALG::BlockSparseMatrixBase& mat, Epetra_Vector& x, Epetra_Vector& b)
{
  const Teuchos::ParameterList& fsidyn = GLOBAL::Problem::Instance()->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
  const bool scaling_infnorm = (bool)CORE::UTILS::IntegralValue<int>(fsimono, "INFNORMSCALING");

  if (scaling_infnorm)
  {
    Teuchos::RCP<Epetra_Vector> ay = Extractor().ExtractVector(x, 1);

    if (ay->Multiply(1.0, *acolsum_, *ay, 0.0)) FOUR_C_THROW("ale scaling failed");

    Extractor().InsertVector(*ay, 1, x);

    Teuchos::RCP<Epetra_Vector> ax = Extractor().ExtractVector(b, 1);

    if (ax->ReciprocalMultiply(1.0, *arowsum_, *ax, 0.0)) FOUR_C_THROW("ale scaling failed");

    Extractor().InsertVector(*ax, 1, b);

    Teuchos::RCP<Epetra_CrsMatrix> A = mat.Matrix(1, 1).EpetraMatrix();
    arowsum_->Reciprocal(*arowsum_);
    acolsum_->Reciprocal(*acolsum_);
    if (A->LeftScale(*arowsum_) or A->RightScale(*acolsum_) or
        mat.Matrix(1, 0).EpetraMatrix()->LeftScale(*arowsum_) or
        mat.Matrix(0, 1).EpetraMatrix()->RightScale(*acolsum_))
      FOUR_C_THROW("ale scaling failed");
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFS::setup_vector(
    Epetra_Vector& f, Teuchos::RCP<const Epetra_Vector> fv, Teuchos::RCP<const Epetra_Vector> av)
{
  // extract the inner dofs of the ale field
  Teuchos::RCP<Epetra_Vector> aov = ale_field()->Interface()->ExtractOtherVector(av);

  Extractor().InsertVector(*fv, 0, f);

  Extractor().InsertVector(*aov, 1, f);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<::NOX::Epetra::LinearSystem> FSI::MonolithicFS::create_linear_system(
    Teuchos::ParameterList& nlParams, ::NOX::Epetra::Vector& noxSoln,
    Teuchos::RCP<::NOX::Utils> utils)
{
  Teuchos::RCP<::NOX::Epetra::LinearSystem> linSys;

  Teuchos::ParameterList& printParams = nlParams.sublist("Printing");
  Teuchos::ParameterList& dirParams = nlParams.sublist("Direction");
  Teuchos::ParameterList& newtonParams = dirParams.sublist("Newton");
  Teuchos::ParameterList* lsParams = nullptr;

  // in case of nonlinCG the linear solver list is somewhere else
  if (dirParams.get("Method", "User Defined") == "User Defined")
    lsParams = &(newtonParams.sublist("Linear Solver"));
  else if (dirParams.get("Method", "User Defined") == "NonlinearCG")
    lsParams = &(dirParams.sublist("Nonlinear CG").sublist("Linear Solver"));
  else
    FOUR_C_THROW("Unknown nonlinear method");

  ::NOX::Epetra::Interface::Jacobian* iJac = this;
  ::NOX::Epetra::Interface::Preconditioner* iPrec = this;
  const Teuchos::RCP<Epetra_Operator> J = systemmatrix_;
  const Teuchos::RCP<Epetra_Operator> M = systemmatrix_;

  switch (linearsolverstrategy_)
  {
    case INPAR::FSI::PreconditionedKrylov:
      linSys = Teuchos::rcp(new FSI::MonolithicLinearSystem(printParams, *lsParams,
          Teuchos::rcp(iJac, false), J, Teuchos::rcp(iPrec, false), M, noxSoln));
      break;
    default:
      FOUR_C_THROW("unsupported linear block solver strategy: %d", linearsolverstrategy_);
      break;
  }

  return linSys;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<::NOX::StatusTest::Combo> FSI::MonolithicFS::create_status_test(
    Teuchos::ParameterList& nlParams, Teuchos::RCP<::NOX::Epetra::Group> grp)
{
  // Create the convergence tests
  Teuchos::RCP<::NOX::StatusTest::Combo> combo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::OR));
  Teuchos::RCP<::NOX::StatusTest::Combo> converged =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::AND));

  Teuchos::RCP<::NOX::StatusTest::MaxIters> maxiters =
      Teuchos::rcp(new ::NOX::StatusTest::MaxIters(nlParams.get("Max Iterations", 100)));
  Teuchos::RCP<::NOX::StatusTest::FiniteValue> fv =
      Teuchos::rcp(new ::NOX::StatusTest::FiniteValue);

  combo->addStatusTest(fv);
  combo->addStatusTest(converged);
  combo->addStatusTest(maxiters);

  // require one solve
  converged->addStatusTest(Teuchos::rcp(new NOX::FSI::MinIters(1)));

  // setup tests for interface

  std::vector<Teuchos::RCP<const Epetra_Map>> interface;
  interface.push_back(fluid_field()->Interface()->FSCondMap());
  interface.push_back(Teuchos::null);
  CORE::LINALG::MultiMapExtractor interfaceextract(*dof_row_map(), interface);

  Teuchos::RCP<::NOX::StatusTest::Combo> interfacecombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> interfaceTest = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "interface", interfaceextract, 0, nlParams.get("Norm abs vel", 1.0e-6),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> interfaceTestUpdate =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("interface update", interfaceextract, 0,
          nlParams.get("Norm abs vel", 1.0e-6), NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(interfaceTest);
  interfacecombo->addStatusTest(interfaceTest);
  // interfacecombo->addStatusTest(interfaceTestUpdate);

  converged->addStatusTest(interfacecombo);

  // setup tests for fluid velocities

  std::vector<Teuchos::RCP<const Epetra_Map>> fluidvel;
  fluidvel.push_back(fluid_field()->InnerVelocityRowMap());
  fluidvel.push_back(Teuchos::null);
  CORE::LINALG::MultiMapExtractor fluidvelextract(*dof_row_map(), fluidvel);

  Teuchos::RCP<::NOX::StatusTest::Combo> fluidvelcombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> innerFluidVel = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "velocity", fluidvelextract, 0, nlParams.get("Norm abs vel", 1.0e-6),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> innerFluidVelUpdate =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("velocity update", fluidvelextract, 0,
          nlParams.get("Norm abs vel", 1.0e-6), NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(innerFluidVel);
  fluidvelcombo->addStatusTest(innerFluidVel);
  // fluidvelcombo->addStatusTest(innerFluidVelUpdate);

  converged->addStatusTest(fluidvelcombo);

  // setup tests for fluid pressure

  std::vector<Teuchos::RCP<const Epetra_Map>> fluidpress;
  fluidpress.push_back(fluid_field()->PressureRowMap());
  fluidpress.push_back(Teuchos::null);
  CORE::LINALG::MultiMapExtractor fluidpressextract(*dof_row_map(), fluidpress);

  Teuchos::RCP<::NOX::StatusTest::Combo> fluidpresscombo =
      Teuchos::rcp(new ::NOX::StatusTest::Combo(::NOX::StatusTest::Combo::OR));

  Teuchos::RCP<NOX::FSI::PartialNormF> fluidPress = Teuchos::rcp(new NOX::FSI::PartialNormF(
      "pressure", fluidpressextract, 0, nlParams.get("Norm abs pres", 1.0e-6),
      ::NOX::Abstract::Vector::TwoNorm, NOX::FSI::PartialNormF::Scaled));
  Teuchos::RCP<NOX::FSI::PartialNormUpdate> fluidPressUpdate =
      Teuchos::rcp(new NOX::FSI::PartialNormUpdate("pressure update", fluidpressextract, 0,
          nlParams.get("Norm abs pres", 1.0e-6), NOX::FSI::PartialNormUpdate::Scaled));

  AddStatusTest(fluidPress);
  fluidpresscombo->addStatusTest(fluidPress);
  // fluidpresscombo->addStatusTest(fluidPressUpdate);

  converged->addStatusTest(fluidpresscombo);

  return combo;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FSI::MonolithicFS::extract_field_vectors(Teuchos::RCP<const Epetra_Vector> x,
    Teuchos::RCP<const Epetra_Vector>& fx, Teuchos::RCP<const Epetra_Vector>& ax)
{
  TEUCHOS_FUNC_TIME_MONITOR("FSI::MonolithicFS::extract_field_vectors");

  fx = Extractor().ExtractVector(x, 0);

  Teuchos::RCP<Epetra_Vector> fcx = fluid_field()->Interface()->ExtractFSCondVector(fx);
  fluid_field()->free_surf_velocity_to_displacement(fcx);

  // process ale unknowns

  Teuchos::RCP<const Epetra_Vector> aox = Extractor().ExtractVector(x, 1);
  Teuchos::RCP<Epetra_Vector> acx = icoupfa_->MasterToSlave(fcx);

  Teuchos::RCP<Epetra_Vector> a = ale_field()->Interface()->InsertOtherVector(aox);
  ale_field()->Interface()->InsertFSCondVector(acx, a);
  ax = a;
}



/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
NOX::FSI::GroupFS::GroupFS(FourC::FSI::MonolithicMainFS& mfsi, Teuchos::ParameterList& printParams,
    const Teuchos::RCP<::NOX::Epetra::Interface::Required>& i, const ::NOX::Epetra::Vector& x,
    const Teuchos::RCP<::NOX::Epetra::LinearSystem>& linSys)
    : ::NOX::Epetra::Group(printParams, i, x, linSys), mfsi_(mfsi)
{
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void NOX::FSI::GroupFS::CaptureSystemState()
{
  // we know we already have the first linear system calculated

  mfsi_.setup_rhs(RHSVector.getEpetraVector(), true);
  mfsi_.setup_system_matrix();

  sharedLinearSystem.getObject(this);
  isValidJacobian = true;
  isValidRHS = true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::FSI::GroupFS::computeF()
{
  ::NOX::Abstract::Group::ReturnType ret = ::NOX::Epetra::Group::computeF();
  if (ret == ::NOX::Abstract::Group::Ok)
  {
    if (not isValidJacobian)
    {
      mfsi_.setup_system_matrix();
      sharedLinearSystem.getObject(this);
      isValidJacobian = true;
    }
  }
  return ret;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::FSI::GroupFS::computeJacobian()
{
  ::NOX::Abstract::Group::ReturnType ret = ::NOX::Epetra::Group::computeJacobian();
  if (ret == ::NOX::Abstract::Group::Ok)
  {
    if (not isValidRHS)
    {
      mfsi_.setup_rhs(RHSVector.getEpetraVector());
      isValidRHS = true;
    }
  }
  return ret;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
::NOX::Abstract::Group::ReturnType NOX::FSI::GroupFS::computeNewton(Teuchos::ParameterList& p)
{
  mfsi_.scale_system(RHSVector.getEpetraVector());
  ::NOX::Abstract::Group::ReturnType status = ::NOX::Epetra::Group::computeNewton(p);
  mfsi_.unscale_solution(NewtonVector.getEpetraVector(), RHSVector.getEpetraVector());
  return status;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::BlockPreconditioningMatrixFS::BlockPreconditioningMatrixFS(
    const CORE::LINALG::MultiMapExtractor& maps, ADAPTER::Fluid& fluid,
    ADAPTER::AleFluidWrapper& ale, int symmetric, double omega, int iterations, double fomega,
    int fiterations, FILE* err)
    : CORE::LINALG::BlockSparseMatrix<CORE::LINALG::DefaultBlockMatrixStrategy>(
          maps, maps, 81, false, true),
      symmetric_(symmetric),
      omega_(omega),
      iterations_(iterations),
      fomega_(fomega),
      fiterations_(fiterations),
      err_(err)
{
  fluidsolver_ = Teuchos::rcp(new CORE::LINALG::Preconditioner(fluid.LinearSolver()));

#ifndef BLOCKMATRIXMERGE
  alesolver_ = Teuchos::rcp(new CORE::LINALG::Preconditioner(ale.LinearSolver()));
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int FSI::BlockPreconditioningMatrixFS::ApplyInverse(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  if (UseTranspose()) FOUR_C_THROW("no transpose preconditioning");

#ifdef BLOCKMATRIXMERGE
  MergeSolve(X, Y);
#else
  SGS(X, Y);
#endif

  return 0;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::BlockPreconditioningMatrixFS::MergeSolve(
    const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
#ifdef BLOCKMATRIXMERGE
  const Epetra_Vector& x = Teuchos::dyn_cast<const Epetra_Vector>(X);
  Epetra_Vector& y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  fluidsolver_->Solve(
      sparse_->EpetraMatrix(), Teuchos::rcp(&y, false), Teuchos::rcp(new Epetra_Vector(x)), true);
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::BlockPreconditioningMatrixFS::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE
  // this is really evil :)
  sparse_ = Merge();
  fluidsolver_->Setup(sparse_->EpetraMatrix());

#else
  const CORE::LINALG::SparseMatrix& fluidInnerOp = Matrix(0, 0);
  const CORE::LINALG::SparseMatrix& aleInnerOp = Matrix(1, 1);

  fluidsolver_->Setup(fluidInnerOp.EpetraMatrix());
  if (constalesolver_ == Teuchos::null) alesolver_->Setup(aleInnerOp.EpetraMatrix());
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::BlockPreconditioningMatrixFS::local_block_richardson(
    Teuchos::RCP<CORE::LINALG::Preconditioner> solver, const CORE::LINALG::SparseMatrix& innerOp,
    Teuchos::RCP<Epetra_Vector> x, Teuchos::RCP<Epetra_Vector> y, Teuchos::RCP<Epetra_Vector> tmpx,
    int iterations, double omega, FILE* err, const Epetra_Comm& comm)
{
  if (iterations > 0)
  {
    y->Scale(omega);
    Teuchos::RCP<Epetra_Vector> tmpy = Teuchos::rcp(new Epetra_Vector(y->Map()));
    if (err != nullptr)
      if (comm.MyPID() == 0) fprintf(err, "    fluid richardson (%d,%f):", iterations, omega);
    for (int i = 0; i < iterations; ++i)
    {
      innerOp.EpetraMatrix()->Multiply(false, *y, *tmpx);
      tmpx->Update(1.0, *x, -1.0);

      if (err != nullptr)
      {
        double n;
        tmpx->Norm2(&n);
        if (comm.MyPID() == 0) fprintf(err, " %e", n);
      }

      solver->Solve(innerOp.EpetraMatrix(), tmpy, tmpx, false);
      y->Update(omega, *tmpy, 1.0);
    }
    if (err != nullptr)
      if (comm.MyPID() == 0) fprintf(err, "\n");
  }
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
FSI::OverlappingBlockMatrixFS::OverlappingBlockMatrixFS(const CORE::LINALG::MultiMapExtractor& maps,
    ADAPTER::Fluid& fluid, ADAPTER::AleFluidWrapper& ale, bool structuresplit, int symmetric,
    double omega, int iterations, double fomega, int fiterations, FILE* err)
    : BlockPreconditioningMatrixFS(
          maps, fluid, ale, symmetric, omega, iterations, fomega, fiterations, err),
      structuresplit_(structuresplit),
      fluid_(fluid),
      ale_(ale)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFS::SetupPreconditioner()
{
#ifdef BLOCKMATRIXMERGE

  // this is really evil :)
  sparse_ = Merge();
  fluidsolver_->Setup(sparse_->EpetraMatrix());

#else
  const CORE::LINALG::SparseMatrix& fluidInnerOp = Matrix(0, 0);
  const CORE::LINALG::SparseMatrix& aleInnerOp = Matrix(1, 1);

  Teuchos::RCP<CORE::LINALG::MapExtractor> fsidofmapex = Teuchos::null;
  Teuchos::RCP<Epetra_Map> irownodes = Teuchos::null;

  fluidsolver_->Setup(fluidInnerOp.EpetraMatrix(), fsidofmapex, fluid_.Discretization(), irownodes,
      structuresplit_);
  if (constalesolver_ == Teuchos::null) alesolver_->Setup(aleInnerOp.EpetraMatrix());
#endif
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FSI::OverlappingBlockMatrixFS::SGS(const Epetra_MultiVector& X, Epetra_MultiVector& Y) const
{
  const CORE::LINALG::SparseMatrix& fluidInnerOp = Matrix(0, 0);
  const CORE::LINALG::SparseMatrix& fluidMeshOp = Matrix(0, 1);
  const CORE::LINALG::SparseMatrix& aleInnerOp = Matrix(1, 1);
  const CORE::LINALG::SparseMatrix& aleBoundOp = Matrix(1, 0);

  // Extract vector blocks
  // RHS

  const Epetra_Vector& x = Teuchos::dyn_cast<const Epetra_Vector>(X);

  // initial guess

  Epetra_Vector& y = Teuchos::dyn_cast<Epetra_Vector>(Y);

  Teuchos::RCP<Epetra_Vector> fy = RangeExtractor().ExtractVector(y, 0);
  Teuchos::RCP<Epetra_Vector> ay = RangeExtractor().ExtractVector(y, 1);

  Teuchos::RCP<Epetra_Vector> fz = Teuchos::rcp(new Epetra_Vector(fy->Map()));
  Teuchos::RCP<Epetra_Vector> az = Teuchos::rcp(new Epetra_Vector(ay->Map()));

  Teuchos::RCP<Epetra_Vector> tmpfx = Teuchos::rcp(new Epetra_Vector(DomainMap(0)));
  Teuchos::RCP<Epetra_Vector> tmpax = Teuchos::rcp(new Epetra_Vector(DomainMap(1)));

  // outer Richardson loop
  for (int run = 0; run < iterations_; ++run)
  {
    Teuchos::RCP<Epetra_Vector> fx = DomainExtractor().ExtractVector(x, 0);
    Teuchos::RCP<Epetra_Vector> ax = DomainExtractor().ExtractVector(x, 1);

    // ----------------------------------------------------------------
    // lower GS
    {
      // Solve ale equations for ay
      if (run > 0)
      {
        aleInnerOp.Multiply(false, *ay, *tmpax);
        ax->Update(-1.0, *tmpax, 1.0);

        aleBoundOp.Multiply(false, *fy, *tmpax);
        ax->Update(-1.0, *tmpax, 1.0);
      }

      alesolver_->Solve(aleInnerOp.EpetraMatrix(), az, ax, true);

      if (run > 0)
      {
        ay->Update(omega_, *az, 1.0);
      }
      else
      {
        ay->Update(omega_, *az, 0.0);
      }
    }

    {
      // Solve fluid equations for fy
      if (run > 0)
      {
        fluidInnerOp.Multiply(false, *fy, *tmpfx);
        fx->Update(-1.0, *tmpfx, 1.0);
      }

      fluidMeshOp.Multiply(false, *ay, *tmpfx);
      fx->Update(-1.0, *tmpfx, 1.0);
      fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(), fz, fx, true);

      local_block_richardson(
          fluidsolver_, fluidInnerOp, fx, fz, tmpfx, fiterations_, fomega_, err_, Comm());

      if (run > 0)
      {
        fy->Update(omega_, *fz, 1.0);
      }
      else
      {
        fy->Update(omega_, *fz, 0.0);
      }
    }

    // ----------------------------------------------------------------
    // the symmetric part of the pc can be skipped

    if (symmetric_)
    {
      fx = DomainExtractor().ExtractVector(x, 0);
      ax = DomainExtractor().ExtractVector(x, 1);

      // ----------------------------------------------------------------
      // upper GS

      {
        fluidInnerOp.Multiply(false, *fy, *tmpfx);
        fx->Update(-1.0, *tmpfx, 1.0);
        fluidMeshOp.Multiply(false, *ay, *tmpfx);
        fx->Update(-1.0, *tmpfx, 1.0);

        fluidsolver_->Solve(fluidInnerOp.EpetraMatrix(), fz, fx, true);

        local_block_richardson(
            fluidsolver_, fluidInnerOp, fx, fz, tmpfx, fiterations_, fomega_, err_, Comm());
        fy->Update(omega_, *fz, 1.0);
      }

      {
        aleInnerOp.Multiply(false, *ay, *tmpax);
        ax->Update(-1.0, *tmpax, 1.0);
        aleBoundOp.Multiply(false, *fy, *tmpax);
        ax->Update(-1.0, *tmpax, 1.0);

        alesolver_->Solve(aleInnerOp.EpetraMatrix(), az, ax, true);
        ay->Update(omega_, *az, 1.0);
      }
    }
  }

  RangeExtractor().InsertVector(*fy, 0, y);
  RangeExtractor().InsertVector(*ay, 1, y);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
const char* FSI::OverlappingBlockMatrixFS::Label() const { return "FSI::OverlappingBlockMatrixFS"; }

FOUR_C_NAMESPACE_CLOSE
