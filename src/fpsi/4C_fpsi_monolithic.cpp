/*----------------------------------------------------------------------*/
/*! \file

\brief General framework for monolithic fpsi solution schemes

\level 3

*/

/*----------------------------------------------------------------------*/
// GENERAL includes
#include "4C_fpsi_monolithic.hpp"

#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fpsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_manipulation.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_method_parameters.hpp"
#include "4C_poroelast_monolithic.hpp"
#include "4C_structure_aux.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

FPSI::MonolithicBase::MonolithicBase(const Epetra_Comm& comm,
    const Teuchos::ParameterList& fpsidynparams, const Teuchos::ParameterList& poroelastdynparams)
    : FpsiBase(comm, fpsidynparams)
{
  // Creation of the subproblems
  // As for general FPSI problems overlapping FSI/FPSI - interfaces can occur, the maps in the
  // different MapExtractors of the Fields are built overlapping! This is necessary to get the right
  // coupling objects, and be able to extract and insert whole FPSI vectors! For building the block
  // matrices of the different fields the following procedure happens to get non overlapping blocks
  // (as wanted):
  // 1. the maps of the FSI and FPSI block are overlapping
  // 2. during the assembling procedure all dofs which belong to the overlapping interface will be
  // assembled into the FSI block (is before the FPSI
  //    condition in the MapExtractors)
  // 3. as there are no entries in of the overlapping interface in the FPSI block, fill_complete()
  // removes them and the FPSI block has no overlapping
  //    entries!

  // --> It's clear that this is not really easily comprehensible, but as the alternatives to this
  // approach requires quite some modifications in 4C
  //     (e.g. extra interface for FSI and FPSI for every Field ...), therefore this version should
  //     be used at the moment.


  // create instance of poroelast subproblem
  poroelast_subproblem_ =
      Teuchos::rcp(new PoroElast::Monolithic(comm, fpsidynparams, Teuchos::null));
  // ask base algorithm for the fluid time integrator
  Global::Problem* problem = Global::Problem::Instance();
  const Teuchos::ParameterList& fluiddynparams = problem->FluidDynamicParams();
  Teuchos::RCP<Adapter::FluidBaseAlgorithm> fluid =
      Teuchos::rcp(new Adapter::FluidBaseAlgorithm(fpsidynparams, fluiddynparams, "fluid", true));
  fluid_subproblem_ = Teuchos::rcp_dynamic_cast<Adapter::FluidFPSI>(fluid->fluid_field());
  // ask base algorithm for the ale time integrator
  Teuchos::RCP<Adapter::AleBaseAlgorithm> ale = Teuchos::rcp(
      new Adapter::AleBaseAlgorithm(fpsidynparams, Global::Problem::Instance()->GetDis("ale")));
  ale_ = Teuchos::rcp_dynamic_cast<Adapter::AleFpsiWrapper>(ale->ale_field());
  if (ale_ == Teuchos::null)
    FOUR_C_THROW("cast from Adapter::Ale to Adapter::AleFpsiWrapper failed");

  coupfa_ = Teuchos::rcp(new Core::Adapter::Coupling());

  coupsf_fsi_ = Teuchos::rcp(new Core::Adapter::Coupling());
  coupsa_fsi_ = Teuchos::rcp(new Core::Adapter::Coupling());
  coupfa_fsi_ = Teuchos::rcp(new Core::Adapter::Coupling());
  icoupfa_fsi_ = Teuchos::rcp(new Core::Adapter::Coupling());

  Teuchos::RCP<FPSI::Utils> FPSI_UTILS = FPSI::Utils::Instance();

  Fluid_PoroFluid_InterfaceMap = FPSI_UTILS->get_fluid_poro_fluid_interface_map();
  PoroFluid_Fluid_InterfaceMap = FPSI_UTILS->get_poro_fluid_fluid_interface_map();

  // build a proxy of the fluid discretization for the ale field
  if (fluid_field()->discretization()->AddDofSet(
          ale_field()->write_access_discretization()->GetDofSetProxy()) != 1)
  {
    FOUR_C_THROW("unexpected numbers of dofsets in fluid field");
  }
  fluid_field()->discretization()->fill_complete(true, false, false);
}  // MonolithicBase
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::read_restart(int step)
{
  poro_field()->read_restart(step);
  fluid_field()->read_restart(step);
  ale_field()->read_restart(step);

  SetTimeStep(fluid_field()->Time(), fluid_field()->Step());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::prepare_time_step()
{
  increment_time_and_step();
  print_header();

  poro_field()->prepare_time_step();
  poro_field()->SetupNewton();
  fluid_field()->prepare_time_step();
  ale_field()->prepare_time_step();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::update()
{
  poro_field()->Update();
  fluid_field()->Update();
  ale_field()->Update();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::prepare_output(bool force_prepare)
{
  poro_field()->prepare_output(force_prepare);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::MonolithicBase::output()
{
  poro_field()->Output();
  fluid_field()->Output();
  ale_field()->Output();
}

/*----------------------------------------------------------------------*/
/*                          Coupling Methods                            */
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::fluid_to_ale(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_->MasterToSlave(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::ale_to_fluid(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_->SlaveToMaster(iv);
}
/// Just in use for problems with FSI-interface ///
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::struct_to_fluid_fsi(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_fsi_->MasterToSlave(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::fluid_to_struct_fsi(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsf_fsi_->SlaveToMaster(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::struct_to_ale_fsi(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_fsi_->MasterToSlave(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::ale_to_struct_fsi(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupsa_fsi_->SlaveToMaster(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::fluid_to_ale_fsi(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_fsi_->MasterToSlave(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::ale_to_fluid_fsi(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return coupfa_fsi_->SlaveToMaster(iv);
}
Teuchos::RCP<Epetra_Vector> FPSI::MonolithicBase::ale_to_fluid_interface_fsi(
    Teuchos::RCP<const Epetra_Vector> iv) const
{
  return icoupfa_fsi_->SlaveToMaster(iv);
}
/// ---------------------------------------------- ///

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<<<<<<<<<  MonolithicBase -> Monolithic  >>>>>>>>>>>>>>>>>>>>>
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

FPSI::Monolithic::Monolithic(const Epetra_Comm& comm, const Teuchos::ParameterList& fpsidynparams,
    const Teuchos::ParameterList& poroelastdynparams)
    : MonolithicBase(comm, fpsidynparams, poroelastdynparams),
      directsolve_(true),
      printscreen_(true),
      printiter_(true),
      timer_("FPSI Monolithic", true),
      isfirsttimestep_(true),
      islinesearch_(false),
      firstcall_(true)
{
  const Teuchos::ParameterList& sdynparams =
      Global::Problem::Instance()->structural_dynamic_params();
  solveradapttol_ = (Core::UTILS::IntegralValue<int>(sdynparams, "ADAPTCONV") == 1);
  solveradaptolbetter_ = (sdynparams.get<double>("ADAPTCONV_BETTER"));

  // hydraulic conductivity (needed for coupling in case of probtype fps3i)
  // is overwritten in class fs3i
  conductivity_ = 0.0;

  // Check if FSI-Interface exists and set flag
  // Will be used to jump over all sections, which are just for FSI condensation procedure required!

  if (fluid_field()->Interface()->Map(FLD::UTILS::MapExtractor::cond_fsi)->NumGlobalElements())
  {
    FSI_Interface_exists_ = true;
    if (comm.MyPID() == 0)
      std::cout << "FPSI Calculation will be performed with FSI - Interface!" << std::endl;
  }
  else
  {
    FSI_Interface_exists_ = false;
    if (comm.MyPID() == 0)
      std::cout << "FPSI Calculation will skip all FSI parts as there is no FSI - Interface!"
                << std::endl;
  }

  // Check for valid predictors
  if (sdynparams.get<std::string>("PREDICT") != "ConstDis")
    FOUR_C_THROW(
        "No Structural Predictor for FPSI implemented at the moment, choose <PREDICT = ConstDis> "
        "in you .dat file! \n --> Or feel free to add the missing terms coming from the predictors "
        "to 4C!");

  const Teuchos::ParameterList& fdynparams = Global::Problem::Instance()->FluidDynamicParams();
  if (fdynparams.get<std::string>("PREDICTOR") != "steady_state")
    FOUR_C_THROW(
        "No Fluid Predictor for FPSI implemented at the moment, choose <PREDICTOR = steady_state> "
        "in you .dat file! \n --> Or feel free to add the missing terms coming from the predictors "
        "to 4C!");

  active_FD_check_ = false;  // to avoid adding RHS of firstiter moreoften!
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::SetupSystem()
{
  const int ndim = Global::Problem::Instance()->NDim();

  Core::Adapter::Coupling& coupfa = fluid_ale_coupling();

  const Epetra_Map* fluidnodemap = fluid_field()->discretization()->NodeRowMap();
  const Epetra_Map* alenodemap = ale_field()->discretization()->NodeRowMap();

  coupfa.setup_coupling(*fluid_field()->discretization(), *ale_field()->discretization(),
      *fluidnodemap, *alenodemap, ndim, false);
  fluid_field()->SetMeshMap(coupfa.MasterDofMap());

  if (FSI_Interface_exists_) SetupSystem_FSI();

  // Setup the FPSI Coupling Adapter
  FPSICoupl() = Teuchos::rcp(new FPSI::FpsiCoupling(poroelast_subproblem_, fluid_subproblem_, ale_,
      Fluid_PoroFluid_InterfaceMap, PoroFluid_Fluid_InterfaceMap));
  FPSICoupl()->SetConductivity(conductivity_);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::SetupSystem_FSI()
{
  // create local FPSI MapExtractors as the MapExtractors of the single Fields are nonoverlapping
  // and therefore not the whole FPSI Interface Map is availlable for FSI, FPSI Interface overlap!!!

  // right now we use matching meshes at the interface

  const int ndim = Global::Problem::Instance()->NDim();

  Core::Adapter::Coupling& coupsf_fsi = structure_fluid_coupling_fsi();
  Core::Adapter::Coupling& coupsa_fsi = structure_ale_coupling_fsi();
  Core::Adapter::Coupling& icoupfa_fsi = interface_fluid_ale_coupling_fsi();

  // structure to fluid
  coupsf_fsi.setup_condition_coupling(*poro_field()->structure_field()->discretization(),
      poro_field()->structure_field()->Interface()->FSICondMap(), *fluid_field()->discretization(),
      fluid_field()->Interface()->FSICondMap(), "FSICoupling", ndim);
  // structure to ale
  coupsa_fsi.setup_condition_coupling(*poro_field()->structure_field()->discretization(),
      poro_field()->structure_field()->Interface()->FSICondMap(), *ale_field()->discretization(),
      ale_field()->Interface()->FSICondMap(), "FSICoupling", ndim);

  // fluid to ale at the interface

  icoupfa_fsi.setup_condition_coupling(*fluid_field()->discretization(),
      fluid_field()->Interface()->FSICondMap(), *ale_field()->discretization(),
      ale_field()->Interface()->FSICondMap(), "FSICoupling", ndim);


  // In the following we assume that both couplings find the same dof
  // map at the structural side. This enables us to use just one
  // interface dof map for all fields and have just one transfer
  // operator from the interface map to the full field map.
  if (not coupsf_fsi.MasterDofMap()->SameAs(*coupsa_fsi.MasterDofMap()))
    FOUR_C_THROW("fsi structure interface dof maps do not match");

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::Timeloop()
{
  PrepareTimeloop();

  while (NotFinished())  // while step < maxsteps and time < maxtime
  {
    prepare_time_step();
    SetupNewton();
    TimeStep();
    constexpr bool force_prepare = false;
    prepare_output(force_prepare);
    update();
    output();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::PrepareTimeloop()
{
  // check if maps were destroyed before entring the timeloop
  Extractor().check_for_valid_map_extractor();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FPSI::Monolithic::TimeStep()
{
  //////////////////////////////////////////////////////////////////////////////////////////////
  ///////////////                                                                ///////////////
  ///////////////                                 LOOP                           ///////////////
  ///////////////                                                                ///////////////
  //////////////////////////////////////////////////////////////////////////////////////////////
  while ((((not Converged()) and (iter_ <= maximumiterations_)) or (iter_ <= minimumiterations_)) or
         islinesearch_ == true)
  {
    // start time measurement
    timer_.reset();
    Teuchos::Time timer("FPSI Time Step", true);
    evaluate(iterinc_);
    // create full monolithic FPSI right-hand-side vector
    // moved to evaluate()

    // create full monolithic FPSI tangent stiffness matrix and check if it is filled
    setup_system_matrix();

    if (not systemmatrix_->Filled())
    {
      FOUR_C_THROW("Effective tangent matrix must be filled here !");
    }
    // (Newton-ready) residual with blanked Dirichlet DOFs (see adapter_timint!)
    // is done in prepare_system_for_newton_solve() within evaluate(iterinc_)
    linear_solve();
    // build norms
    build_convergence_norms();

    // print stuff
    if (islinesearch_ == false) print_newton_iter();

    // reset solver tolerance
    solver_->ResetTolerance();

    // increment equilibrium loop index
    if (islinesearch_ == false)
    {
      iter_ += 1;
      poro_field()->IncrementPoroIter();
    }
  }  // end loop

  // correct iteration counter
  iter_ -= 1;

  // test whether max iterations was hit
  if ((Converged()) and (Comm().MyPID() == 0))
  {
    if (linesearch_counter > 0.5)
      std::cout << "            Evaluation of residual with scaled increment yields: " << normofrhs_
                << std::endl;
    islinesearch_ = false;
    linesearch_counter = 0.;
  }
  else if (iter_ >= maximumiterations_)
  {
    FOUR_C_THROW("Newton found no convergence in %d iterations", iter_);
  }

  poro_field()->recover_lagrange_multiplier_after_time_step();

  // recover Lagrange multiplier \lambda_{\Gamma} at the interface at the end of each time step
  // (i.e. condensed traction/forces onto the structure) needed for rhs in next time step
  if (FSI_Interface_exists_)
    recover_lagrange_multiplier();  // LagrangeMultiplier of the FSI interface!
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/

void FPSI::Monolithic::evaluate(Teuchos::RCP<const Epetra_Vector> stepinc)
{
  TEUCHOS_FUNC_TIME_MONITOR("FPSI::Monolithic::Evaluate");

  if (linesearch_ and islinesearch_ == false)
  {
    fluid_field()->discretization()->ClearState();
    linesearch_counter = 0.;
    fluid_field()->discretization()->set_state(0, "dispnp", fluid_field()->Dispnp());
    meshdispold_ = ale_to_fluid(ale_field()->Dispnp());
    porointerfacedisplacementsold_ = FPSICoupl()->iPorostructToAle(
        poro_field()->structure_field()->extract_interface_dispnp(true));
  }


  Teuchos::RCP<const Epetra_Vector> sx;
  Teuchos::RCP<const Epetra_Vector> fx;
  Teuchos::RCP<const Epetra_Vector> pfx;
  Teuchos::RCP<const Epetra_Vector> ax;

  if (stepinc != Teuchos::null)
  {
    extract_field_vectors(stepinc, sx, pfx, fx, ax, (iter_ == 1 and !active_FD_check_));
  }
  else
  {
    FOUR_C_THROW("No existing increment vector !");
  }

  poro_field()->update_state_incrementally(sx, pfx);
  poro_field()->EvaluateFields(Teuchos::null);
  poro_field()->setup_system_matrix();

  Teuchos::RCP<Epetra_Vector> porointerfacedisplacements_FPSI = FPSICoupl()->iPorostructToAle(
      poro_field()->structure_field()->extract_interface_dispnp(true));
  ale_field()->apply_interface_displacements(porointerfacedisplacements_FPSI);

  if (FSI_Interface_exists_)
  {
    Teuchos::RCP<Epetra_Vector> porointerfacedisplacements_FSI =
        struct_to_ale_fsi(poro_field()->structure_field()->extract_interface_dispnp(false));
    ale_field()->apply_fsi_interface_displacements(porointerfacedisplacements_FSI);
  }

  ale_field()->WriteAccessDispnp()->Update(
      1.0, *ax, 1.0);  // displacement increments on the interfaces are zero!!!
  ale_field()->evaluate(Teuchos::null);

  Teuchos::RCP<const Epetra_Vector> aledisplacements = ale_to_fluid(ale_field()->Dispnp());
  fluid_field()->apply_mesh_displacement(aledisplacements);

  fluid_field()->UpdateNewton(fx);

  fluid_field()->evaluate(Teuchos::null);

  // Evaluate FPSI Coupling Matrixes and RHS
  FPSICoupl()->evaluate_coupling_matrixes_rhs();

  setup_rhs(iter_ == 1);

}  // Evaluate

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::TestResults(const Epetra_Comm& comm)
{
  Global::Problem::Instance()->AddFieldTest(poro_field()->structure_field()->CreateFieldTest());
  Global::Problem::Instance()->AddFieldTest(poro_field()->fluid_field()->CreateFieldTest());
  Global::Problem::Instance()->AddFieldTest(fluid_field()->CreateFieldTest());
  Global::Problem::Instance()->TestAll(comm);
}

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<              Methods concerning          >>>>>>>>>>>>>>>>>>>>
//<<<<<<<<<<<<<<                    solver                >>>>>>>>>>>>>>>>>>>>
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

void FPSI::Monolithic::SetupSolver()
{
  const Teuchos::ParameterList& fpsidynamicparams =
      Global::Problem::Instance()->FPSIDynamicParams();

  const int linsolvernumber = fpsidynamicparams.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "No linear solver defined for FPSI problem. Please set LINEAR_SOLVER in FPSI DYNAMIC to a "
        "valid number !");

  const Teuchos::ParameterList& solverparams =
      Global::Problem::Instance()->SolverParams(linsolvernumber);
  const auto solvertype =
      Teuchos::getIntegralValue<Core::LinearSolver::SolverType>(solverparams, "SOLVER");

  directsolve_ = (solvertype == Core::LinearSolver::SolverType::umfpack or
                  solvertype == Core::LinearSolver::SolverType::superlu);

  if (directsolve_)
    solver_ = Teuchos::rcp(new Core::LinAlg::Solver(solverparams, Comm(),
        Global::Problem::Instance()->solver_params_callback(),
        Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
            Global::Problem::Instance()->IOParams(), "VERBOSITY")));
  else
    // create a linear solver
    create_linear_solver();

  // Get the parameters for the Newton iteration
  maximumiterations_ = fpsidynamicparams.get<int>("ITEMAX");
  minimumiterations_ = fpsidynamicparams.get<int>("ITEMIN");
  normtypeinc_ =
      Core::UTILS::IntegralValue<Inpar::FPSI::ConvergenceNorm>(fpsidynamicparams, "NORM_INC");
  normtypefres_ =
      Core::UTILS::IntegralValue<Inpar::FPSI::ConvergenceNorm>(fpsidynamicparams, "NORM_RESF");
  combinedconvergence_ =
      Core::UTILS::IntegralValue<Inpar::FPSI::BinaryOp>(fpsidynamicparams, "NORMCOMBI_RESFINC");

  {
    std::istringstream tolresstream(
        Teuchos::getNumericStringParameter(fpsidynamicparams, "RESTOL"));
    std::string word;
    while (tolresstream >> word)
    {
      toleranceresidualforceslist_.push_back(std::atof(word.c_str()));
    }
    toleranceresidualforces_ = toleranceresidualforceslist_[0];

    std::istringstream tolincstream(
        Teuchos::getNumericStringParameter(fpsidynamicparams, "INCTOL"));
    while (tolincstream >> word)
    {
      toleranceiterinclist_.push_back(std::atof(word.c_str()));
    }
    toleranceiterinc_ = toleranceiterinclist_[0];
  }

  Global::Problem* problem = Global::Problem::Instance();
  const Teuchos::ParameterList& fpsidynparams = problem->FPSIDynamicParams();
  linesearch_ = Core::UTILS::IntegralValue<int>(fpsidynparams, "LineSearch");
  if (linesearch_ == 1)
    FOUR_C_THROW(
        "Parameter 'LineSearch' is set to 'Yes' in the FPSI Dynamic section in your input-file.  \n"
        "Though the framework for a line search algorithm is implemented in fpsi_monolithic.cpp, \n"
        "a proper routine to reset the participating single fields is still required. In Chuck's \n"
        "experimental 4C this was solved by performing an evaluate with the negative increment.\n"
        "However this has not yet been committed.\n");
  linesearch_counter = 0.;

  return;
}

/*----------------------------------------------------------------------*
 | create linear solver                                   vuong 08/15 |
 *----------------------------------------------------------------------*/
void FPSI::Monolithic::create_linear_solver()
{
  // get dynamic section
  const Teuchos::ParameterList& fpsidyn = Global::Problem::Instance()->FPSIDynamicParams();

  // get the linear solver number
  const int linsolvernumber = fpsidyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "No linear solver defined for FPSI problem. Please set LINEAR_SOLVER in FPSI DYNAMIC to a "
        "valid number !");

  // get parameter list of structural dynamics
  const Teuchos::ParameterList& sdyn = Global::Problem::Instance()->structural_dynamic_params();
  // use solver blocks for structure
  // get the solver number used for structural solver
  const int slinsolvernumber = sdyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (slinsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for structural field. Please set LINEAR_SOLVER in STRUCTURAL "
        "DYNAMIC to a valid number!");

  // get parameter list of fluid dynamics
  const Teuchos::ParameterList& fdyn = Global::Problem::Instance()->FluidDynamicParams();
  // use solver blocks for fluid
  // get the solver number used for fluid solver
  const int flinsolvernumber = fdyn.get<int>("LINEAR_SOLVER");
  // check if the fluid solver has a valid solver number
  if (flinsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for fluid field. Please set LINEAR_SOLVER in FLUID DYNAMIC to a "
        "valid number!");

  // get parameter list of structural dynamics
  const Teuchos::ParameterList& aledyn = Global::Problem::Instance()->AleDynamicParams();
  // use solver blocks for structure
  // get the solver number used for structural solver
  const int alinsolvernumber = aledyn.get<int>("LINEAR_SOLVER");
  // check if the structural solver has a valid solver number
  if (alinsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for ALE field. Please set LINEAR_SOLVER in ALE DYNAMIC to a "
        "valid number!");

  // get solver parameter list of linear Poroelasticity solver
  const Teuchos::ParameterList& fpsisolverparams =
      Global::Problem::Instance()->SolverParams(linsolvernumber);

  const auto solvertype =
      Teuchos::getIntegralValue<Core::LinearSolver::SolverType>(fpsisolverparams, "SOLVER");

  if (solvertype != Core::LinearSolver::SolverType::belos)
  {
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    std::cout << " Note: the BGS2x2 preconditioner now " << std::endl;
    std::cout << " uses the structural solver and fluid solver blocks" << std::endl;
    std::cout << " for building the internal inverses" << std::endl;
    std::cout << " Remove the old BGS PRECONDITIONER BLOCK entries " << std::endl;
    std::cout << " in the dat files!" << std::endl;
    std::cout << "!!!!!!!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!!!!!!" << std::endl;
    FOUR_C_THROW("Iterative solver expected");
  }
  const auto azprectype =
      Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(fpsisolverparams, "AZPREC");

  // plausibility check
  switch (azprectype)
  {
    case Core::LinearSolver::PreconditionerType::multigrid_nxn:
    {
      // no plausibility checks here
      // if you forget to declare an xml file you will get an error message anyway
    }
    break;
    default:
      FOUR_C_THROW("AMGnxn preconditioner expected");
      break;
  }

  solver_ = Teuchos::rcp(new Core::LinAlg::Solver(fpsisolverparams, Comm(),
      Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY")));

  // use solver blocks for structure and fluid
  const Teuchos::ParameterList& ssolverparams =
      Global::Problem::Instance()->SolverParams(slinsolvernumber);
  const Teuchos::ParameterList& fsolverparams =
      Global::Problem::Instance()->SolverParams(flinsolvernumber);
  const Teuchos::ParameterList& asolverparams =
      Global::Problem::Instance()->SolverParams(alinsolvernumber);

  // for now, use same solver parameters for poro fluid and free fluid

  // poro/structure
  solver_->put_solver_params_to_sub_params("Inverse1", ssolverparams,
      Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY"));
  // poro fluid
  solver_->put_solver_params_to_sub_params("Inverse2", fsolverparams,
      Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY"));
  // fluid
  solver_->put_solver_params_to_sub_params("Inverse3", fsolverparams,
      Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY"));
  // ale
  solver_->put_solver_params_to_sub_params("Inverse4", asolverparams,
      Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY"));

  // prescribe rigid body modes
  poro_field()->structure_field()->discretization()->compute_null_space_if_necessary(
      solver_->Params().sublist("Inverse1"));
  poro_field()->fluid_field()->discretization()->compute_null_space_if_necessary(
      solver_->Params().sublist("Inverse2"));
  fluid_field()->discretization()->compute_null_space_if_necessary(
      solver_->Params().sublist("Inverse3"));
  ale_field()->write_access_discretization()->compute_null_space_if_necessary(
      solver_->Params().sublist("Inverse4"));

  // fixing length of Inverse1 nullspace (solver/preconditioner ML)
  {
    std::string inv = "Inverse1";
    const Epetra_Map& oldmap = *(Global::Problem::Instance()->GetDis("structure")->dof_row_map());
    const Epetra_Map& newmap =
        systemmatrix_->Matrix(structure_block_, structure_block_).EpetraMatrix()->RowMap();
    Core::LinearSolver::Parameters::FixNullSpace(
        inv.data(), oldmap, newmap, solver_->Params().sublist("Inverse1"));
  }
  // fixing length of Inverse2 nullspace (solver/preconditioner ML)
  {
    std::string inv = "Inverse2";
    const Epetra_Map& oldmap = *(poro_field()->fluid_field()->dof_row_map());
    const Epetra_Map& newmap =
        systemmatrix_->Matrix(porofluid_block_, porofluid_block_).EpetraMatrix()->RowMap();
    Core::LinearSolver::Parameters::FixNullSpace(
        inv.data(), oldmap, newmap, solver_->Params().sublist("Inverse2"));
  }
  // fixing length of Inverse3 nullspace (solver/preconditioner ML)
  {
    std::string inv = "Inverse3";
    const Epetra_Map& oldmap = *(fluid_field()->dof_row_map());
    const Epetra_Map& newmap =
        systemmatrix_->Matrix(fluid_block_, fluid_block_).EpetraMatrix()->RowMap();
    Core::LinearSolver::Parameters::FixNullSpace(
        inv.data(), oldmap, newmap, solver_->Params().sublist("Inverse3"));
  }
  // fixing length of Inverse4 nullspace (solver/preconditioner ML)
  {
    std::string inv = "Inverse4";
    const Epetra_Map& oldmap = *(ale_field()->dof_row_map());
    const Epetra_Map& newmap =
        systemmatrix_->Matrix(ale_i_block_, ale_i_block_).EpetraMatrix()->RowMap();
    Core::LinearSolver::Parameters::FixNullSpace(
        inv.data(), oldmap, newmap, solver_->Params().sublist("Inverse4"));
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::linear_solve()
{
  Core::LinAlg::SolverParams solver_params;
  if (solveradapttol_ and (iter_ > 1))
  {
    solver_params.nonlin_tolerance = normofrhs_;
    solver_params.nonlin_residual = toleranceresidualforces_;
    solver_params.lin_tol_better = solveradaptolbetter_;
  }
  Global::Problem* problem = Global::Problem::Instance();
  const Teuchos::ParameterList& fpsidynparams = problem->FPSIDynamicParams();
  if (Teuchos::getIntegralValue<int>(fpsidynparams, "FDCheck"))
  {
    FPSIFDCheck();
  }

  iterinc_->PutScalar(0.0);  // Useful? depends on solver and more
  poro_field()->ClearPoroIterinc();

  if (directsolve_)
  {
    Teuchos::RCP<Core::LinAlg::SparseMatrix> sparse = systemmatrix_->Merge();

    if (FSI_Interface_exists_)
    {
      // remove entries in condensed dofs from matrix and rhs...
      Core::LinAlg::apply_dirichlet_to_system(
          *sparse, *iterinc_, *rhs_, *zeros_, *fluid_field()->Interface()->FSICondMap());
    }

    Core::LinAlg::apply_dirichlet_to_system(
        *sparse, *iterinc_, *rhs_, *zeros_, *combined_dbc_map());

    // line search
    if (linesearch_)
      LineSearch(sparse);
    else
    {
      // standard solver call
      solver_params.refactor = true;
      solver_params.reset = iter_ == 1;
      solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, solver_params);
    }
  }
  else
  {
    if (FSI_Interface_exists_)
    {
      // remove entries in condensed dofs from matrix and rhs...
      Core::LinAlg::apply_dirichlet_to_system(
          *systemmatrix_, *iterinc_, *rhs_, *zeros_, *fluid_field()->Interface()->FSICondMap());
    }

    Core::LinAlg::apply_dirichlet_to_system(
        *systemmatrix_, *iterinc_, *rhs_, *zeros_, *combined_dbc_map());

    // standard solver call
    Core::LinAlg::SolverParams solver_params;
    solver_params.refactor = true;
    solver_params.reset = iter_ == 1;
    solver_->Solve(systemmatrix_->EpetraOperator(), iterinc_, rhs_, solver_params);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::LineSearch(Teuchos::RCP<Core::LinAlg::SparseMatrix>& sparse)
{
  // Note: the line search code seems to be experimental and is not
  // working properly (perhaps just a sign is wrong somewhere ...)
  // I would not recommend to use it without thorough testing
  // vuong 08/15

  if (iter_ > 1)
  {
    rhs_->Norm2(&normofrhs_);
    if (normofrhs_ - normofrhsold_ > 1e-13)
    {
      if (linesearch_counter > 0.5)
        std::cout << "            Evaluation of residual with bisected increment yields: "
                  << normofrhs_ << std::endl;

      islinesearch_ = true;
      iterinc_->Update(pow(0.5, (linesearch_counter)), *iterincold_, 0.0);
      linesearch_counter = linesearch_counter + 1.0;
      std::cout << "linesearch_ : " << std::setprecision(1)
                << static_cast<int>(linesearch_counter + 0.5) << " iterinc_ multiplied by "
                << std::setprecision(4) << pow(0.5, linesearch_counter)
                << "   residual = " << normofrhs_ << " > " << normofrhsold_ << std::endl;

      // substract the old interinc_ from all fields (undo the update)
      Teuchos::RCP<Epetra_Vector> sx;
      Teuchos::RCP<Epetra_Vector> pfx;
      Teuchos::RCP<Epetra_Vector> fx;
      Teuchos::RCP<const Epetra_Vector> constsx;
      Teuchos::RCP<const Epetra_Vector> constfpx;
      Teuchos::RCP<const Epetra_Vector> constfx;
      Teuchos::RCP<const Epetra_Vector> ax;

      sx = Teuchos::rcp(new Epetra_Vector(*poro_field()->structure_field()->dof_row_map(), true));
      pfx = Teuchos::rcp(new Epetra_Vector(*poro_field()->fluid_field()->dof_row_map(), true));
      fx = Teuchos::rcp(new Epetra_Vector(*fluid_field()->dof_row_map(), true));

      extract_field_vectors(iterinc_, constsx, constfpx, constfx, ax, iter_ == 1);
      iterinc_->Norm2(&normofiterinc_);
      std::cout << "            Norm of step back: " << normofiterinc_ << std::endl;
      // poro_field()  ->ResetNewton(sx);
      // fluid_field() ->ResetNewton(fx);
      sx->Update(1.0, *constsx, 0.0);
      pfx->Update(1.0, *constfpx, 0.0);
      fx->Update(1.0, *constfx, 0.0);
      sx->Scale(-1.0);
      pfx->Scale(-1.0);
      fx->Scale(-1.0);

      poro_field()->update_state_incrementally(sx, pfx);
      poro_field()->EvaluateFields(Teuchos::null);
      poro_field()->setup_system_matrix();

      fluid_field()->UpdateNewton(Teuchos::rcp_dynamic_cast<const Epetra_Vector>(fx));
      // ale_field()   ->ResetNewton(ax);

      fluid_field()->apply_mesh_displacement(meshdispold_);
      ale_field()->apply_interface_displacements(porointerfacedisplacementsold_);


      // set iterinc_ to a fraction of the old iterinc_
      iterinc_->Update(pow(0.5, linesearch_counter), *iterincold_, 0.0);
      iterinc_->Norm2(&normofiterinc_);
      std::cout << "            Norm of old increment: " << normofiterincold_
                << "  Norm of bisected increment: " << normofiterinc_ << std::endl;
    }
    else
    {
      islinesearch_ = false;
      if (linesearch_counter > 0.5)
        std::cout << "            Evaluation of residual with bisected increment yields: "
                  << normofrhs_ << std::endl;
      linesearch_counter = 0.0;
    }
  }

  // prepare linesearch_
  // copy the old iterinc_ before new solve
  if (linesearch_ and islinesearch_ == false)
  {
    rhsold_ = Core::LinAlg::CreateVector(*dof_row_map(), true);
    rhsold_->Update(1.0, *rhs_, 0.0);
    rhsold_->Norm2(&normofrhsold_);
    if (abs(normofrhs_ - normofrhsold_) > 1.e-12 and iter_ > 1)
      FOUR_C_THROW(" wrong copy of rhs_ ");
  }
  // end prepare linesearch_

  // standard solver call
  if (islinesearch_ == false)
  {
    Core::LinAlg::SolverParams solver_params;
    solver_params.refactor = true;
    solver_params.reset = iter_ == 1;
    solver_->Solve(sparse->EpetraOperator(), iterinc_, rhs_, solver_params);
  }

  if (islinesearch_ == false)
  {
    // check whether iterinc_ points in right direction
    Teuchos::RCP<Epetra_Vector> tempvec = Core::LinAlg::CreateVector(*dof_row_map(), true);
    sparse->Multiply(true, *rhs_, *tempvec);
    double climb = 0.0;
    tempvec->Dot(*iterinc_, &climb);
    climb = -climb;

    if (climb > 0.0)
    {
      std::cout << "########################################################################"
                << std::endl;
      std::cout << "##                                                                    ##"
                << std::endl;
      std::cout << "## WARNING: A*x-b=0 ; A^T*b*x > 0 ; increment vector multiplied by -1 ##"
                << std::endl;
      std::cout << "##                                                                    ##"
                << std::endl;
      std::cout << "##                       Value = " << std::setprecision(9) << climb << "    ##"
                << std::endl;
      std::cout << "##                                                                    ##"
                << std::endl;
      std::cout << "########################################################################"
                << std::endl;
      iterinc_->Update(-1.0, *iterinc_, 0.0);
    }
  }

  if (linesearch_ and islinesearch_ == false)
  {
    iterincold_ = Core::LinAlg::CreateVector(*dof_row_map(), true);
    iterincold_->Update(1.0, *iterinc_, 0.0);
    iterincold_->Norm2(&normofiterincold_);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_Map> FPSI::Monolithic::combined_dbc_map()
{
  const Teuchos::RCP<const Epetra_Map> scondmap =
      poro_field()->structure_field()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> pfcondmap =
      poro_field()->fluid_field()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> fcondmap = fluid_field()->GetDBCMapExtractor()->CondMap();
  const Teuchos::RCP<const Epetra_Map> acondmap = ale_field()->GetDBCMapExtractor()->CondMap();
  Teuchos::RCP<Epetra_Map> tempmap = Core::LinAlg::MergeMap(scondmap, pfcondmap, false);
  Teuchos::RCP<Epetra_Map> condmap_0 = Core::LinAlg::MergeMap(tempmap, fcondmap, false);
  Teuchos::RCP<Epetra_Map> condmap = Core::LinAlg::MergeMap(condmap_0, acondmap, false);

  return condmap;
}  // combined_dbc_map()

//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
//<<<<<<<<<<<<<<                 Newton Loop              >>>>>>>>>>>>>>>>>>>>
//<<<<<<<<<<<<<<                   Methods                >>>>>>>>>>>>>>>>>>>>
//<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
bool FPSI::Monolithic::Converged()
{
  // check for single norms
  bool convinc = false;
  bool convfres = false;

  // residual increments
  switch (normtypeinc_)
  {
    case Inpar::FPSI::absoluteconvergencenorm:
      convinc = normofiterinc_ < toleranceiterinc_;
      break;
    case Inpar::FPSI::absoluteconvergencenorm_sys_split:
      FOUR_C_THROW(
          "Check for convergence of primary variables with type <absolut sys split> not "
          "implemented yet!");
      break;
    case Inpar::FPSI::relativconvergencenorm_sys:
      // increment convergence is checked relative to the average dofs of the different fields
      std::cout << " |ps| " << norm1_ps_ << " |pfv| " << norm1_pfv_ << " |pfp| " << norm1_pfp_
                << " |fv| " << norm1_fv_ << " |fp| " << norm1_fp_ << " |a| " << norm1_a_
                << std::endl;
      convinc =
          ((normofiterincporostruct_ / norm1_ps_ * sqrtnps_ < toleranceiterinc_) and  // poro struct
              (normofiterincporofluidvelocity_ / norm1_pfv_ * sqrtnpfv_ <
                  toleranceiterinc_) and  // poro fluid velocity
              (normofiterincporofluidpressure_ / norm1_pfp_ * sqrtnpfp_ <
                  toleranceiterinc_) and  // poro fluid pressure
              (normofiterincale_ / norm1_a_ * sqrtna_ < toleranceiterinc_) and  // ale
              (normofiterincfluidvelocity_ / norm1_fv_ * sqrtnfv_ <
                  toleranceiterinc_) and  // fluid velocity
              (normofiterincfluidpressure_ / norm1_fp_ * sqrtnfp_ <
                  toleranceiterinc_)  // fluid pressure
          );
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of primary variables for any reason :-p !");
      break;
  }
  // residual forces
  switch (normtypefres_)
  {
    case Inpar::FPSI::absoluteconvergencenorm:
      convfres = normofrhs_ < toleranceresidualforces_;
      break;
    case Inpar::FPSI::absoluteconvergencenorm_sys_split:
      convfres = ((normrhsporostruct_ / sqrtnps_ < toleranceresidualforceslist_[2]) and
                  (normrhsfluidvelocity_ / sqrtnpfv_ < toleranceresidualforceslist_[0]) and
                  (normrhsporofluidpressure_ / sqrtnpfp_ < toleranceresidualforceslist_[1]) and
                  (normrhsale_ / sqrtna_ < toleranceresidualforceslist_[5]) and
                  (normrhsfluidvelocity_ / sqrtnfv_ < toleranceresidualforceslist_[3]) and
                  (normrhsfluidpressure_ / sqrtnfp_ < toleranceresidualforceslist_[4]));
      break;
    case Inpar::FPSI::relativconvergencenorm_sys:
      FOUR_C_THROW(
          "Check for convergence of residual forces with type <relativ_sys> not implemented yet!");
      break;
    default:
      FOUR_C_THROW("Cannot check for convergence of residual forces for any reason :-P !");
      break;
  }

  // combine increments and forces
  bool converged = false;
  if (combinedconvergence_ == Inpar::FPSI::bop_and)
    converged = convinc and convfres;
  else if (combinedconvergence_ == Inpar::FPSI::bop_or)
    converged = convinc or convfres;
  else
    FOUR_C_THROW("Something went terribly wrong with binary operator!");

  // return things
  return converged;
}  // Converged()


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::build_convergence_norms()
{
  rhs_->Norm2(&normofrhs_);
  Teuchos::RCP<const Epetra_Vector> rhs_fluid;
  Teuchos::RCP<const Epetra_Vector> rhs_fluidvelocity;
  Teuchos::RCP<const Epetra_Vector> rhs_fluidpressure;
  Teuchos::RCP<const Epetra_Vector> rhs_porofluidvelocity;
  Teuchos::RCP<const Epetra_Vector> rhs_porofluidpressure;
  Teuchos::RCP<const Epetra_Vector> rhs_porointerface;
  Teuchos::RCP<const Epetra_Vector> rhs_fluidinterface;
  Teuchos::RCP<const Epetra_Vector> rhs_porofluid;
  Teuchos::RCP<const Epetra_Vector> rhs_porostruct;
  Teuchos::RCP<const Epetra_Vector> rhs_ale;


  rhs_porostruct = Extractor().ExtractVector(rhs_, structure_block_);
  rhs_porofluid = Extractor().ExtractVector(rhs_, porofluid_block_);
  rhs_porofluidvelocity = poro_field()->fluid_field()->ExtractVelocityPart(rhs_porofluid);
  rhs_porofluidpressure = poro_field()->fluid_field()->ExtractPressurePart(rhs_porofluid);
  rhs_porointerface =
      FPSICoupl()->poro_fluid_fpsi_vel_pres_extractor()->ExtractCondVector(rhs_porofluid);

  rhs_fluid = Extractor().ExtractVector(rhs_, fluid_block_);
  //  Teuchos::RCP<const Epetra_Vector> rhs_fullfluid = Teuchos::rcp(new
  //  Epetra_Vector(*fluid_field()->dof_row_map())); Teuchos::RCP<const Epetra_Vector> rhs_fsi =
  //  Teuchos::rcp(new
  //  Epetra_Vector(*fluid_field()->Interface()->Map(FLD::UTILS::MapExtractor::cond_fsi),true));
  //  rhs_fullfluid = Core::LinAlg::MergeVector(rhs_fluid,rhs_fsi,false);

  rhs_fluidvelocity = fluid_field()->ExtractVelocityPart(rhs_fluid);
  rhs_fluidpressure = fluid_field()->ExtractPressurePart(rhs_fluid);
  rhs_fluidinterface = FPSICoupl()->fluid_fpsi_vel_pres_extractor()->ExtractCondVector(rhs_fluid);

  rhs_ale = Extractor().ExtractVector(rhs_, ale_i_block_);  // Extractor().ExtractVector(rhs_, 2);

  rhs_porostruct->Norm2(&normrhsporostruct_);
  rhs_fluid->Norm2(&normrhsfluid_);
  rhs_fluidvelocity->Norm2(&normrhsfluidvelocity_);
  rhs_fluidpressure->Norm2(&normrhsfluidpressure_);
  rhs_porofluidvelocity->Norm2(&normrhsporofluidvelocity_);
  rhs_porofluidpressure->Norm2(&normrhsporofluidpressure_);
  rhs_porointerface->Norm2(&normrhsporointerface_);
  rhs_fluidinterface->Norm2(&normrhsfluidinterface_);
  rhs_ale->Norm2(&normrhsale_);

  // get length of the porostructural, porofluid, fluid and ale vector
  sqrtnfv_ = rhs_fluidvelocity->GlobalLength();  // correct length here
  sqrtnfp_ = rhs_fluidpressure->GlobalLength();
  sqrtnpfv_ = rhs_porofluidvelocity->GlobalLength();
  sqrtnpfp_ = rhs_porofluidpressure->GlobalLength();
  sqrtnps_ = rhs_porostruct->GlobalLength();
  sqrtna_ = rhs_ale->GlobalLength();
  sqrtnall_ = sqrtnfv_ + sqrtnfp_ + sqrtnpfv_ + sqrtnpfp_ + sqrtnps_ + sqrtna_;

  sqrtnfv_ = sqrt(sqrtnfv_);
  sqrtnfp_ = sqrt(sqrtnfp_);
  sqrtnpfv_ = sqrt(sqrtnpfv_);
  sqrtnpfp_ = sqrt(sqrtnpfp_);
  sqrtnps_ = sqrt(sqrtnps_);
  sqrtna_ = sqrt(sqrtna_);
  sqrtnall_ = sqrt(sqrtnall_);

  if (islinesearch_ == false) iterinc_->Norm2(&normofiterinc_);

  Teuchos::RCP<const Epetra_Vector> iterincporostruct;
  Teuchos::RCP<const Epetra_Vector> iterincporofluid;
  Teuchos::RCP<const Epetra_Vector> iterincfluid;
  Teuchos::RCP<const Epetra_Vector> iterincfluidvelocity;
  Teuchos::RCP<const Epetra_Vector> iterincfluidpressure;
  Teuchos::RCP<const Epetra_Vector> iterincporofluidvelocity;
  Teuchos::RCP<const Epetra_Vector> iterincporofluidpressure;
  Teuchos::RCP<const Epetra_Vector> iterincale;
  Teuchos::RCP<const Epetra_Vector> iterincporointerface;
  Teuchos::RCP<const Epetra_Vector> iterincfluidinterface;

  iterincporostruct = Extractor().ExtractVector(iterinc_, structure_block_);
  iterincporofluid = Extractor().ExtractVector(iterinc_, porofluid_block_);
  iterincporofluidvelocity = poro_field()->fluid_field()->ExtractVelocityPart(iterincporofluid);
  iterincporofluidpressure = poro_field()->fluid_field()->ExtractPressurePart(iterincporofluid);
  iterincporointerface =
      FPSICoupl()->poro_fluid_fpsi_vel_pres_extractor()->ExtractCondVector(iterincporofluid);

  iterincfluid = Extractor().ExtractVector(iterinc_, fluid_block_);

  //  Teuchos::RCP<const Epetra_Vector> iterincfullfluid = Teuchos::rcp(new
  //  Epetra_Vector(*fluid_field()->dof_row_map())); Teuchos::RCP<const Epetra_Vector> iterincfsi =
  //  Teuchos::rcp(new
  //  Epetra_Vector(*fluid_field()->Interface()->Map(FLD::UTILS::MapExtractor::cond_fsi),true));
  //  iterincfullfluid = Core::LinAlg::MergeVector(iterincfluid,iterincfsi,false);

  iterincfluidvelocity = fluid_field()->ExtractVelocityPart(iterincfluid);
  iterincfluidpressure = fluid_field()->ExtractPressurePart(iterincfluid);
  iterincfluidinterface =
      FPSICoupl()->fluid_fpsi_vel_pres_extractor()->ExtractCondVector(iterincfluid);

  iterincale = Extractor().ExtractVector(iterinc_, ale_i_block_);

  iterincporostruct->Norm2(&normofiterincporostruct_);
  iterincporofluid->Norm2(&normofiterincporofluid_);
  iterincporofluidvelocity->Norm2(&normofiterincporofluidvelocity_);
  iterincporofluidpressure->Norm2(&normofiterincporofluidpressure_);
  iterincporointerface->Norm2(&normofiterincporointerface_);

  iterincfluid->Norm2(&normofiterincfluid_);
  iterincfluidvelocity->Norm2(&normofiterincfluidvelocity_);
  iterincfluidpressure->Norm2(&normofiterincfluidpressure_);
  iterincfluidinterface->Norm2(&normofiterincfluidinterface_);

  iterincale->Norm2(&normofiterincale_);

  // Get Norm1 of dof values for each field

  Teuchos::RCP<const Epetra_Vector> porofluiddof;
  Teuchos::RCP<const Epetra_Vector> porofluidvelocity;
  Teuchos::RCP<const Epetra_Vector> porofluidpressure;
  Teuchos::RCP<const Epetra_Vector> porostructdisplacements;
  Teuchos::RCP<const Epetra_Vector> fluiddof;
  Teuchos::RCP<const Epetra_Vector> fluidvelocity;
  Teuchos::RCP<const Epetra_Vector> fluidpressure;
  Teuchos::RCP<const Epetra_Vector> aledisplacements;
  Teuchos::RCP<const Epetra_Vector> porointerface;
  Teuchos::RCP<const Epetra_Vector> fluidinterface;

  porofluiddof = poro_field()->fluid_field()->Velnp();
  porofluidvelocity = poro_field()->fluid_field()->ExtractVelocityPart(porofluiddof);
  porofluidpressure = poro_field()->fluid_field()->ExtractPressurePart(porofluiddof);
  porostructdisplacements = poro_field()->structure_field()->Dispnp();
  fluiddof = fluid_field()->Velnp();
  fluidvelocity = fluid_field()->ExtractVelocityPart(fluiddof);
  fluidpressure = fluid_field()->ExtractPressurePart(fluiddof);
  aledisplacements = ale_field()->Dispnp();


  porofluidvelocity->Norm1(&norm1_pfv_);
  porofluidpressure->Norm1(&norm1_pfp_);
  porostructdisplacements->Norm1(&norm1_ps_);
  fluidvelocity->Norm1(&norm1_fv_);
  fluidpressure->Norm1(&norm1_fp_);
  aledisplacements->Norm1(&norm1_a_);
  norm1_alldof_ = norm1_pfv_ + norm1_pfp_ + norm1_ps_ + norm1_fv_ + norm1_fp_ + norm1_a_;

  // add small number to avoid division by 0 --> division by 10^-10 results anyway in 'not
  // converged'
  norm1_alldof_ += 1e-10;
  norm1_pfv_ += 1e-10;
  norm1_pfp_ += 1e-10;
  norm1_ps_ += 1e-10;
  norm1_fv_ += 1e-10;
  norm1_fp_ += 1e-10;
  norm1_a_ += 1e-10;

  return;
}  // build_convergence_norms

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::print_newton_iter()
{
  // print to standard out
  // replace myrank_ here general by Comm().MyPID()
  if ((Comm().MyPID() == 0) and printscreen_ and (Step() % printscreen_ == 0) and printiter_)
  {
    if (iter_ == 1) print_newton_iter_header(stdout);
    print_newton_iter_text(stdout);
  }

}  // print_newton_iter()

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::print_newton_iter_header(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter converged state etc
  oss << std::setw(8) << "numiter";

  // different style due relative or absolute error checking
  // displacement
  switch (normtypefres_)
  {
    case Inpar::FPSI::absoluteconvergencenorm:
      oss << std::setw(11) << "abs-res";
      break;
    case Inpar::FPSI::absoluteconvergencenorm_sys_split:
      oss << std::setw(11) << "abs-res_s";
      break;
    case Inpar::FPSI::relativconvergencenorm_sys:
      oss << std::setw(11) << "rel-res_s";
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case Inpar::FPSI::absoluteconvergencenorm:
      oss << std::setw(11) << "abs-inc";
      break;
    case Inpar::FPSI::absoluteconvergencenorm_sys_split:
      oss << std::setw(11) << "abs-inc_s";
      break;
    case Inpar::FPSI::relativconvergencenorm_sys:
      oss << std::setw(11) << "rel-inc_s";
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }

  switch (normtypefres_)
  {
    case Inpar::FPSI::absoluteconvergencenorm:
    case Inpar::FPSI::absoluteconvergencenorm_sys_split:
    case Inpar::FPSI::relativconvergencenorm_sys:
      oss << std::setw(16) << "poro-s-res";
      // oss <<std::setw(15)<< "abs-f-res";
      oss << std::setw(15) << "poro-fvel-res";
      oss << std::setw(15) << "poro-fpres-res";
      oss << std::setw(15) << "fld-fvel-res";
      oss << std::setw(15) << "fld-fpres-res";
      // oss <<std::setw(15)<< "pinterface-res";
      // oss <<std::setw(15)<< "finterface-res";
      break;
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case Inpar::FPSI::absoluteconvergencenorm:
    case Inpar::FPSI::absoluteconvergencenorm_sys_split:
      oss << std::setw(15) << "poro-s-inc";
      // oss <<std::setw(15)<< "abs-f-inc";
      oss << std::setw(16) << "poro-fvel-inc";
      oss << std::setw(16) << "poro-fpres-inc";
      oss << std::setw(15) << "fld-fvel-inc";
      oss << std::setw(15) << "fld-fpres-inc";
      oss << std::setw(11) << "ale-inc";
      oss << std::setw(14) << "poro-int-inc";
      oss << std::setw(14) << "fld-int-inc";
      break;
    case Inpar::FPSI::relativconvergencenorm_sys:
      oss << std::setw(15) << "poro-s-inc";
      // oss <<std::setw(15)<< "abs-f-inc";
      oss << std::setw(16) << "poro-fvel-inc";
      oss << std::setw(16) << "poro-fpres-inc";
      oss << std::setw(15) << "fld-fvel-inc";
      oss << std::setw(15) << "fld-fpres-inc";
      oss << std::setw(11) << "ale-inc";
      break;
    default:
      FOUR_C_THROW("Begin at the beginning and go on till you come to the end. Then stop.");
      break;
  }

  //  // add solution time
  //  oss << std::setw(14)<< "wct";

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile == nullptr) FOUR_C_THROW("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;
}  // print_newton_iter_header()


/*----------------------------------------------------------------------*
 | print Newton-Raphson iteration to screen                       |
 | originally by lw 12/07, tk 01/08                                     |
 *----------------------------------------------------------------------*/
void FPSI::Monolithic::print_newton_iter_text(FILE* ofile)
{
  // open outstringstream
  std::ostringstream oss;

  // enter  state etc
  oss << std::setw(4) << iter_;

  // different style due relative or absolute error checking
  // displacement
  switch (normtypefres_)
  {
    case Inpar::FPSI::absoluteconvergencenorm:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normofrhs_;
      break;
    case Inpar::FPSI::absoluteconvergencenorm_sys_split:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normofrhs_ / sqrtnall_;
      break;
    case Inpar::FPSI::relativconvergencenorm_sys:
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case Inpar::FPSI::absoluteconvergencenorm:
      oss << std::setw(12) << std::setprecision(5) << std::scientific << normofiterinc_;
      break;
    case Inpar::FPSI::relativconvergencenorm_sys:
      oss << std::setw(12) << std::setprecision(5) << std::scientific
          << normofiterinc_ / norm1_alldof_ * sqrtnall_;
      break;
    case Inpar::FPSI::absoluteconvergencenorm_sys_split:
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }

  switch (normtypefres_)
  {
    case Inpar::FPSI::absoluteconvergencenorm:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporostruct_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporofluidvelocity_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporofluidpressure_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsfluidvelocity_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsfluidpressure_;
      // oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporointerface_;
      // oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsfluidinterface_;
      break;
    case Inpar::FPSI::absoluteconvergencenorm_sys_split:
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normrhsporostruct_ / sqrtnps_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normrhsporofluidvelocity_ / sqrtnpfv_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normrhsporofluidpressure_ / sqrtnpfp_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normrhsfluidvelocity_ / sqrtnfv_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normrhsfluidpressure_ / sqrtnfp_;
      // oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsporointerface_;
      // oss << std::setw(15) << std::setprecision(5) << std::scientific << normrhsfluidinterface_;
      break;
    case Inpar::FPSI::relativconvergencenorm_sys:
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }

  switch (normtypeinc_)
  {
    case Inpar::FPSI::absoluteconvergencenorm:
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normofiterincporostruct_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincporofluidvelocity_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincporofluidpressure_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincfluidvelocity_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincfluidpressure_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific << normofiterincale_;
      oss << std::setw(13) << std::setprecision(5) << std::scientific
          << normofiterincfluidinterface_;
      oss << std::setw(14) << std::setprecision(5) << std::scientific
          << normofiterincporointerface_;
      break;
    case Inpar::FPSI::relativconvergencenorm_sys:
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincporostruct_ / norm1_ps_ * sqrtnps_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincporofluidvelocity_ / norm1_pfv_ * sqrtnpfv_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincporofluidpressure_ / norm1_pfp_ * sqrtnpfp_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincfluidvelocity_ / norm1_fv_ * sqrtnfv_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincfluidpressure_ / norm1_fp_ * sqrtnfp_;
      oss << std::setw(15) << std::setprecision(5) << std::scientific
          << normofiterincale_ / norm1_a_ * sqrtna_;
      break;
    case Inpar::FPSI::absoluteconvergencenorm_sys_split:
    default:
      FOUR_C_THROW("You should not turn up here.");
      break;
  }

  // finish oss
  oss << std::ends;

  // print to screen (could be done differently...)
  if (ofile == nullptr) FOUR_C_THROW("no ofile available");
  fprintf(ofile, "%s\n", oss.str().c_str());

  // print it, now
  fflush(ofile);

  // nice to have met you
  return;

}  // PrintN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::SetupNewton()
{
  // initialise equilibrium loop and norms
  iter_ = 1;
  normofrhs_ = 0.0;
  normofiterinc_ = 0.0;
  normrhsfluid_ = 0.0;
  normofiterincfluid_ = 0.0;
  normrhsfluidvelocity_ = 0.0;
  normofiterincfluidvelocity_ = 0.0;
  normrhsporostruct_ = 0.0;
  normofiterincporostruct_ = 0.0;
  normofiterincporofluid_ = 0.0;
  normrhsfluidpressure_ = 0.0;
  normofiterincfluidpressure_ = 0.0;
  normrhsporofluidvelocity_ = 0.0;
  normofiterincporofluidvelocity_ = 0.0;
  normrhsporofluidpressure_ = 0.0;
  normofiterincporofluidpressure_ = 0.0;
  normrhsporointerface_ = 0.0;
  normofiterincporointerface_ = 0.0;
  normrhsfluidinterface_ = 0.0;
  normofiterincfluidinterface_ = 0.0;
  sqrtnall_ = 1;
  sqrtnfv_ = 1;
  sqrtnfp_ = 1;
  sqrtnpfv_ = 1;
  sqrtnpfp_ = 1;
  sqrtnps_ = 1;
  sqrtna_ = 1;
  norm1_alldof_ = 1.0;
  norm1_fv_ = 1.0;
  norm1_fp_ = 1.0;
  norm1_pfv_ = 1.0;
  norm1_pfp_ = 1.0;
  norm1_ps_ = 1.0;
  norm1_a_ = 1.0;


  // incremental solution vector with length of all dofs
  iterinc_ = Core::LinAlg::CreateVector(*dof_row_map(), true);
  iterinc_->PutScalar(0.0);

  // a zero vector of full length
  zeros_ = Core::LinAlg::CreateVector(*dof_row_map(), true);
  zeros_->PutScalar(0.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::FPSIFDCheck()
{
  // FD check is nice to check your linearisations, but be aware that we did not linearize following
  // terms:
  // gridvelocity in convective term, dispnp for the stabilization
  active_FD_check_ = true;
  // Set states
  poro_field()->fluid_field()->discretization()->ClearState();

  poro_field()->fluid_field()->discretization()->set_state(
      0, "dispnp", poro_field()->fluid_field()->Dispnp());
  poro_field()->fluid_field()->discretization()->set_state(
      0, "gridv", poro_field()->fluid_field()->GridVel());
  poro_field()->fluid_field()->discretization()->set_state(
      0, "dispn", poro_field()->fluid_field()->Dispn());
  poro_field()->fluid_field()->discretization()->set_state(
      0, "veln", poro_field()->fluid_field()->Veln());
  poro_field()->fluid_field()->discretization()->set_state(
      0, "velaf", poro_field()->fluid_field()->Velnp());
  poro_field()->fluid_field()->discretization()->set_state(
      0, "velnp", poro_field()->fluid_field()->Velnp());

  fluid_field()->discretization()->ClearState();

  fluid_field()->discretization()->set_state(0, "dispnp", fluid_field()->Dispnp());
  fluid_field()->discretization()->set_state(0, "gridv", fluid_field()->GridVel());
  fluid_field()->discretization()->set_state(0, "dispn", fluid_field()->Dispn());
  fluid_field()->discretization()->set_state(0, "veln", fluid_field()->Veln());
  fluid_field()->discretization()->set_state(0, "velaf", fluid_field()->Velnp());
  fluid_field()->discretization()->set_state(0, "velnp", fluid_field()->Velnp());

  // define and set toggle parameter delta
  const double delta = 1e-8;

  Global::Problem* problem = Global::Problem::Instance();
  const Teuchos::ParameterList& fpsidynparams = problem->FPSIDynamicParams();
  int columntocheck = fpsidynparams.get<int>("FDCheck_column");
  int rowtocheck = fpsidynparams.get<int>("FDCheck_row");
  //////////////////////////////////////////////////////////////7
  // matrices and vectors
  //////////////////////////////////////////////////////////////
  // build artificial iteration increment
  Teuchos::RCP<Epetra_Vector> iterinc = Core::LinAlg::CreateVector(*dof_row_map(), true);
  const int dofs = iterinc->GlobalLength();
  iterinc->PutScalar(0.0);
  iterinc->ReplaceGlobalValue(0, 0, delta);

  // build approximated FD stiffness matrix
  Teuchos::RCP<Epetra_CrsMatrix> stiff_approx = Core::LinAlg::CreateMatrix(*dof_row_map(), 81);

  // store old rhs
  Teuchos::RCP<Epetra_Vector> rhs_old = Teuchos::rcp(new Epetra_Vector(*dof_row_map(), true));
  rhs_old->Update(1.0, *rhs_, 0.0);

  Teuchos::RCP<Epetra_Vector> rhs_copy = Teuchos::rcp(new Epetra_Vector(*dof_row_map(), true));

  Teuchos::RCP<Core::LinAlg::SparseMatrix> sparse = systemmatrix_->Merge();

  Teuchos::RCP<Core::LinAlg::SparseMatrix> sparse_copy =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(sparse->EpetraMatrix(), Core::LinAlg::Copy));


  std::cout << "\n****************** FPSI finite difference check ******************" << std::endl;

  int dof_poro_struct = (poro_field()->structure_field()->discretization()->NumGlobalNodes()) * 3;
  int dof_poro_fluid = (poro_field()->fluid_field()->discretization()->NumGlobalNodes()) * 4;
  int dof_fluid = (fluid_field()->discretization()->NumGlobalNodes()) * 4;
  int dof_ale = (ale_field()->discretization()->NumGlobalNodes()) * 3;

  std::cout << "poro structure field has " << dof_poro_struct << " DOFs" << std::endl;
  std::cout << "poro fluid field has     " << dof_poro_fluid << " DOFs" << std::endl;
  std::cout << "fluid field has          " << dof_fluid << " DOFs" << std::endl;
  std::cout << "ale field has            " << dof_ale << " DOFs" << std::endl;
  std::cout << "in total                 " << dofs << " DOFs" << std::endl;


  for (int i_loc = 0; i_loc < dofs; ++i_loc)  // loop over columns
  {
    int i = dof_row_map()->GID(i_loc);
    int im1 = dof_row_map()->GID(i_loc - 1);
    int ip1 = dof_row_map()->GID(i_loc + 1);
    std::cout << i << "... " << std::flush;
    if (combined_dbc_map()->MyGID(i) or fluid_field()->Interface()->FSICondMap()->MyGID(i))
    {
      iterinc->ReplaceGlobalValue(i, 0, 0.0);
    }
    evaluate(iterinc);  // initial iterinc is varied at first dof (0-th element)
    setup_system_matrix();

    rhs_copy->Update(1.0, *rhs_, 0.0);

    iterinc_->PutScalar(0.0);  // Useful? depends on solver and more
    poro_field()->ClearPoroIterinc();
    Core::LinAlg::apply_dirichlet_to_system(
        *sparse_copy, *iterinc_, *rhs_copy, *zeros_, *fluid_field()->Interface()->FSICondMap());
    Core::LinAlg::apply_dirichlet_to_system(
        *sparse_copy, *iterinc_, *rhs_copy, *zeros_, *combined_dbc_map());

    rhs_copy->Update(-1.0, *rhs_old, 1.0);  // finite difference approximation of partial derivative
    rhs_copy->Scale(-1.0 / delta);

    if (i == columntocheck)
    {
      std::cout << "iterinc:  " << *iterinc << std::endl;
      std::cout << "rhs_old:  " << *rhs_old << std::endl;
      std::cout << "rhs_copy: " << *rhs_copy << std::endl;
      FOUR_C_THROW("Stopped by FPSI - fd_check!");
    }

    int* index = &i;
    for (int j_loc = 0; j_loc < dofs; ++j_loc)  // loop over rows
    {
      int j = dof_row_map()->GID(j_loc);
      double value = (*rhs_copy)[j_loc];
      stiff_approx->InsertGlobalValues(
          j, 1, &value, index);  // int InsertGlobalValues(int GlobalRow, int NumEntries, double*
                                 // Values, int* Indices);

    }  // j-loop (rows)

    if (not combined_dbc_map()->MyGID(i) and not fluid_field()->Interface()->FSICondMap()->MyGID(i))
      iterinc->ReplaceGlobalValue(i, 0, -delta);

    iterinc->ReplaceGlobalValue(im1, 0, 0.0);

    if (i_loc != dofs - 1) iterinc->ReplaceGlobalValue(ip1, 0, delta);

  }  // i-loop (columns)
  evaluate(iterinc);
  setup_system_matrix();

  int err = stiff_approx->FillComplete();
  if (err) FOUR_C_THROW("FD_Check: FillComplete failed with err-code: %d", err);

  Teuchos::RCP<Core::LinAlg::SparseMatrix> temp =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(stiff_approx, Core::LinAlg::Copy));

  Teuchos::RCP<Epetra_CrsMatrix> stiff_approx_sparse = temp->EpetraMatrix();

  Teuchos::RCP<Epetra_CrsMatrix> sparse_crs = sparse_copy->EpetraMatrix();

  // calc error (subtraction of sparse_crs and stiff_approx_sparse)
  for (int i_loc = 0; i_loc < dofs; i_loc++)
  {
    int i = dof_row_map()->GID(i_loc);
    int length;
    int numentries = sparse_crs->NumGlobalEntries(i);
    std::vector<double> values(numentries);
    std::vector<int> indices(numentries);
    sparse_crs->ExtractGlobalRowCopy(i, numentries, length, values.data(), indices.data());

    for (int k = 0; k < numentries; k++)
    {
      values[k] = -values[k];
    }

    stiff_approx_sparse->SumIntoGlobalValues(i, numentries, values.data(), indices.data());
  }
  stiff_approx_sparse->FillComplete();
  sparse_crs->FillComplete();

  bool success = true;
  double error_max = 0.0;
  for (int i_loc = 0; i_loc < dofs; ++i_loc)
  {
    int i = dof_row_map()->GID(i_loc);
    if (not combined_dbc_map()->MyGID(i) and
        not fluid_field()->Interface()->FSICondMap()->MyGID(
            i))  // only if there is no dirichlet value on dof and dof is not condensed
    {
      for (int j_loc = 0; j_loc < dofs; ++j_loc)
      {
        int j = dof_row_map()->GID(j_loc);
        if (not combined_dbc_map()->MyGID(j) and
            not fluid_field()->Interface()->FSICondMap()->MyGID(j))
        {
          double stiff_approx_ij = 0.0;
          double sparse_ij = 0.0;
          double error_ij = 0.0;

          // get error_crs entry ij
          int errornumentries;
          int errorlength = stiff_approx_sparse->NumGlobalEntries(i);
          std::vector<double> errorvalues(errorlength);
          std::vector<int> errorindices(errorlength);
          stiff_approx_sparse->ExtractGlobalRowCopy(
              i, errorlength, errornumentries, errorvalues.data(), errorindices.data());
          for (int k = 0; k < errorlength; ++k)
          {
            if (errorindices[k] == j)
            {
              error_ij = errorvalues[k];
              break;
            }
            else
              error_ij = 0.0;
          }

          // get sparse_ij entry ij
          int sparsenumentries;
          int sparselength = sparse_crs->NumGlobalEntries(i);
          std::vector<double> sparsevalues(sparselength);
          std::vector<int> sparseindices(sparselength);
          sparse_crs->ExtractGlobalRowCopy(
              i, sparselength, sparsenumentries, sparsevalues.data(), sparseindices.data());
          for (int k = 0; k < sparselength; ++k)
          {
            if (sparseindices[k] == j)
            {
              sparse_ij = sparsevalues[k];
              break;
            }
            else
              sparse_ij = 0.0;
          }
          // get stiff_approx entry ijs
          int approxnumentries;
          int approxlength = stiff_approx->NumGlobalEntries(i);
          std::vector<double> approxvalues(approxlength);
          std::vector<int> approxindices(approxlength);
          stiff_approx->ExtractGlobalRowCopy(
              i, approxlength, approxnumentries, approxvalues.data(), approxindices.data());
          for (int k = 0; k < approxlength; ++k)
          {
            if (approxindices[k] == j)
            {
              stiff_approx_ij = approxvalues[k];
              break;
            }
            else
              stiff_approx_ij = 0.0;
          }

          double error = 0.0;
          if (abs(stiff_approx_ij) > 1e-3)
            error = error_ij / stiff_approx_ij;
          else if (abs(sparse_ij) > 1e-3)
            error = error_ij / sparse_ij;

          if (i == rowtocheck and j == columntocheck)
          {
            std::cout << "K_approx value: " << stiff_approx_ij << std::endl;
            std::cout << "K value : " << sparse_ij << std::endl;
            std::cout << "error : " << error << std::endl;
            std::cout << "error_ij : " << error_ij << std::endl;
            std::cout << "i : " << i << std::endl;
            std::cout << "j : " << j << std::endl;
          }

          if (abs(error) > abs(error_max)) error_max = abs(error);

          if ((abs(error) > 1e-4))
          {
            if (abs(error_ij) > 1e-4)
            {
              std::cout << "finite difference check failed entry (" << i << "," << j
                        << ")! stiff: " << sparse_ij << ", approx: " << stiff_approx_ij
                        << " ,abs. error: " << error_ij << " , rel. error: " << error << std::endl;

              success = false;
            }
          }
        }
      }
    }
  }  // loop over dofs of succes check

  if (success)
  {
    std::cout << "finite difference check successful, max. rel. error: " << error_max << std::endl;
    std::cout << "****************** finite difference check done ***************\n\n" << std::endl;
  }
  else
    std::cout << "FPSIFDCheck failed" << std::endl;
  // FOUR_C_THROW("FPSIFDCheck failed");
  poro_field()->fluid_field()->discretization()->ClearState();
  fluid_field()->discretization()->ClearState();

  active_FD_check_ = false;
  return;
}
/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::extract_columnsfrom_sparse(Teuchos::RCP<Epetra_CrsMatrix> src,
    const Teuchos::RCP<const Epetra_Map>& colmap, Teuchos::RCP<Epetra_CrsMatrix> dst)
{
  dst->PutScalar(0.0);  // clear matrix
  int rows = src->NumGlobalRows();
  for (int row = 0; row < rows; ++row)
  {
    int g_row = src->RangeMap().GID(row);
    int numentries;
    int length = src->NumGlobalEntries(g_row);
    std::vector<double> values(length);
    std::vector<int> indices(length);
    src->ExtractGlobalRowCopy(g_row, length, numentries, values.data(), indices.data());
    for (int col = 0; col < length; ++col)  // loop over non-zero columns in active row
    {
      if (colmap->LID(indices[col]) != -1)
      {
        dst->InsertGlobalValues(
            g_row, 1, &values[col], &indices[col]);  // add column value of active row!
      }
    }
  }
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FPSI::Monolithic::SetConductivity(double conduct)
{
  if (FPSICoupl() != Teuchos::null) FPSICoupl()->SetConductivity(conduct);
  conductivity_ = conduct;  // remove me...
}

FOUR_C_NAMESPACE_CLOSE
