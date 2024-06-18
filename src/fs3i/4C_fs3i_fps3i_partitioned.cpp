/*----------------------------------------------------------------------*/
/*! \file
\brief General algorithmic routines for partitioned solution approaches
       to fluid-porous-structure-scalar-scalar interaction (FPS3I), that is,
       algorithmic routines not specifically related to partitioned
       solution approaches to one -or two-way-coupled problem
       configurations, respectively.

\level 3



*----------------------------------------------------------------------*/

#include "4C_fs3i_fps3i_partitioned.hpp"

#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_fem_condition_selector.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_fluid_ele_action.hpp"
#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_fpsi_monolithic_plain.hpp"
#include "4C_fpsi_utils.hpp"
#include "4C_fsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linear_solver_method.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_poroelast_monolithic.hpp"
#include "4C_poroelast_scatra_utils.hpp"
#include "4C_poroelast_scatra_utils_clonestrategy.hpp"
#include "4C_scatra_algorithm.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_utils_clonestrategy.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  Constructor                                           hemmler 07/14 |
 *----------------------------------------------------------------------*/
FS3I::PartFPS3I::PartFPS3I(const Epetra_Comm& comm) : FS3IBase(), comm_(comm)
{
  // keep empty
  return;
}


/*----------------------------------------------------------------------*
 |  Init                                                    rauch 09/16 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::init()
{
  FS3I::FS3IBase::init();

  if (comm_.MyPID() == 0)
  {
    // ##################       0.- Warning          //#########################
    std::cout << std::endl;
    std::cout << "##############################################################################"
              << std::endl;
    std::cout << "################################# WARNING!!! #################################"
              << std::endl;
    std::cout << "##############################################################################"
              << std::endl;
    std::cout << std::endl;
    std::cout << "This version of Fluid-porous-structure-scatra-scatra interaction (FPS3I) does NOT"
              << std::endl;
    std::cout << "account for the convective scalar transport at the fluid-poro interface!"
              << std::endl;
    std::cout << "The conservation of mass at the interface is only guaranteed for purely "
                 "diffusive transport"
              << std::endl;
    std::cout << std::endl;
    std::cout << "##############################################################################"
              << std::endl;
    std::cout << "################################# WARNING!!! #################################"
              << std::endl;
    std::cout << "##############################################################################"
              << std::endl;
    std::cout << std::endl;
  }
  // ##################       1.- Parameter reading          //#########################
  Global::Problem* problem = Global::Problem::Instance();
  const Teuchos::ParameterList& fs3idyn = problem->FS3IDynamicParams();
  const Teuchos::ParameterList& fpsidynparams = problem->FPSIDynamicParams();
  const Teuchos::ParameterList& poroelastdynparams = problem->poroelast_dynamic_params();
  const Teuchos::ParameterList& scatradyn = problem->scalar_transport_dynamic_params();

  double dt_fpsi = fpsidynparams.get<double>("TIMESTEP");
  double dt_poroelast = poroelastdynparams.get<double>("TIMESTEP");
  if (dt_fpsi != dt_poroelast)
  {
    FOUR_C_THROW(
        "Please set \"TIMESTEP\" in \"POROELASTICITY DYNAMIC\" to the same value as in \"FPSI "
        "DYNAMIC\"!");
  }

  Teuchos::RCP<FPSI::Utils> FPSI_UTILS = FPSI::Utils::Instance();

  // ##################    2.- Creation of Poroelastic + Fluid problem. (discretization called
  //  inside)     //##################
  Teuchos::RCP<FPSI::FpsiBase> fpsi_algo = Teuchos::null;

  fpsi_algo = FPSI_UTILS->setup_discretizations(comm_, fpsidynparams, poroelastdynparams);

  // only monolithic coupling of fpsi problem is supported!
  int coupling = Core::UTILS::IntegralValue<int>(fpsidynparams, "COUPALGO");
  if (coupling == fpsi_monolithic_plain)
  {
    // Cast needed because functions such as poro_field() and fluid_field() are just a
    // member-functions of the derived class MonolithicPlain, but not of the base class FPSI_Base
    fpsi_ = Teuchos::rcp_dynamic_cast<FPSI::MonolithicPlain>(fpsi_algo);
  }
  else
  {
    FOUR_C_THROW(
        "Partitioned solution scheme not implemented for FPSI, yet. "
        "Make sure that the parameter COUPALGO is set to 'fpsi_monolithic_plain', "
        "and the parameter PARITIONED is set to 'monolithic'. ");
  }

  // ##################      3. discretization of Scatra problem       //##################
  problem->GetDis("scatra1")->fill_complete();
  problem->GetDis("scatra2")->fill_complete();

  //---------------------------------------------------------------------
  // access discretizations for poro (structure) and fluid as well as fluid-
  // and poro-based scalar transport and get material map for fluid
  // and scalar transport elements
  //---------------------------------------------------------------------
  Teuchos::RCP<Core::FE::Discretization> fluiddis = problem->GetDis("fluid");
  Teuchos::RCP<Core::FE::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<Core::FE::Discretization> fluidscatradis = problem->GetDis("scatra1");
  Teuchos::RCP<Core::FE::Discretization> structscatradis = problem->GetDis("scatra2");

  // determine type of scalar transport
  const Inpar::ScaTra::ImplType impltype_fluid =
      Core::UTILS::IntegralValue<Inpar::ScaTra::ImplType>(
          Global::Problem::Instance()->FS3IDynamicParams(), "FLUIDSCAL_SCATRATYPE");

  //---------------------------------------------------------------------
  // create discretization for fluid-based scalar transport from and
  // according to fluid discretization
  //---------------------------------------------------------------------
  if (fluiddis->NumGlobalNodes() == 0) FOUR_C_THROW("Fluid discretization is empty!");


  // std::map<std::pair<std::string,std::string>,std::map<int,int> > clonefieldmatmap =
  // problem->CloningMaterialMap(); if (clonefieldmatmap.size() < 2)
  //  FOUR_C_THROW("At least two material lists required for partitioned FS3I!");

  // create fluid-based scalar transport elements if fluid-based scalar
  // transport discretization is empty
  if (fluidscatradis->NumGlobalNodes() == 0)
  {
    // fill fluid-based scatra discretization by cloning fluid discretization
    Core::FE::CloneDiscretization<ScaTra::ScatraFluidCloneStrategy>(
        fluiddis, fluidscatradis, Global::Problem::Instance()->CloningMaterialMap());
    fluidscatradis->fill_complete();

    // set implementation type of cloned scatra elements to advanced reactions
    for (int i = 0; i < fluidscatradis->NumMyColElements(); ++i)
    {
      Discret::ELEMENTS::Transport* element =
          dynamic_cast<Discret::ELEMENTS::Transport*>(fluidscatradis->lColElement(i));
      if (element == nullptr)
        FOUR_C_THROW("Invalid element type!");
      else
        element->SetImplType(impltype_fluid);
    }
  }
  else
    FOUR_C_THROW("Fluid AND ScaTra discretization present. This is not supported.");

  //---------------------------------------------------------------------
  // create discretization for poro-based scalar transport from and
  // according to poro (structure) discretization
  //--------------------------------------------------------------------
  if (fluiddis->NumGlobalNodes() == 0) FOUR_C_THROW("Fluid discretization is empty!");

  if (!structscatradis->Filled()) structscatradis->fill_complete();
  if (structscatradis->NumGlobalNodes() == 0)
  {
    // fill poro-based scatra discretization by cloning structure discretization
    Core::FE::CloneDiscretization<PoroElastScaTra::UTILS::PoroScatraCloneStrategy>(
        structdis, structscatradis, Global::Problem::Instance()->CloningMaterialMap());

    // redistribute FPSI interface here, since if done before the PoroScatra cloning does not work
    // fpsi_->redistribute_interface();
    // after redistributing the interface we have to fix the material pointers of the
    // structure-scatra discretisation
    // PoroElast::UTILS::SetMaterialPointersMatchingGrid(structdis,structscatradis);
  }
  else
    FOUR_C_THROW("Structure AND ScaTra discretization present. This is not supported.");

  // ##################      End of discretization       //##################

  //---------------------------------------------------------------------
  // create instances for fluid- and poro (structure)-based scalar transport
  // solver and arrange them in combined vector
  //---------------------------------------------------------------------
  // get the solver number used for structural ScalarTransport solver
  const int linsolver1number = fs3idyn.get<int>("LINEAR_SOLVER1");
  // get the solver number used for structural ScalarTransport solver
  const int linsolver2number = fs3idyn.get<int>("LINEAR_SOLVER2");

  // check if the linear solver has a valid solver number
  if (linsolver1number == (-1))
    FOUR_C_THROW(
        "no linear solver defined for fluid ScalarTransport solver. Please set LINEAR_SOLVER2 in "
        "FS3I DYNAMIC to a valid number!");
  if (linsolver2number == (-1))
    FOUR_C_THROW(
        "no linear solver defined for structural ScalarTransport solver. Please set LINEAR_SOLVER2 "
        "in FS3I DYNAMIC to a valid number!");
  fluidscatra_ = Teuchos::rcp(new Adapter::ScaTraBaseAlgorithm(
      fs3idyn, scatradyn, problem->SolverParams(linsolver1number), "scatra1", true));

  // now we can call init() on the scatra time integrator
  fluidscatra_->init();
  fluidscatra_->ScaTraField()->set_number_of_dof_set_displacement(1);
  fluidscatra_->ScaTraField()->set_number_of_dof_set_velocity(1);
  fluidscatra_->ScaTraField()->set_number_of_dof_set_wall_shear_stress(1);
  fluidscatra_->ScaTraField()->set_number_of_dof_set_pressure(1);

  structscatra_ = Teuchos::rcp(new Adapter::ScaTraBaseAlgorithm(
      fs3idyn, scatradyn, problem->SolverParams(linsolver2number), "scatra2", true));

  // only now we must call init() on the scatra time integrator.
  // all objects relying on the parallel distribution are
  // created and pointers are set.
  structscatra_->init();
  structscatra_->ScaTraField()->set_number_of_dof_set_displacement(1);
  structscatra_->ScaTraField()->set_number_of_dof_set_velocity(1);
  structscatra_->ScaTraField()->set_number_of_dof_set_wall_shear_stress(2);
  structscatra_->ScaTraField()->set_number_of_dof_set_pressure(2);

  scatravec_.push_back(fluidscatra_);
  scatravec_.push_back(structscatra_);

  //---------------------------------------------------------------------
  // check various input parameters
  //---------------------------------------------------------------------
  const Teuchos::ParameterList& structdyn = problem->structural_dynamic_params();
  const Teuchos::ParameterList& fluiddyn = problem->FluidDynamicParams();
  // check consistency of time-integration schemes in input file
  // (including parameter theta itself in case of one-step-theta scheme)
  // and rule out unsupported versions of generalized-alpha time-integration
  // scheme (as well as other inappropriate schemes) for fluid subproblem
  Inpar::ScaTra::TimeIntegrationScheme scatratimealgo =
      Core::UTILS::IntegralValue<Inpar::ScaTra::TimeIntegrationScheme>(scatradyn, "TIMEINTEGR");
  Inpar::FLUID::TimeIntegrationScheme fluidtimealgo =
      Core::UTILS::IntegralValue<Inpar::FLUID::TimeIntegrationScheme>(fluiddyn, "TIMEINTEGR");

  Inpar::STR::DynamicType structtimealgo =
      Core::UTILS::IntegralValue<Inpar::STR::DynamicType>(structdyn, "DYNAMICTYP");

  if (fluidtimealgo == Inpar::FLUID::timeint_one_step_theta)
  {
    if (scatratimealgo != Inpar::ScaTra::timeint_one_step_theta or
        structtimealgo != Inpar::STR::dyna_onesteptheta)
      FOUR_C_THROW(
          "Partitioned FS3I computations should feature consistent time-integration schemes for "
          "the subproblems; in this case, a one-step-theta scheme is intended to be used for the "
          "fluid subproblem, and different schemes are intended to be used for the structure "
          "and/or scalar transport subproblems!");

    if (scatradyn.get<double>("THETA") != fluiddyn.get<double>("THETA") or
        scatradyn.get<double>("THETA") != structdyn.sublist("ONESTEPTHETA").get<double>("THETA"))
      FOUR_C_THROW(
          "Parameter(s) theta for one-step-theta time-integration scheme defined in one or more of "
          "the individual fields do(es) not match for partitioned FS3I computation.");
  }
  else if (fluidtimealgo == Inpar::FLUID::timeint_afgenalpha)
  {
    if (scatratimealgo != Inpar::ScaTra::timeint_gen_alpha or
        structtimealgo != Inpar::STR::dyna_genalpha)
      FOUR_C_THROW(
          "Partitioned FS3I computations should feature consistent time-integration schemes for "
          "the subproblems; in this case, a (alpha_f-based) generalized-alpha scheme is intended "
          "to be used for the fluid subproblem, and different schemes are intended to be used for "
          "the structure and/or scalar transport subproblems!");
  }
  else if (fluidtimealgo == Inpar::FLUID::timeint_npgenalpha)
  {
    FOUR_C_THROW(
        "Partitioned FS3I computations do not support n+1-based generalized-alpha time-integration "
        "schemes for the fluid subproblem!");
  }
  else if (fluidtimealgo == Inpar::FLUID::timeint_bdf2 or
           fluidtimealgo == Inpar::FLUID::timeint_stationary)
  {
    FOUR_C_THROW(
        "Partitioned FS3I computations do not support stationary of BDF2 time-integration schemes "
        "for the fluid subproblem!");
  }

  // check that incremental formulation is used for scalar transport field,
  // according to structure and fluid field
  if (scatravec_[0]->ScaTraField()->IsIncremental() == false)
    FOUR_C_THROW("Incremental formulation required for partitioned FS3I computations!");


  return;
}


/*----------------------------------------------------------------------*
 |  Setup                                                   rauch 09/16 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::setup()
{
  FS3I::FS3IBase::setup();

  // only now we must call setup() on the scatra base algo.
  // all objects relying on the parallel distribution are
  // created and pointers are set.
  // calls setup() on time integrator inside.
  fluidscatra_->setup();
  structscatra_->setup();

  //---------------------------------------------------------------------
  // check existence of scatra coupling conditions for both
  // discretizations and definition of the permeability coefficient
  //---------------------------------------------------------------------
  CheckFS3IInputs();

  // in case of FPS3I we have to handle the conductivity, too
  Teuchos::RCP<Core::FE::Discretization> dis = scatravec_[0]->ScaTraField()->discretization();
  std::vector<Core::Conditions::Condition*> coupcond;
  dis->GetCondition("ScaTraCoupling", coupcond);
  double myconduct = coupcond[0]->parameters().get<double>(
      "hydraulic conductivity");  // here we assume the conductivity to be the same in every BC

  // conductivity is not only needed in scatracoupling but also in FPSI coupling!
  if (myconduct == 0.0)
  {
    FOUR_C_THROW(
        "conductivity of 0.0 is not allowed!!! Should be set in \"DESIGN SCATRA COUPLING SURF "
        "CONDITIONS\"");
  }
  fpsi_->SetConductivity(myconduct);

  return;
}

/*----------------------------------------------------------------------*
 |  Restart                                               hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::read_restart()
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  Global::Problem* problem = Global::Problem::Instance();
  const int restart = problem->restart();

  if (restart)
  {
    // restart of FPSI problem
    fpsi_->read_restart(restart);

    // restart of scatra problem
    for (unsigned i = 0; i < scatravec_.size(); ++i)
    {
      Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
      currscatra->ScaTraField()->read_restart(restart);
    }

    time_ = fpsi_->fluid_field()->Time();
    step_ = fpsi_->fluid_field()->Step();
  }
}

/*----------------------------------------------------------------------*
 | redistribute the FPSI interface                           thon 11/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::redistribute_interface()
{
  fpsi_->redistribute_interface();

  Global::Problem* problem = Global::Problem::Instance();

  if (comm_.NumProc() >
      1)  // if we have more than one processor, we need to redistribute at the FPSI interface
  {
    Teuchos::RCP<FPSI::Utils> FPSI_UTILS = FPSI::Utils::Instance();

    Teuchos::RCP<std::map<int, int>> Fluid_PoroFluid_InterfaceMap =
        FPSI_UTILS->get_fluid_poro_fluid_interface_map();
    Teuchos::RCP<std::map<int, int>> PoroFluid_Fluid_InterfaceMap =
        FPSI_UTILS->get_poro_fluid_fluid_interface_map();

    FPSI_UTILS->redistribute_interface(
        problem->GetDis("scatra1"), Teuchos::null, "", *PoroFluid_Fluid_InterfaceMap);
    FPSI_UTILS->redistribute_interface(
        problem->GetDis("scatra2"), Teuchos::null, "", *Fluid_PoroFluid_InterfaceMap);
  }

  Teuchos::RCP<Core::FE::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<Core::FE::Discretization> structscatradis = problem->GetDis("scatra2");

  // after redistributing the interface we have to fix the material pointers of the structure-scatra
  // discretisation
  PoroElast::UTILS::SetMaterialPointersMatchingGrid(structdis, structscatradis);
}

/*----------------------------------------------------------------------*
 |  System Setup                                          hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::SetupSystem()
{
  // do the coupling setup and create the combined dofmap

  // Setup FPSI system
  fpsi_->SetupSystem();
  fpsi_->SetupSolver();

  /*----------------------------------------------------------------------*/
  /*                  General set up for scalar fields                    */
  /*----------------------------------------------------------------------*/

  // create map extractors needed for scatra condition coupling
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> currscatra = scatravec_[i];
    const int numscal = currscatra->ScaTraField()->NumScal();
    Teuchos::RCP<Core::FE::Discretization> currdis = currscatra->ScaTraField()->discretization();
    Teuchos::RCP<Core::LinAlg::MultiMapExtractor> mapex =
        Teuchos::rcp(new Core::LinAlg::MultiMapExtractor());
    Core::Conditions::MultiConditionSelector mcs;
    mcs.AddSelector(Teuchos::rcp(
        new Core::Conditions::NDimConditionSelector(*currdis, "ScaTraCoupling", 0, numscal)));
    mcs.SetupExtractor(*currdis, *currdis->dof_row_map(), *mapex);
    scatrafieldexvec_.push_back(mapex);
  }

  scatracoup_->setup_condition_coupling(*(scatravec_[0]->ScaTraField()->discretization()),
      scatrafieldexvec_[0]->Map(1), *(scatravec_[1]->ScaTraField()->discretization()),
      scatrafieldexvec_[1]->Map(1), "ScaTraCoupling",
      scatravec_[0]
          ->ScaTraField()
          ->NumScal());  // we assume here that both discretisation have the same number of scalars

  // create map extractor for coupled scatra fields
  // the second field is always split
  std::vector<Teuchos::RCP<const Epetra_Map>> maps;

  // In the limiting case of an infinite permeability of the interface between
  // different scatra fields, the concentrations on both sides of the interface are
  // constrained to be equal. In this case, we keep the fluid scatra dofs at the
  // interface as unknowns in the overall system, whereas the poro (structure) scatra
  // dofs are condensed (cf. "structuresplit" in a monolithic FPSI
  // system). Otherwise, both concentrations are kept in the overall system
  // and the equality of fluxes is considered explicitly.
  if (infperm_)
  {
    maps.push_back(scatrafieldexvec_[0]->FullMap());
    maps.push_back(scatrafieldexvec_[1]->Map(0));
  }
  else
  {
    maps.push_back(scatrafieldexvec_[0]->FullMap());
    maps.push_back(scatrafieldexvec_[1]->FullMap());
  }
  Teuchos::RCP<Epetra_Map> fullmap = Core::LinAlg::MultiMapExtractor::MergeMaps(maps);
  scatraglobalex_->setup(*fullmap, maps);

  // create coupling vectors and matrices (only needed for finite surface permeabilities)
  if (not infperm_)
  {
    for (unsigned i = 0; i < scatravec_.size(); ++i)
    {
      Teuchos::RCP<Epetra_Vector> scatracoupforce =
          Teuchos::rcp(new Epetra_Vector(*(scatraglobalex_->Map(i)), true));
      scatracoupforce_.push_back(scatracoupforce);

      Teuchos::RCP<Core::LinAlg::SparseMatrix> scatracoupmat =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(scatraglobalex_->Map(i)), 27, false, true));
      scatracoupmat_.push_back(scatracoupmat);

      const Epetra_Map* dofrowmap = scatravec_[i]->ScaTraField()->discretization()->dof_row_map();
      Teuchos::RCP<Epetra_Vector> zeros = Core::LinAlg::CreateVector(*dofrowmap, true);
      scatrazeros_.push_back(zeros);
    }
  }
  // create scatra block matrix
  scatrasystemmatrix_ =
      Teuchos::rcp(new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
          *scatraglobalex_, *scatraglobalex_, 27, false, true));
  // create scatra rhs vector
  scatrarhs_ = Teuchos::rcp(new Epetra_Vector(*scatraglobalex_->FullMap(), true));
  // create scatra increment vector
  scatraincrement_ = Teuchos::rcp(new Epetra_Vector(*scatraglobalex_->FullMap(), true));
  // check whether potential Dirichlet conditions at scatra interface are
  // defined for both discretizations
  check_interface_dirichlet_bc();
  // scatra solver
  Teuchos::RCP<Core::FE::Discretization> firstscatradis =
      (scatravec_[0])->ScaTraField()->discretization();
#ifdef SCATRABLOCKMATRIXMERGE
  Teuchos::RCP<Teuchos::ParameterList> scatrasolvparams = Teuchos::rcp(new Teuchos::ParameterList);
  Core::UTILS::AddEnumClassToParameterList<Core::LinearSolver::SolverType>(
      "SOLVER", Core::LinearSolver::SolverType::umfpack, scatrasolvparams);
  scatrasolver_ = Teuchos::rcp(new Core::LinAlg::Solver(scatrasolvparams, firstscatradis->Comm()));
#else
  const Teuchos::ParameterList& fs3idyn = Global::Problem::Instance()->FS3IDynamicParams();
  // get solver number used for fs3i
  const int linsolvernumber = fs3idyn.get<int>("COUPLED_LINEAR_SOLVER");
  // check if LOMA solvers has a valid number
  if (linsolvernumber == (-1))
    FOUR_C_THROW(
        "no linear solver defined for FS3I problems. Please set COUPLED_LINEAR_SOLVER in FS3I "
        "DYNAMIC to a valid number!");
  const Teuchos::ParameterList& coupledscatrasolvparams =
      Global::Problem::Instance()->SolverParams(linsolvernumber);

  const auto solvertype =
      Teuchos::getIntegralValue<Core::LinearSolver::SolverType>(coupledscatrasolvparams, "SOLVER");

  if (solvertype != Core::LinearSolver::SolverType::belos)
    FOUR_C_THROW("Iterative solver expected");

  const auto azprectype = Teuchos::getIntegralValue<Core::LinearSolver::PreconditionerType>(
      coupledscatrasolvparams, "AZPREC");

  if (azprectype != Core::LinearSolver::PreconditionerType::block_gauss_seidel_2x2)
    FOUR_C_THROW("Block Gauss-Seidel preconditioner expected");

  // use coupled scatra solver object
  scatrasolver_ = Teuchos::rcp(new Core::LinAlg::Solver(coupledscatrasolvparams,
      firstscatradis->Comm(), Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY")));
  // get the solver number used for structural ScalarTransport solver
  const int linsolver1number = fs3idyn.get<int>("LINEAR_SOLVER1");
  // get the solver number used for structural ScalarTransport solver
  const int linsolver2number = fs3idyn.get<int>("LINEAR_SOLVER2");

  // check if the linear solver has a valid solver number
  if (linsolver1number == (-1))
    FOUR_C_THROW(
        "no linear solver defined for fluid ScalarTransport solver. Please set LINEAR_SOLVER2 in "
        "FS3I DYNAMIC to a valid number!");
  if (linsolver2number == (-1))
    FOUR_C_THROW(
        "no linear solver defined for structural ScalarTransport solver. Please set LINEAR_SOLVER2 "
        "in FS3I DYNAMIC to a valid number!");
  scatrasolver_->put_solver_params_to_sub_params("Inverse1",
      Global::Problem::Instance()->SolverParams(linsolver1number),
      Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY"));
  scatrasolver_->put_solver_params_to_sub_params("Inverse2",
      Global::Problem::Instance()->SolverParams(linsolver2number),
      Global::Problem::Instance()->solver_params_callback(),
      Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
          Global::Problem::Instance()->IOParams(), "VERBOSITY"));
  (scatravec_[0])
      ->ScaTraField()
      ->discretization()
      ->compute_null_space_if_necessary(scatrasolver_->Params().sublist("Inverse1"));
  (scatravec_[1])
      ->ScaTraField()
      ->discretization()
      ->compute_null_space_if_necessary(scatrasolver_->Params().sublist("Inverse2"));

#endif
}


/*----------------------------------------------------------------------*
 |  Test results                                          hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::TestResults(const Epetra_Comm& comm)
{
  Global::Problem::Instance()->AddFieldTest(fpsi_->fluid_field()->CreateFieldTest());

  fpsi_->poro_field()->structure_field()->CreateFieldTest();
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    Global::Problem::Instance()->AddFieldTest(scatra->create_sca_tra_field_test());
  }
  Global::Problem::Instance()->TestAll(comm);
}


/*----------------------------------------------------------------------*
 |  Transfer FPSI solution                                hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::SetFPSISolution()
{
  // we clear every state, including the states of the secondary dof sets
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    scatravec_[i]->ScaTraField()->discretization()->ClearState(true);
    // we have to manually clear this since this can not be saved directly in the
    // primary dof set (because it is cleared in between)
    scatravec_[i]->ScaTraField()->clear_external_concentrations();
  }

  set_mesh_disp();
  set_velocity_fields();
  set_wall_shear_stresses();
  SetPressureFields();
  set_membrane_concentration();
}

/*----------------------------------------------------------------------*
 |  Transfer scatra solution                              hemmler 07/14 |
 *----------------------------------------------------------------------*/
// only needed for two-way coupling; at the moment function is not used
void FS3I::PartFPS3I::set_struct_scatra_solution()
{
  fpsi_->poro_field()->structure_field()->discretization()->set_state(
      1, "scalarfield", (scatravec_[1])->ScaTraField()->Phinp());
}


/*----------------------------------------------------------------------*
 |  Set displacements                                     hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::set_mesh_disp()
{
  // fluid field
  scatravec_[0]->ScaTraField()->ApplyMeshMovement(fpsi_->fluid_field()->Dispnp());

  // Poro field
  scatravec_[1]->ScaTraField()->ApplyMeshMovement(fpsi_->poro_field()->structure_field()->Dispnp());
}


/*----------------------------------------------------------------------*
 |  Set velocities                                        hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::set_velocity_fields()
{
  Global::Problem* problem = Global::Problem::Instance();
  const Teuchos::ParameterList& scatradyn = problem->scalar_transport_dynamic_params();
  int cdvel = Core::UTILS::IntegralValue<int>(scatradyn, "VELOCITYFIELD");
  switch (cdvel)
  {
    case Inpar::ScaTra::velocity_zero:
    case Inpar::ScaTra::velocity_function:
    {
      for (unsigned i = 0; i < scatravec_.size(); ++i)
      {
        Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
        scatra->ScaTraField()->set_velocity_field();
      }
      break;
    }
    case Inpar::ScaTra::velocity_Navier_Stokes:
    {
      std::vector<Teuchos::RCP<const Epetra_Vector>> convel;
      std::vector<Teuchos::RCP<const Epetra_Vector>> vel;
      ExtractVel(convel, vel);

      for (unsigned i = 0; i < scatravec_.size(); ++i)
      {
        Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
        scatra->ScaTraField()->set_velocity_field(convel[i], Teuchos::null, vel[i], Teuchos::null);
      }
      break;
    }
  }
}

/*----------------------------------------------------------------------*
 |  Set wall shear stresses                               hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::set_wall_shear_stresses()
{
  std::vector<Teuchos::RCP<const Epetra_Vector>> wss;
  ExtractWSS(wss);

  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->set_wall_shear_stresses(wss[i]);
  }
}

/*----------------------------------------------------------------------*
 |  Set presures                                          hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::SetPressureFields()
{
  std::vector<Teuchos::RCP<const Epetra_Vector>> pressure;
  ExtractPressure(pressure);

  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->SetPressureField(pressure[i]);
  }
}

/*----------------------------------------------------------------------*
 |  Evaluate scatra fields                                hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::evaluate_scatra_fields()
{
  // membrane concentration at the interface needed for membrane equation of Kedem and Katchalsky.
  // NOTE: needs to be set here, since it depends on the scalar interface values on both
  // discretisations changing with each Newton iteration
  set_membrane_concentration();

  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra_adap = scatravec_[i];
    Teuchos::RCP<ScaTra::ScaTraTimIntImpl> scatra = scatra_adap->ScaTraField();

    scatra->PrepareLinearSolve();

    // add contributions due to finite interface permeability
    if (!infperm_)
    {
      Teuchos::RCP<Epetra_Vector> coupforce = scatracoupforce_[i];
      Teuchos::RCP<Core::LinAlg::SparseMatrix> coupmat = scatracoupmat_[i];

      coupforce->PutScalar(0.0);
      coupmat->Zero();

      // evaluate interface; second Kedem-Katchalsky equation for coupling of solute flux
      scatra->KedemKatchalsky(coupmat, coupforce);

      // apply Dirichlet boundary conditions to coupling matrix and vector
      const Teuchos::RCP<const Epetra_Map> dbcmap = scatra->DirichMaps()->CondMap();
      coupmat->ApplyDirichlet(*dbcmap, false);
      Core::LinAlg::apply_dirichlet_to_system(*coupforce, *scatrazeros_[i], *dbcmap);
    }
  }
}


/*----------------------------------------------------------------------*
 |  Extract velocities                                    hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::ExtractVel(std::vector<Teuchos::RCP<const Epetra_Vector>>& convel,
    std::vector<Teuchos::RCP<const Epetra_Vector>>& vel)
{
  // ############ Fluid Field ###############
  convel.push_back(fpsi_->fluid_field()->ConvectiveVel());
  vel.push_back(fpsi_->fluid_field()->Velnp());

  // ############ Poro Field ###############
  convel.push_back(fpsi_->poro_field()->fluid_field()->ConvectiveVel());
  vel.push_back(fpsi_->poro_field()->fluid_field()->Velnp());
}


/*----------------------------------------------------------------------*
 |  Extract wall shear stresses                           hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::ExtractWSS(std::vector<Teuchos::RCP<const Epetra_Vector>>& wss)
{
  // ############ Fluid Field ###############

  Teuchos::RCP<Adapter::FluidFSI> fluid =
      Teuchos::rcp_dynamic_cast<Adapter::FluidFSI>(fpsi_->fluid_field());
  if (fluid == Teuchos::null) FOUR_C_THROW("Dynamic cast to Adapter::FluidFSI failed!");

  Teuchos::RCP<Epetra_Vector> WallShearStress =
      fluid->calculate_wall_shear_stresses();  // CalcWallShearStress();
  wss.push_back(WallShearStress);

  // ############ Poro Field ###############

  // Hint: The Wall shear stresses in the fluid field at the Interface are equal to the ones of the
  // poro structure
  //      Therefore, we map the results of the wss of the fluid field to the dofs of the poro field
  //      without computing them explicitly in the poro field

  // extract FPSI-Interface from fluid field
  WallShearStress =
      fpsi_->FPSICoupl()->fluid_fpsi_vel_pres_extractor()->ExtractCondVector(WallShearStress);

  // replace global fluid interface dofs through porofluid interface dofs
  WallShearStress = fpsi_->FPSICoupl()->iFluidToPorofluid(WallShearStress);

  // insert porofluid interface entries into vector with full porofluid length
  Teuchos::RCP<Epetra_Vector> porofluid =
      Core::LinAlg::CreateVector(*(fpsi_->poro_field()->fluid_field()->dof_row_map()), true);

  // Parameter int block of function InsertVector:
  fpsi_->FPSICoupl()->poro_fluid_fpsi_vel_pres_extractor()->InsertVector(
      WallShearStress, 1, porofluid);

  wss.push_back(porofluid);
}

/*----------------------------------------------------------------------*
 |  Extract pressures                                     hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFPS3I::ExtractPressure(std::vector<Teuchos::RCP<const Epetra_Vector>>& pressure)
{
  // ############ Fluid Field ###############
  pressure.push_back(
      fpsi_->fluid_field()->Velnp());  // we extract the velocities as well. We sort them out later.

  // ############ Poro Field ###############
  pressure.push_back(fpsi_->poro_field()
                         ->fluid_field()
                         ->Velnp());  // we extract the velocities as well. We sort them out later.
}

FOUR_C_NAMESPACE_CLOSE
