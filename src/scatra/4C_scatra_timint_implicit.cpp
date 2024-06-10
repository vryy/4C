/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for convection-diffusion (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme

     o generalized-alpha time-integration scheme

     and stationary solver.

\level 1


*/
/*----------------------------------------------------------------------*/
#include "4C_scatra_timint_implicit.hpp"

#include "4C_comm_utils_gid_vector.hpp"
#include "4C_contact_nitsche_strategy_ssi.hpp"
#include "4C_fem_condition_periodic.hpp"
#include "4C_fem_condition_selector.hpp"
#include "4C_fem_general_assemblestrategy.hpp"
#include "4C_fluid_rotsym_periodicbc_utils.hpp"
#include "4C_fluid_turbulence_dyn_vreman.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_pstream.hpp"
#include "4C_linalg_krylov_projector.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_linear_solver_method_parameters.hpp"
#include "4C_mat_elchmat.hpp"
#include "4C_mat_electrode.hpp"
#include "4C_mat_list.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_scatra.hpp"
#include "4C_nurbs_discret.hpp"
#include "4C_nurbs_discret_apply_nurbs_initial_condition.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_boundary_calc_elch_electrode_utils.hpp"
#include "4C_scatra_ele_parameter_timint.hpp"
#include "4C_scatra_resulttest.hpp"
#include "4C_scatra_timint_heterogeneous_reaction_strategy.hpp"
#include "4C_scatra_timint_meshtying_strategy_artery.hpp"
#include "4C_scatra_timint_meshtying_strategy_fluid.hpp"
#include "4C_scatra_timint_meshtying_strategy_s2i.hpp"
#include "4C_scatra_timint_meshtying_strategy_std.hpp"
#include "4C_scatra_turbulence_hit_initial_scalar_field.hpp"
#include "4C_scatra_turbulence_hit_scalar_forcing.hpp"
#include "4C_scatra_utils.hpp"
#include "4C_ssi_contact_strategy.hpp"
#include "4C_utils_function.hpp"
#include "4C_utils_parameter_list.hpp"

#include <unordered_set>
#include <utility>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
ScaTra::ScaTraTimIntImpl::ScaTraTimIntImpl(Teuchos::RCP<Core::FE::Discretization> actdis,
    Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> extraparams,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output, const int probnum)
    : problem_(Global::Problem::Instance(probnum)),
      probnum_(probnum),
      solver_(std::move(solver)),
      params_(params),
      extraparams_(extraparams),
      myrank_(actdis->Comm().MyPID()),
      splitter_(Teuchos::null),
      strategy_(Teuchos::null),
      additional_model_evaluator_(nullptr),
      isale_(extraparams->get<bool>("isale")),
      solvtype_(Core::UTILS::IntegralValue<Inpar::ScaTra::SolverType>(*params, "SOLVERTYPE")),
      equilibrationmethod_(
          Teuchos::getIntegralValue<Core::LinAlg::EquilibrationMethod>(*params, "EQUILIBRATION")),
      matrixtype_(Teuchos::getIntegralValue<Core::LinAlg::MatrixType>(*params, "MATRIXTYPE")),
      incremental_(true),
      fssgd_(Core::UTILS::IntegralValue<Inpar::ScaTra::FSSUGRDIFF>(*params, "FSSUGRDIFF")),
      turbmodel_(Inpar::FLUID::no_model),
      s2ikinetics_(actdis->GetCondition("S2IKinetics") != nullptr),
      s2imeshtying_(actdis->GetCondition("S2IMeshtying") != nullptr),
      arterycoupling_(Core::UTILS::IntegralValue<int>(
                          problem_->poro_multi_phase_scatra_dynamic_params(), "ARTERY_COUPLING") &&
                      actdis->Name() == "scatra"),
      heteroreaccoupling_(actdis->GetCondition("ScatraHeteroReactionSlave") != nullptr),
      macro_scale_(
          problem_->Materials()->FirstIdByType(Core::Materials::m_scatra_multiscale) != -1 or
          problem_->Materials()->FirstIdByType(Core::Materials::m_newman_multiscale) != -1),
      micro_scale_(probnum != 0),
      isemd_(extraparams->get<bool>("ELECTROMAGNETICDIFFUSION", false)),
      emd_source_(extraparams->get<int>("EMDSOURCE", -1)),
      has_external_force_(
          Core::UTILS::IntegralValue<bool>(params_->sublist("EXTERNAL FORCE"), "EXTERNAL_FORCE")),
      calcflux_domain_(
          Core::UTILS::IntegralValue<Inpar::ScaTra::FluxType>(*params, "CALCFLUX_DOMAIN")),
      calcflux_domain_lumped_(Core::UTILS::IntegralValue<bool>(*params, "CALCFLUX_DOMAIN_LUMPED")),
      calcflux_boundary_(
          Core::UTILS::IntegralValue<Inpar::ScaTra::FluxType>(*params, "CALCFLUX_BOUNDARY")),
      calcflux_boundary_lumped_(
          Core::UTILS::IntegralValue<bool>(*params, "CALCFLUX_BOUNDARY_LUMPED")),
      writefluxids_(Teuchos::rcp(new std::vector<int>)),
      flux_domain_(Teuchos::null),
      flux_boundary_(Teuchos::null),
      flux_boundary_maps_(Teuchos::null),
      sumnormfluxintegral_(Teuchos::null),
      lastfluxoutputstep_(-1),
      outputscalars_(
          Core::UTILS::IntegralValue<Inpar::ScaTra::OutputScalarType>(*params, "OUTPUTSCALARS")),
      outputgmsh_(Core::UTILS::IntegralValue<int>(*params, "OUTPUT_GMSH")),
      output_state_matlab_(Core::UTILS::IntegralValue<int>(*params, "MATLAB_STATE_OUTPUT")),
      fdcheck_(Core::UTILS::IntegralValue<Inpar::ScaTra::FdCheck>(*params, "FDCHECK")),
      fdcheckeps_(params->get<double>("FDCHECKEPS")),
      fdchecktol_(params->get<double>("FDCHECKTOL")),
      computeintegrals_(
          Core::UTILS::IntegralValue<Inpar::ScaTra::ComputeIntegrals>(*params, "COMPUTEINTEGRALS")),
      calcerror_(Core::UTILS::IntegralValue<Inpar::ScaTra::CalcError>(*params, "CALCERROR")),
      time_(0.0),
      maxtime_(params->get<double>("MAXTIME")),
      step_(0),
      stepmax_(params->get<int>("NUMSTEP")),
      dta_(params->get<double>("TIMESTEP")),
      dtele_(0.0),
      dtsolve_(0.0),
      iternum_(0),
      iternum_outer_(0),
      timealgo_(
          Core::UTILS::IntegralValue<Inpar::ScaTra::TimeIntegrationScheme>(*params, "TIMEINTEGR")),
      nsd_(problem_->NDim()),
      scalarhandler_(Teuchos::null),
      outputscalarstrategy_(Teuchos::null),
      outputdomainintegralstrategy_(Teuchos::null),
      // Initialization of degrees of freedom variables
      phin_(Teuchos::null),
      phinp_(Teuchos::null),
      phinp_inc_(Teuchos::null),
      phinp_inc_old_(Teuchos::null),
      omega_(0, 0.),
      phidtn_(Teuchos::null),
      phidtnp_(Teuchos::null),
      hist_(Teuchos::null),
      densafnp_(Teuchos::null),
      velocity_field_type_(
          Core::UTILS::IntegralValue<Inpar::ScaTra::VelocityField>(*params, "VELOCITYFIELD")),
      mean_conc_(Teuchos::null),
      membrane_conc_(Teuchos::null),
      phinp_micro_(Teuchos::null),
      nds_disp_(-1),
      nds_growth_(-1),
      nds_micro_(-1),
      nds_pres_(-1),
      nds_scatra_(-1),
      nds_thermo_(-1),
      nds_two_tensor_quantitiy_(-1),
      nds_vel_(-1),
      nds_wss_(-1),
      densific_(0, 0.0),
      c0_(0, 0.0),
      macro_micro_rea_coeff_(0.0),
      discret_(actdis),
      output_(std::move(output)),
      convform_(Core::UTILS::IntegralValue<Inpar::ScaTra::ConvForm>(*params, "CONVFORM")),
      sysmat_(Teuchos::null),
      zeros_(Teuchos::null),
      dbcmaps_(Teuchos::null),
      neumann_loads_(Teuchos::null),
      normals_(Teuchos::null),
      residual_(Teuchos::null),
      trueresidual_(Teuchos::null),
      increment_(Teuchos::null),
      msht_(Core::UTILS::IntegralValue<Inpar::FLUID::MeshTying>(*params, "MESHTYING")),
      // Initialization of AVM3 variables
      sysmat_sd_(Teuchos::null),
      Sep_(Teuchos::null),
      Mnsv_(Teuchos::null),
      // Initialization of turbulent flow variables
      DynSmag_(Teuchos::null),
      Vrem_(Teuchos::null),
      samstart_(0),
      samstop_(0),
      dumperiod_(0),
      turbinflow_(Core::UTILS::IntegralValue<int>(
          extraparams->sublist("TURBULENT INFLOW"), "TURBULENTINFLOW")),
      numinflowsteps_(extraparams->sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")),
      special_flow_("initialization"),
      forcing_(Teuchos::null),
      homisoturb_forcing_(Teuchos::null),
      // Initialization of Krylov
      updateprojection_(false),
      projector_(Teuchos::null),
      // Initialization of
      upres_(params->get<int>("RESULTSEVRY")),
      uprestart_(params->get<int>("RESTARTEVRY")),
      neumanninflow_(Core::UTILS::IntegralValue<int>(*params, "NEUMANNINFLOW")),
      convheatrans_(Core::UTILS::IntegralValue<int>(*params, "CONV_HEAT_TRANS")),
      phinp_macro_(0, 0.),
      q_(0.0),
      dq_dphi_(0, 0.),
      // Initialization of Biofilm specific stuff
      scfldgrdisp_(Teuchos::null),
      scstrgrdisp_(Teuchos::null),
      outintegrreac_(Core::UTILS::IntegralValue<int>(*params, "OUTINTEGRREAC")),
      skipinitder_(Core::UTILS::IntegralValue<int>(*params, "SKIPINITDER")),
      timestepadapted_(false),
      issetup_(false),
      isinit_(false)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
}


/*------------------------------------------------------------------------*
 *------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::Init()
{
  set_is_setup(false);

  // -------------------------------------------------------------------
  // safety check for spherical coordinates
  // -------------------------------------------------------------------
  if (Core::UTILS::IntegralValue<bool>(*params_, "SPHERICALCOORDS") and nsd_ > 1)
    FOUR_C_THROW("Spherical coordinates only available for 1D problems!");

  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  // note: pbcs have to be correctly set up before extended ghosting is applied
  auto pbc = Teuchos::rcp(new Core::Conditions::PeriodicBoundaryConditions(discret_, false));
  if (pbc->HasPBC() and not isinit_)
  {
    pbc->update_dofs_for_periodic_boundary_conditions();
  }

  // -------------------------------------------------------------------
  // determine whether linear incremental or nonlinear solver
  // -------------------------------------------------------------------
  switch (solvtype_)
  {
    case Inpar::ScaTra::solvertype_nonlinear:
    case Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro:
    case Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken:
    case Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit:
    case Inpar::ScaTra::solvertype_nonlinear_multiscale_microtomacro:
    case Inpar::ScaTra::solvertype_linear_incremental:
    {
      incremental_ = true;
    }
    break;
    case Inpar::ScaTra::solvertype_linear_full:
    {
      incremental_ = false;
    }
    break;
    default:
      FOUR_C_THROW("Received illegal scatra solvertype enum.");
      break;
  }

  // -----------------------------------------------------------------------
  // determine number of degrees of freedom and transported scalars per node
  // -----------------------------------------------------------------------
  create_scalar_handler();

  // -------------------------------------------------------------------
  // check compatibility of boundary conditions
  // -------------------------------------------------------------------
  if (neumanninflow_ and convheatrans_)
  {
    FOUR_C_THROW(
        "Neumann inflow and convective heat transfer boundary conditions must not appear "
        "simultaneously for the same problem!");
  }

  // -----------------------------------------------------------------------------
  // initialize meshtying strategy (including standard strategy without meshtying)
  // -----------------------------------------------------------------------------
  // safety checks
  if (msht_ != Inpar::FLUID::no_meshtying and s2imeshtying_)
  {
    FOUR_C_THROW(
        "Fluid-fluid meshtying in combination with scatra-scatra interface mesh tying is not "
        "implemented yet!");
  }
  if (s2imeshtying_ and !incremental_)
  {
    FOUR_C_THROW(
        "Scatra-scatra interface mesh tying only working for incremental solve so far!\n"
        "Set the parameter SOLVERTYPE in SCALAR TRANSPORT DYNAMIC section to 'nonlinear' or "
        "'linear_incremental'!");
  }

  ScaTraUtils::CheckConsistencyOfS2IConditions(discretization());

  // create strategy
  create_meshtying_strategy();

  // initialize strategy
  strategy_->InitMeshtying();

  // we have successfully initialized this class
  set_is_init(true);
}  // ScaTraTimIntImpl::Init()


/*------------------------------------------------------------------------*
 *------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::Setup()
{
  // we have to call Init() first
  check_is_init();

  // compute Null Space
  compute_null_space_if_necessary();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices: local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  // initialize the scalar handler
  if (scalarhandler_ == Teuchos::null)
    FOUR_C_THROW("Make sure you construct the scalarhandler_ in initialization.");
  else
    scalarhandler_->Setup(this);

  // setup splitter (needed to solve initialization problems before setup_meshtying())
  SetupSplitter();

  // setup the matrix block maps and the meshtying strategy
  setup_matrix_block_maps_and_meshtying();

  // -------------------------------------------------------------------
  // create empty system matrix (27 adjacent nodes as 'good' guess)
  // -------------------------------------------------------------------
  if (fssgd_ != Inpar::ScaTra::fssugrdiff_no and not incremental_)
  {
    // do not save the graph if fine-scale subgrid diffusivity is used in non-incremental case (very
    // special case)
    sysmat_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(discret_->dof_row_map()), 27));
  }
  else
    sysmat_ = init_system_matrix();

  // for some special meshtying cases we need to override the information with the information from
  // the meshtying strategy, in such a case the meshtying strategy implements the method below
  if (strategy_->system_matrix_initialization_needed()) sysmat_ = strategy_->init_system_matrix();

  // -------------------------------------------------------------------
  // create vectors containing problem variables
  // -------------------------------------------------------------------
  // solutions at time n+1 and n
  phinp_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  phin_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  if (NdsMicro() != -1)
    phinp_micro_ = Core::LinAlg::CreateVector(*discret_->dof_row_map(NdsMicro()));

  if (solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro or
      solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken or
      solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit or
      solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_microtomacro)
  {
    phinp_inc_ = Core::LinAlg::CreateVector(*dofrowmap, true);
    if (solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken or
        solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit)
    {
      phinp_inc_old_ = Core::LinAlg::CreateVector(*dofrowmap, true);
      if (solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken)
        omega_.resize(1, 1.);
      else
        omega_.resize(NumDofPerNode(), 1.);
    }
  }

  // temporal solution derivative at time n+1
  phidtnp_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  // temporal solution derivative at time n
  phidtn_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // history vector (a linear combination of phinm, phin (BDF)
  // or phin, phidtn (One-Step-Theta, Generalized-alpha))
  hist_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // -------------------------------------------------------------------
  // create vectors associated to boundary conditions
  // -------------------------------------------------------------------
  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new Core::LinAlg::MapExtractor());
  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time", time_);
    eleparams.set<const Core::UTILS::FunctionManager*>(
        "function_manager", &Global::Problem::Instance()->FunctionManager());
    const Core::ProblemType problem_type = Core::ProblemType::scatra;
    eleparams.set<const Core::ProblemType*>("problem_type", &problem_type);
    discret_->evaluate_dirichlet(
        eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0);  // just in case of change
  }

  // -------------------------------------------------------------------
  // create vectors associated to solution process
  // -------------------------------------------------------------------
  // the vector containing body and surface forces
  neumann_loads_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // the residual vector --- more or less the rhs
  residual_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // residual vector containing the normal boundary fluxes
  trueresidual_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // incremental solution vector
  increment_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // subgrid-diffusivity(-scaling) vector
  // (used either for AVM3 approach or temperature equation
  //  with all-scale subgrid-diffusivity model)
  if (fssgd_ != Inpar::ScaTra::fssugrdiff_no)
    subgrdiff_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // -------------------------------------------------------------------
  // set parameters associated to potential statistical flux evaluations
  // -------------------------------------------------------------------
  // initialize vector for statistics (assume a maximum of 10 conditions)
  sumnormfluxintegral_ = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(10));

  if (calcflux_domain_ != Inpar::ScaTra::flux_none or
      calcflux_boundary_ != Inpar::ScaTra::flux_none)
  {
    // safety check
    if (not scalarhandler_->EqualNumDof())
    {
      FOUR_C_THROW(
          "Flux output only implement for equal number of DOFs per node within ScaTra "
          "discretization!");
    }

    // if the writefluxids vector has not been set yet
    if (writefluxids_->empty())
    {
      // write one by one of scalars (as flux output in the input file defined)
      // to the temporary variable word1
      int word1 = 0;
      std::istringstream mystream(Teuchos::getNumericStringParameter(*params_, "WRITEFLUX_IDS"));
      while (mystream >> word1)
        // get desired scalar id's for flux output
        writefluxids_->push_back(word1);
    }

    // default value (-1): flux is written for all dof's
    // scalar transport: numdofpernode_ = numscal_
    // elch:             numdofpernode_ = numscal_+1
    // -> current flux for potential only if div i is used to close the system otherwise zero
    if (writefluxids_->size() and
        (*writefluxids_)[0] == (-1))  // default is to perform flux output for ALL scalars
    {
      writefluxids_->resize(NumDofPerNode());
      for (int k = 0; k < NumDofPerNode(); ++k) (*writefluxids_)[k] = k + 1;
    }

    // flux_ vector is initialized when CalcFlux() is called

    // screen output
    if (myrank_ == 0)
    {
      Core::IO::cout << "Flux output is performed for " << writefluxids_->size() << " scalars: ";
      for (unsigned int i = 0; i < writefluxids_->size(); i++)
      {
        const int id = (*writefluxids_)[i];
        Core::IO::cout << id << " ";
        if ((id < 1) or (id > NumDofPerNode()))  // check validity of these numbers as well !
          FOUR_C_THROW("Received illegal scalar id for flux output: %d", id);
      }
      Core::IO::cout << Core::IO::endl;
    }

    // initialize map extractor associated with boundary segments for flux calculation
    if (calcflux_boundary_ != Inpar::ScaTra::flux_none)
    {
      // extract conditions for boundary flux calculation
      std::vector<Core::Conditions::Condition*> conditions;
      discret_->GetCondition("ScaTraFluxCalc", conditions);

      // set up map extractor
      flux_boundary_maps_ = Teuchos::rcp(new Core::LinAlg::MultiMapExtractor());
      Core::Conditions::MultiConditionSelector mcs;
      mcs.SetOverlapping(true);
      for (auto& condition : conditions)
        mcs.AddSelector(Teuchos::rcp(new Core::Conditions::ConditionSelector(
            *discret_, std::vector<Core::Conditions::Condition*>(1, condition))));
      mcs.SetupExtractor(*discret_, *discret_->dof_row_map(), *flux_boundary_maps_);
    }
  }

  // -------------------------------------------------------------------
  // preparations for turbulence models
  // -------------------------------------------------------------------
  const Epetra_Map* noderowmap = discret_->NodeRowMap();
  InitTurbulenceModel(dofrowmap, noderowmap);

  // -------------------------------------------------------------------
  // set initial field
  // -------------------------------------------------------------------
  SetInitialField(Core::UTILS::IntegralValue<Inpar::ScaTra::InitialField>(*params_, "INITIALFIELD"),
      params_->get<int>("INITFUNCNO"));

  // -------------------------------------------------------------------
  // preparations for natural convection
  // -------------------------------------------------------------------
  if (static_cast<bool>(Core::UTILS::IntegralValue<int>(*params_, "NATURAL_CONVECTION")))
  {
    // allocate global density vector and initialize
    densafnp_ = Core::LinAlg::CreateVector(*discret_->dof_row_map(), true);
    densafnp_->PutScalar(1.);
  }

  // -------------------------------------------------------------------
  // preparations for total and mean values of transported scalars
  // -------------------------------------------------------------------
  if (outputscalars_ != Inpar::ScaTra::outputscalars_none)
  {
    // input check
    if (outputscalars_ == Inpar::ScaTra::outputscalars_entiredomain)
    {
      std::vector<Core::Conditions::Condition*> conditions;
      // extract conditions for calculation of total and mean values of transported scalars
      discret_->GetCondition("TotalAndMeanScalar", conditions);
      // input check
      if (conditions.size())
      {
        FOUR_C_THROW(
            "Found 'DESIGN TOTAL AND MEAN SCALAR' condition on ScaTra discretization, but "
            "'OUTPUTSCALAR' \n"
            "in 'SCALAR TRANSPORT DYNAMIC' is set to 'entire domain'. Either switch on the output "
            "of mean and total scalars\n"
            "on conditions or remoove the 'DESIGN TOTAL AND MEAN SCALAR' condition from your input "
            "file!");
      }
    }

    // build helper class for total and mean scalar output depending on input parameter
    switch (outputscalars_)
    {
      case Inpar::ScaTra::outputscalars_entiredomain:
        outputscalarstrategy_ = Teuchos::rcp(new OutputScalarsStrategyDomain);
        break;
      case Inpar::ScaTra::outputscalars_condition:
        outputscalarstrategy_ = Teuchos::rcp(new OutputScalarsStrategyCondition);
        break;
      case Inpar::ScaTra::outputscalars_entiredomain_condition:
        outputscalarstrategy_ = Teuchos::rcp(new OutputScalarsStrategyDomainAndCondition);
        break;
      default:
        FOUR_C_THROW("Unknown option for output of total and mean scalars!");
        break;
    }

    // initialize scalar output strategy
    outputscalarstrategy_->Init(this);
  }
  else
  {
    // input check

    std::vector<Core::Conditions::Condition*> conditions;
    // extract conditions for calculation of total and mean values of transported scalars
    discret_->GetCondition("TotalAndMeanScalar", conditions);
    // input check
    if (conditions.size())
    {
      FOUR_C_THROW(
          "Found 'DESIGN TOTAL AND MEAN SCALAR' condition on ScaTra discretization, but "
          "'OUTPUTSCALAR' \n"
          "in 'SCALAR TRANSPORT DYNAMIC' is set to 'none'. Either switch on the output of mean and "
          "total scalars\n"
          "or remove the 'DESIGN TOTAL AND MEAN SCALAR' condition from your input file!");
    }
  }

  // -------------------------------------------------------------------
  // preparations for domain integrals
  // -------------------------------------------------------------------
  if (computeintegrals_ != Inpar::ScaTra::computeintegrals_none)
  {
    // initialize domain integral output strategy
    outputdomainintegralstrategy_ = Teuchos::rcp(new OutputDomainIntegralStrategy);
    outputdomainintegralstrategy_->Init(this);
  }
  else
  {
    // input check
    std::vector<Core::Conditions::Condition*> conditions_boundary;
    // extract conditions for calculation of total and mean values of transported scalars
    discret_->GetCondition("BoundaryIntegral", conditions_boundary);
    std::vector<Core::Conditions::Condition*> conditions_domain;
    // extract conditions for calculation of total and mean values of transported scalars
    discret_->GetCondition("DomainIntegral", conditions_domain);
    // input check
    if (conditions_boundary.size() > 0 || conditions_domain.size() > 0)
    {
      FOUR_C_THROW(
          "Found 'DESIGN DOMAIN INTEGRAL SURF CONDITIONS' or 'DESIGN DOMAIN INTEGRAL VOL "
          "CONDITIONS' condition on ScaTra discretization, but COMPUTEINTEGRALS\n"
          "in 'SCALAR TRANSPORT DYNAMIC' is set to 'none'. Either switch on the output of domain "
          "integrals "
          "or remove the 'DESIGN DOMAIN INTEGRAL * CONDITIONS' condition from your input file!");
    }
  }

  // -------------------------------------------------------------------
  // preparations for error evaluation
  // -------------------------------------------------------------------
  if (calcerror_ != Inpar::ScaTra::calcerror_no)
  {
    if (calcerror_ == Inpar::ScaTra::calcerror_bycondition)
    {
      std::vector<Core::Conditions::Condition*> relerrorconditions;
      discret_->GetCondition("ScatraRelError", relerrorconditions);
      const unsigned ncond = relerrorconditions.size();
      if (!ncond)
      {
        FOUR_C_THROW(
            "Calculation of relative error based on conditions desired, but no conditions "
            "specified!");
      }
      relerrors_ =
          Teuchos::rcp(new std::vector<double>(2 * NumDofPerNode() * relerrorconditions.size()));
    }
    else if (calcerror_ == Inpar::ScaTra::calcerror_AnalyticSeries)
      relerrors_ = Teuchos::rcp(new std::vector<double>(2));  // TODO: Update two n species
    else
    {
      // It is important to make a distinction as HDG always have NumDofPerNode = 0
      // The vector is therefore sized to contain the errors of one scalar and its gradient
      if (Global::Problem::Instance()->spatial_approximation_type() ==
          Core::FE::ShapeFunctionType::hdg)
        relerrors_ = Teuchos::rcp(new std::vector<double>(2));  // TODO: update to n species
      else
        relerrors_ = Teuchos::rcp(new std::vector<double>(2 * NumDofPerNode()));
    }
  }

  // we have successfully set up this class
  set_is_setup(true);
}  // ScaTraTimIntImpl::Setup()


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::SetupNatConv()
{
  // calculate the initial mean concentration value
  if (NumScal() < 1) FOUR_C_THROW("Error since numscal = %d. Not allowed since < 1", NumScal());
  c0_.resize(NumScal());

  discret_->set_state("phinp", phinp_);

  // set action for elements
  Teuchos::ParameterList eleparams;
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::calc_total_and_mean_scalars, eleparams);
  eleparams.set("inverting", false);
  eleparams.set("calc_grad_phi", false);

  // evaluate integrals of concentrations and domain
  Teuchos::RCP<Core::LinAlg::SerialDenseVector> scalars =
      Teuchos::rcp(new Core::LinAlg::SerialDenseVector(NumScal() + 1));
  discret_->EvaluateScalars(eleparams, scalars);

  // calculate mean concentrations
  const double domint = (*scalars)[NumScal()];
  if (std::abs(domint) < 1e-15) FOUR_C_THROW("Domain has zero volume!");
  for (int k = 0; k < NumScal(); ++k) c0_[k] = (*scalars)[k] / domint;

  // initialization of the densification coefficient vector
  densific_.resize(NumScal());
  Core::Elements::Element* element = discret_->lRowElement(0);
  Teuchos::RCP<Core::Mat::Material> mat = element->Material();

  if (mat->MaterialType() == Core::Materials::m_matlist or
      mat->MaterialType() == Core::Materials::m_matlist_reactions)
  {
    Teuchos::RCP<const Mat::MatList> actmat = Teuchos::rcp_static_cast<const Mat::MatList>(mat);

    for (int k = 0; k < NumScal(); ++k)
    {
      const int matid = actmat->MatID(k);
      Teuchos::RCP<const Core::Mat::Material> singlemat = actmat->MaterialById(matid);

      if (singlemat->MaterialType() == Core::Materials::m_scatra)
      {
        Teuchos::RCP<const Mat::ScatraMat> actsinglemat =
            Teuchos::rcp_static_cast<const Mat::ScatraMat>(singlemat);

        densific_[k] = actsinglemat->Densification();

        if (densific_[k] < 0.0) FOUR_C_THROW("received negative densification value");
      }
      else
        FOUR_C_THROW("Material type is not allowed!");
    }
  }

  // for a single species calculation
  else if (mat->MaterialType() == Core::Materials::m_scatra)
  {
    Teuchos::RCP<const Mat::ScatraMat> actmat = Teuchos::rcp_static_cast<const Mat::ScatraMat>(mat);

    densific_[0] = actmat->Densification();

    if (densific_[0] < 0.0) FOUR_C_THROW("received negative densification value");
    if (NumScal() > 1) FOUR_C_THROW("Single species calculation but numscal = %d > 1", NumScal());
  }
  else
    FOUR_C_THROW("Material type is not allowed!");
}  // ScaTraTimIntImpl::SetupNatConv


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::InitTurbulenceModel(
    const Epetra_Map* dofrowmap, const Epetra_Map* noderowmap)
{
  // get fluid turbulence sublist
  Teuchos::ParameterList* turbparams = &(extraparams_->sublist("TURBULENCE MODEL"));

  // parameters for statistical evaluation of normal fluxes
  samstart_ = turbparams->get<int>("SAMPLING_START");
  samstop_ = turbparams->get<int>("SAMPLING_STOP");
  dumperiod_ = turbparams->get<int>("DUMPING_PERIOD");
  if (dumperiod_ < 0) FOUR_C_THROW("dumperiod_ is negative!");

  // -------------------------------------------------------------------
  // necessary only for AVM3 approach:
  // initialize subgrid-diffusivity matrix + respective output
  // -------------------------------------------------------------------
  if (fssgd_ != Inpar::ScaTra::fssugrdiff_no and
      Core::UTILS::IntegralValue<int>(*turbparams, "TURBMODEL_LS"))
  {
    sysmat_sd_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dofrowmap, 27));

    // Output
    if (myrank_ == 0)
    {
      std::cout << "SCATRA: Fine-scale subgrid-diffusivity approach based on AVM3: ";
      std::cout << fssgd_;
      std::cout << " with turbulent Prandtl number Prt= ";
      std::cout << extraparams_->sublist("SUBGRID VISCOSITY").get<double>("C_TURBPRANDTL");
      std::cout << &std::endl << &std::endl;
    }

    if (turbparams->get<std::string>("PHYSICAL_MODEL") != "Multifractal_Subgrid_Scales")
    {
      if (fssgd_ == Inpar::ScaTra::fssugrdiff_smagorinsky_small and
          turbparams->get<std::string>("FSSUGRVISC") != "Smagorinsky_small")
        FOUR_C_THROW("Same subgrid-viscosity approach expected!");
      if (fssgd_ == Inpar::ScaTra::fssugrdiff_smagorinsky_all and
          turbparams->get<std::string>("FSSUGRVISC") != "Smagorinsky_all")
        FOUR_C_THROW("Same subgrid-viscosity approach expected!");
    }
  }
  else
    fssgd_ = Inpar::ScaTra::fssugrdiff_no;  // in case of not "TURBMODEL_LS"

  // -------------------------------------------------------------------
  // get turbulence model and parameters
  // -------------------------------------------------------------------
  turbmodel_ = Inpar::FLUID::no_model;

  if (Core::UTILS::IntegralValue<int>(*turbparams, "TURBMODEL_LS"))
  {
    // set turbulence model
    if (turbparams->get<std::string>("PHYSICAL_MODEL") == "Smagorinsky")
    {
      turbmodel_ = Inpar::FLUID::smagorinsky;

      // Output
      if (turbmodel_ and myrank_ == 0)
      {
        std::cout << "All-scale subgrid-diffusivity model: ";
        std::cout << turbparams->get<std::string>("PHYSICAL_MODEL");
        std::cout << &std::endl << &std::endl;
      }
    }
    else if (turbparams->get<std::string>("PHYSICAL_MODEL") == "Dynamic_Smagorinsky")
    {
      turbmodel_ = Inpar::FLUID::dynamic_smagorinsky;
      // access to the dynamic Smagorinsky class will provided by the
      // scatra fluid couling algorithm
    }
    else if (turbparams->get<std::string>("PHYSICAL_MODEL") == "Multifractal_Subgrid_Scales")
    {
      turbmodel_ = Inpar::FLUID::multifractal_subgrid_scales;

      // initalize matrix used to build the scale separation operator
      sysmat_sd_ = Teuchos::rcp(new Core::LinAlg::SparseMatrix(*dofrowmap, 27));

      Teuchos::ParameterList* mfsparams = &(extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES"));
      if (mfsparams->get<std::string>("SCALE_SEPARATION") != "algebraic_multigrid_operator")
        FOUR_C_THROW("Only scale separation by plain algebraic multigrid available in scatra!");

      // Output
      if (turbmodel_ and myrank_ == 0)
      {
        std::cout << "Multifractal subgrid-scale model: ";
        std::cout << turbparams->get<std::string>("PHYSICAL_MODEL");
        std::cout << &std::endl << &std::endl;
      }
    }
    else if (turbparams->get<std::string>("PHYSICAL_MODEL") == "Dynamic_Vreman")
    {
      // equivalent to dynamic smagorinsky
      turbmodel_ = Inpar::FLUID::dynamic_vreman;
    }

    // warning No. 1: if classical (all-scale) turbulence model other than
    // Smagorinsky or multifractal subrgid-scale modeling
    // is intended to be used
    if (turbparams->get<std::string>("PHYSICAL_MODEL") != "Smagorinsky" and
        turbparams->get<std::string>("PHYSICAL_MODEL") != "Dynamic_Smagorinsky" and
        turbparams->get<std::string>("PHYSICAL_MODEL") != "Multifractal_Subgrid_Scales" and
        turbparams->get<std::string>("PHYSICAL_MODEL") != "Dynamic_Vreman" and
        turbparams->get<std::string>("PHYSICAL_MODEL") != "no_model")
    {
      FOUR_C_THROW(
          "No classical (all-scale) turbulence model other than constant-coefficient Smagorinsky "
          "model and multifractal subrgid-scale modeling currently possible!");
    }

    // warning No. 2: if classical (all-scale) turbulence model and fine-scale
    // subgrid-viscosity approach are intended to be used simultaneously
    if (turbmodel_ == Inpar::FLUID::smagorinsky and fssgd_ != Inpar::ScaTra::fssugrdiff_no)
    {
      FOUR_C_THROW(
          "No combination of classical turbulence model and fine-scale subgrid-diffusivity "
          "approach currently possible!");
    }
  }

  if (turbmodel_ != Inpar::FLUID::no_model and NumScal() > 1)
    FOUR_C_THROW("Turbulent passive scalar transport not supported for more than one scalar!");

  // -------------------------------------------------------------------
  // initialize forcing for homogeneous isotropic turbulence
  // -------------------------------------------------------------------
  // flag for special flow
  special_flow_ =
      extraparams_->sublist("TURBULENCE MODEL").get<std::string>("CANONICAL_FLOW", "no");
  if (special_flow_ == "scatra_forced_homogeneous_isotropic_turbulence")
  {
    if (extraparams_->sublist("TURBULENCE MODEL").get<std::string>("SCALAR_FORCING") == "isotropic")
    {
      forcing_ = Core::LinAlg::CreateVector(*dofrowmap, true);
      forcing_->PutScalar(0.0);
    }
  }
}  // ScaTraTimIntImpl::InitTurbulenceModel()


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::prepare_krylov_projection()
{
  // sysmat might be singular (some modes are defined only up to a constant)
  // in this case, we need basis vectors for the nullspace/kernel

  // get condition "KrylovSpaceProjection" from discretization
  std::vector<Core::Conditions::Condition*> KSPCond;
  discret_->GetCondition("KrylovSpaceProjection", KSPCond);
  std::size_t numcond = KSPCond.size();
  int numscatra = 0;

  Core::Conditions::Condition* kspcond = nullptr;
  // check if for scatra Krylov projection is required
  for (std::size_t icond = 0; icond < numcond; icond++)
  {
    const auto& name = KSPCond[icond]->parameters().Get<std::string>("discretization");
    if (name == "scatra")
    {
      numscatra++;
      kspcond = KSPCond[icond];
    }
  }

  // initialize variables for Krylov projection if necessary
  if (numscatra == 1)
  {
    setup_krylov_space_projection(kspcond);
    if (myrank_ == 0)
      std::cout << "\nSetup of KrylovSpaceProjection in scatra field\n" << std::endl;
  }
  else if (numscatra == 0)
  {
    projector_ = Teuchos::null;
  }
  else
    FOUR_C_THROW("Received more than one KrylovSpaceCondition for scatra field");
}

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::set_element_nodeset_parameters() const
{
  Teuchos::ParameterList eleparams;

  // set action
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::set_nodeset_parameter, eleparams);

  eleparams.set<int>("ndsdisp", NdsDisp());
  eleparams.set<int>("ndsgrowth", NdsGrowth());
  eleparams.set<int>("ndspres", NdsPressure());
  eleparams.set<int>("ndsscatra", NdsScaTra());
  eleparams.set<int>("ndsthermo", NdsThermo());
  eleparams.set<int>("ndsTwoTensorQuantity", nds_two_tensor_quantity());
  eleparams.set<int>("ndsvel", NdsVel());
  eleparams.set<int>("ndswss", NdsWallShearStress());

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}

/*--------------------------------------------------------------------------------*
 *--------------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::set_element_general_parameters(bool calcinitialtimederivative) const
{
  Teuchos::ParameterList eleparams;

  // set action
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::set_general_scatra_parameter, eleparams);

  // set problem number
  eleparams.set<int>("probnum", probnum_);

  eleparams.set<int>("convform", convform_);

  eleparams.set<bool>("isale", isale_);

  // set flag for writing the flux vector fields
  eleparams.set<int>("calcflux_domain", calcflux_domain_);

  //! set vector containing ids of scalars for which flux vectors are calculated
  eleparams.set<Teuchos::RCP<std::vector<int>>>("writeflux_ids", writefluxids_);

  // parameters for stabilization
  eleparams.sublist("stabilization") = params_->sublist("STABILIZATION");
  if (calcinitialtimederivative)
  {
    // deactivate stabilization when calculating initial time derivative
    eleparams.sublist("stabilization").set<std::string>("STABTYPE", "no_stabilization");
    eleparams.sublist("stabilization").set<std::string>("DEFINITION_TAU", "Zero");
    // deactivate subgrid-scale velocity
    eleparams.sublist("stabilization").set<std::string>("SUGRVEL", "no");
    // deactivate subgrid diffusivity
    eleparams.sublist("stabilization").set<std::string>("ASSUGRDIFF", "no");
  }

  // parameters for finite difference check
  if (calcinitialtimederivative)
  {  // deactivate finite difference check when calculating initial
     // time derivative
    eleparams.set<int>("fdcheck", Inpar::ScaTra::fdcheck_none);
  }
  else
    eleparams.set<int>("fdcheck", fdcheck_);

  eleparams.set<double>("fdcheckeps", fdcheckeps_);
  eleparams.set<double>("fdchecktol", fdchecktol_);

  // flag for spherical coordinates
  eleparams.set<bool>(
      "sphericalcoords", Core::UTILS::IntegralValue<bool>(*params_, "SPHERICALCOORDS"));

  // flag for truly partitioned multi-scale simulation
  eleparams.set<bool>("partitioned_multiscale",
      solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro or
          solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken or
          solvtype_ ==
              Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit or
          solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_microtomacro);

  // flag for electromagnetic diffusion
  eleparams.set<bool>("electromagnetic_diffusion", isemd_);
  // current source function
  if (isemd_) eleparams.set<int>("electromagnetic_diffusion_source", emd_source_);

  // flag for external force
  eleparams.set<bool>("has_external_force", has_external_force_);

  // add parameters associated with meshtying strategy
  strategy_->set_element_general_parameters(eleparams);

  // additional problem-specific parameters for non-standard scalar transport problems
  // (electrochemistry etc.)
  set_element_specific_sca_tra_parameters(eleparams);

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::set_element_turbulence_parameters(
    bool calcinitialtimederivative) const
{
  Teuchos::ParameterList eleparams;

  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::set_turbulence_scatra_parameter, eleparams);

  eleparams.sublist("TURBULENCE MODEL") = extraparams_->sublist("TURBULENCE MODEL");
  if (calcinitialtimederivative)
  {
    // deactivate turbulence model when calculating initial time
    // derivative
    Teuchos::setStringToIntegralParameter<int>("PHYSICAL_MODEL", "no_model",
        "Classical LES approaches require an additional model for\nthe turbulent viscosity.",
        Teuchos::tuple<std::string>("no_model"),
        Teuchos::tuple<std::string>("If classical LES is our turbulence approach, this is a "
                                    "contradiction and should cause a FOUR_C_THROW."),
        Teuchos::tuple<int>(Inpar::FLUID::no_model), &eleparams.sublist("TURBULENCE MODEL"));
  }

  // set model-dependent parameters
  eleparams.sublist("SUBGRID VISCOSITY") = extraparams_->sublist("SUBGRID VISCOSITY");

  // set parameters for multifractal subgrid-scale modeling
  eleparams.sublist("MULTIFRACTAL SUBGRID SCALES") =
      extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES");

  eleparams.set<bool>("turbulent inflow", turbinflow_);

  if (calcinitialtimederivative)
    eleparams.set<int>("fs subgrid diffusivity", Inpar::ScaTra::fssugrdiff_no);
  else
    eleparams.set<int>("fs subgrid diffusivity", fssgd_);

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::prepare_time_loop()
{
  // provide information about initial field (do not do for restarts!)
  if (step_ == 0)
  {
    // write out initial state
    check_and_write_output_and_restart();

    // compute error for problems with analytical solution (initial field!)
    evaluate_error_compared_to_analytical_sol();

    // calculate mean concentration of micro discretization and set state to nds_micro_
    if (macro_scale_ and NdsMicro() != -1) calc_mean_micro_concentration();
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::prepare_time_step()
{
  // time measurement: prepare time step
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + prepare time step");

  // -------------------------------------------------------------------
  //                       initialization
  // -------------------------------------------------------------------
  if (step_ == 0) prepare_first_time_step();

  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  // adapt time step size if desired
  adapt_time_step_size();

  // tell micro scale about updated time step
  if (macro_scale_ and TimeStepAdapted()) set_time_stepping_to_micro_scale();

  // note the order of the following three functions is important
  increment_time_and_step();

  // -------------------------------------------------------------------
  // set part of the rhs vector belonging to the old timestep
  // -------------------------------------------------------------------
  set_old_part_of_righthandside();
  // TODO (Thon): We do not really want to call set_element_time_parameter() every time step.
  // But for now we just do it since "total time" has to be changed in the parameter class..
  set_element_time_parameter();

  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  // TODO: Dirichlet auch im Fall von genalpha phinp
  // Neumann(n + alpha_f)
  apply_dirichlet_bc(time_, phinp_, Teuchos::null);
  apply_neumann_bc(neumann_loads_);

  // By definition: Applying DC on the slave side of an internal interface is not allowed
  //                since it leads to an over-constraint system
  // Therefore, nodes belonging to the slave side of an internal interface have to be excluded from
  // the DC. However, a velocity value (projected from the Dirichlet condition on the master side)
  // has to be assigned to the DOF's on the slave side in order to evaluate the system matrix
  // completely

  // Preparation for including DC on the master side in the condensation process
  strategy_->include_dirichlet_in_condensation();

  // -------------------------------------------------------------------
  //     update velocity field if given by function (it might depend on time)
  // -------------------------------------------------------------------
  if (velocity_field_type_ == Inpar::ScaTra::velocity_function) set_velocity_field();

  // -------------------------------------------------------------------
  //     update external force given by function (it might depend on time)
  // -------------------------------------------------------------------
  if (has_external_force_) SetExternalForce();

  // -------------------------------------------------------------------
  //           preparation of AVM3-based scale separation
  // -------------------------------------------------------------------
  if ((step_ == 1 or (turbinflow_ and step_ == numinflowsteps_ + 1)) and
      (fssgd_ != Inpar::ScaTra::fssugrdiff_no or
          turbmodel_ == Inpar::FLUID::multifractal_subgrid_scales))
    av_m3_preparation();

  // -------------------------------------------------------------------
  // compute values at intermediate time steps (only for gen.-alpha)
  // -------------------------------------------------------------------
  compute_intermediate_values();

  // -------------------------------------------------------------------
  // prepare time step on micro scale if necessary
  // -------------------------------------------------------------------
  if (macro_scale_)
  {
    // create parameter list for macro-scale elements
    Teuchos::ParameterList eleparams;

    // set action
    Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
        "action", ScaTra::Action::micro_scale_prepare_time_step, eleparams);

    // add state vectors
    add_time_integration_specific_vectors();

    // loop over macro-scale elements
    discret_->Evaluate(
        eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  }
}

/*------------------------------------------------------------------------------*
 *------------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::prepare_first_time_step()
{
  if (not skipinitder_)
  {
    if (NdsVel() != -1)  // if some velocity field has been set
    {
      // TODO: Restructure enforcement of Dirichlet boundary conditions on phin_
      // A clean solution would incorporate apply_dirichlet_bc(...) into
      // calc_initial_time_derivative(). However, this would make a number of test cases fail. We
      // should have a closer look at this problem and fix it eventually.
      apply_dirichlet_bc(time_, phin_, Teuchos::null);
      calc_initial_time_derivative();
    }

    // if initial velocity field has not been set here, the initial time derivative of phi will be
    // calculated wrongly for some time integration schemes
    else
      FOUR_C_THROW("Initial velocity field has not been set!");
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::PrepareLinearSolve()
{
  // special preparations for multifractal subgrid-scale model
  if (turbmodel_ == Inpar::FLUID::multifractal_subgrid_scales) recompute_mean_csgs_b();

  // call elements to calculate system matrix and rhs and assemble
  assemble_mat_and_rhs();

  // apply Dirichlet boundary conditions
  apply_dirichlet_to_system();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::set_velocity_field()
{
  // safety check
  if (NdsVel() >= discret_->NumDofSets()) FOUR_C_THROW("Too few dofsets on scatra discretization!");

  // initialize velocity vectors
  Teuchos::RCP<Epetra_Vector> convel =
      Core::LinAlg::CreateVector(*discret_->dof_row_map(NdsVel()), true);
  Teuchos::RCP<Epetra_Vector> vel =
      Core::LinAlg::CreateVector(*discret_->dof_row_map(NdsVel()), true);

  switch (velocity_field_type_)
  {
    case Inpar::ScaTra::velocity_zero:
    {
      // no action needed in case for zero velocity field
      break;
    }

    case Inpar::ScaTra::velocity_function:
    {
      int err(0);
      const int velfuncno = params_->get<int>("VELFUNCNO");

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        Core::Nodes::Node* lnode = discret_->lRowNode(lnodeid);

        // get dofs associated with current node
        std::vector<int> nodedofs = discret_->Dof(NdsVel(), lnode);

        for (int index = 0; index < nsd_; ++index)
        {
          double value = problem_->FunctionById<Core::UTILS::FunctionOfSpaceTime>(velfuncno - 1)
                             .Evaluate(lnode->X().data(), time_, index);

          // get global and local dof IDs
          const int gid = nodedofs[index];
          const int lid = convel->Map().LID(gid);

          if (lid < 0) FOUR_C_THROW("Local ID not found in map for given global ID!");
          err = convel->ReplaceMyValue(lid, 0, value);
          if (err != 0) FOUR_C_THROW("error while inserting a value into convel");
          err = vel->ReplaceMyValue(lid, 0, value);
          if (err != 0) FOUR_C_THROW("error while inserting a value into vel");
        }
      }

      break;
    }

    default:
    {
      FOUR_C_THROW("Wrong SetVelocity() action for velocity field type %d!", velocity_field_type_);
      break;
    }
  }

  // provide scatra discretization with convective velocity
  discret_->set_state(NdsVel(), "convective velocity field", convel);

  // provide scatra discretization with velocity
  discret_->set_state(NdsVel(), "velocity field", vel);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::SetExternalForce()
{
  const auto input_params_external_force = params_->sublist("EXTERNAL FORCE");
  const int external_force_function_id = input_params_external_force.get<int>("FORCE_FUNCTION_ID");
  const int intrinsic_mobility_function_id =
      input_params_external_force.get<int>("INTRINSIC_MOBILITY_FUNCTION_ID");

  if (NdsVel() >= discret_->NumDofSets()) FOUR_C_THROW("Too few dofsets on scatra discretization!");

  // vector for the external force
  auto external_force = Core::LinAlg::CreateVector(*discret_->dof_row_map(NdsVel()), true);

  // vector for the intrinsic mobility
  auto intrinsic_mobility = Core::LinAlg::CreateVector(*discret_->dof_row_map(NdsVel()), true);

  // vector for the velocity due to the external force:
  // force_velocity = intrinsic_mobility * external_force
  auto force_velocity = Core::LinAlg::CreateVector(*discret_->dof_row_map(NdsVel()), true);

  for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
  {
    auto* const current_node = discret_->lRowNode(lnodeid);
    const auto nodedofs = discret_->Dof(NdsVel(), current_node);

    for (int spatial_dimension = 0; spatial_dimension < nsd_; ++spatial_dimension)
    {
      const double external_force_value =
          problem_->FunctionById<Core::UTILS::FunctionOfSpaceTime>(external_force_function_id - 1)
              .Evaluate(current_node->X().data(), time_, spatial_dimension);

      const double intrinsic_mobility_value =
          problem_
              ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(intrinsic_mobility_function_id - 1)
              .Evaluate(current_node->X().data(), time_, spatial_dimension);
      const double force_velocity_value = external_force_value * intrinsic_mobility_value;

      const int gid = nodedofs[spatial_dimension];
      const int lid = force_velocity->Map().LID(gid);

      if (lid < 0) FOUR_C_THROW("Local ID not found in map for given global ID!");
      const int error_force_velocity = force_velocity->ReplaceMyValue(lid, 0, force_velocity_value);
      if (error_force_velocity != 0)
        FOUR_C_THROW("Error while inserting a force_velocity_value into force_velocity.");

      const int error_external_force = external_force->ReplaceMyValue(lid, 0, external_force_value);
      if (error_external_force != 0)
        FOUR_C_THROW("Error while inserting a external_force_value into external_force.");

      const int error_intrinsic_mobility =
          intrinsic_mobility->ReplaceMyValue(lid, 0, intrinsic_mobility_value);
      if (error_intrinsic_mobility != 0)
        FOUR_C_THROW("Error while inserting a intrinsic_mobility_value into intrinsic_mobility.");
    }
  }

  discret_->set_state(NdsVel(), "external_force", external_force);
  discret_->set_state(NdsVel(), "intrinsic_mobility", intrinsic_mobility);
  discret_->set_state(NdsVel(), "force_velocity", force_velocity);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::set_wall_shear_stresses(Teuchos::RCP<const Epetra_Vector> wss)
{
  if (wss == Teuchos::null) FOUR_C_THROW("WSS state is Teuchos::null");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // We rely on the fact, that the nodal distribution of both fields is the same.
  // Although Scatra discretization was constructed as a clone of the fluid or
  // structure mesh, respectively, at the beginning, the nodal distribution may
  // have changed meanwhile (e.g., due to periodic boundary conditions applied only
  // to the fluid field)!
  // We have to be sure that everything is still matching.
  if (not wss->Map().SameAs(*discret_->dof_row_map(NdsWallShearStress())))
    FOUR_C_THROW("Maps are NOT identical. Emergency!");
#endif

  discret_->set_state(NdsWallShearStress(), "WallShearStress", wss);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::set_old_part_of_righthandside()
{
  // compute history values associated with meshtying strategy
  strategy_->SetOldPartOfRHS();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::SetPressureField(Teuchos::RCP<const Epetra_Vector> pressure)
{
  if (pressure == Teuchos::null) FOUR_C_THROW("Pressure state is Teuchos::null");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // We rely on the fact, that the nodal distribution of both fields is the same.
  // Although Scatra discretization was constructed as a clone of the fluid or
  // structure mesh, respectively, at the beginning, the nodal distribution may
  // have changed meanwhile (e.g., due to periodic boundary conditions applied only
  // to the fluid field)!
  // We have to be sure that everything is still matching.
  if (not pressure->Map().SameAs(*discret_->dof_row_map(NdsPressure())))
    FOUR_C_THROW("Maps are NOT identical. Emergency!");
#endif

  discret_->set_state(NdsPressure(), "Pressure", pressure);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::set_membrane_concentration(
    Teuchos::RCP<const Epetra_Vector> MembraneConc)
{
  if (MembraneConc == Teuchos::null) FOUR_C_THROW("MeanConc state is Teuchos::null");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // We rely on the fact, that the nodal distribution of both fields is the same.
  // Although Scatra discretization was constructed as a clone of the fluid or
  // structure mesh, respectively, at the beginning, the nodal distribution may
  // have changed meanwhile (e.g., due to periodic boundary conditions applied only
  // to the fluid field)!
  // We have to be sure that everything is still matching.
  if (not MembraneConc->Map().SameAs(*discret_->dof_row_map(0)))
    FOUR_C_THROW("Maps are NOT identical. Emergency!");
#endif

  // Note: we can not simply write this into the secondary discretisation here
  // since it is a variable of the primary dofset and is hence cleared
  // in between
  membrane_conc_ = MembraneConc;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::set_mean_concentration(Teuchos::RCP<const Epetra_Vector> MeanConc)
{
  if (MeanConc == Teuchos::null) FOUR_C_THROW("MeanConc state is Teuchos::null");

#ifdef FOUR_C_ENABLE_ASSERTIONS
  // We rely on the fact, that the nodal distribution of both fields is the same.
  // Although Scatra discretization was constructed as a clone of the fluid or
  // structure mesh, respectively, at the beginning, the nodal distribution may
  // have changed meanwhile (e.g., due to periodic boundary conditions applied only
  // to the fluid field)!
  // We have to be sure that everything is still matching.
  if (not MeanConc->Map().SameAs(*discret_->dof_row_map(0)))
    FOUR_C_THROW("Maps are NOT identical. Emergency!");
#endif

  // Note: we can not simply write this into the secondary discretisation here
  // since it is a variable of the primary dofset and is hence cleared
  // in between
  mean_conc_ = MeanConc;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::set_velocity_field(Teuchos::RCP<const Epetra_Vector> convvel,
    Teuchos::RCP<const Epetra_Vector> acc, Teuchos::RCP<const Epetra_Vector> vel,
    Teuchos::RCP<const Epetra_Vector> fsvel, const bool setpressure)
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA: set convective velocity field");

  //---------------------------------------------------------------------------
  // preliminaries
  //---------------------------------------------------------------------------
  if (convvel == Teuchos::null) FOUR_C_THROW("Velocity state is Teuchos::null");

  if (velocity_field_type_ != Inpar::ScaTra::velocity_Navier_Stokes)
    FOUR_C_THROW(
        "Wrong set_velocity_field() called for velocity field type %d!", velocity_field_type_);

  if (NdsVel() >= discret_->NumDofSets()) FOUR_C_THROW("Too few dofsets on scatra discretization!");

  // boolean indicating whether fine-scale velocity vector exists
  // -> if yes, multifractal subgrid-scale modeling is applied
  bool fsvelswitch = (fsvel != Teuchos::null);

  // some thing went wrong if we want to use multifractal subgrid-scale modeling
  // and have not got the fine-scale velocity
  if (step_ >= 1 and
      (turbmodel_ == Inpar::FLUID::multifractal_subgrid_scales or
          fssgd_ == Inpar::ScaTra::fssugrdiff_smagorinsky_small) and
      not fsvelswitch)
    FOUR_C_THROW("Fine-scale velocity expected for multifractal subgrid-scale modeling!");
  // as fsvelswitch is also true for smagorinsky_all, we have to reset fsvelswitch
  // as the corresponding vector, which is not necessary, is not provided in scatra
  if (fssgd_ == Inpar::ScaTra::fssugrdiff_smagorinsky_all and fsvelswitch) fsvelswitch = false;
  // as fsvelswitch is true in case of turned-off model in scalar field,
  // we have to ensure false
  if (turbmodel_ == Inpar::FLUID::no_model and fssgd_ == Inpar::ScaTra::fssugrdiff_no)
    fsvelswitch = false;

  // provide scatra discretization with convective velocity
  discret_->set_state(NdsVel(), "convective velocity field", convvel);

  // provide scatra discretization with velocity
  if (vel != Teuchos::null)
    discret_->set_state(NdsVel(), "velocity field", vel);
  else
  {
    // if velocity vector is not provided by the respective algorithm, we
    // assume that it equals the given convective velocity:
    discret_->set_state(NdsVel(), "velocity field", convvel);
  }

  // provide scatra discretization with acceleration field if required
  if (acc != Teuchos::null) discret_->set_state(NdsVel(), "acceleration field", acc);

  // provide scatra discretization with fine-scale convective velocity if required
  if (fsvelswitch) discret_->set_state(NdsVel(), "fine-scale velocity field", fsvel);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::TimeLoop()
{
  // safety checks
  check_is_init();
  check_is_setup();

  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:  + time loop");

  // prepare time loop
  prepare_time_loop();

  while (NotFinished())
  {
    // -------------------------------------------------------------------
    // prepare time step
    // -------------------------------------------------------------------
    prepare_time_step();

    // -------------------------------------------------------------------
    //                  solve nonlinear / linear equation
    // -------------------------------------------------------------------
    // store time before calling nonlinear solver
    double time = Teuchos::Time::wallTime();

    PreSolve();
    Solve();
    PostSolve();

    // determine time spent by nonlinear solver and take maximum over all processors via
    // communication
    double mydtnonlinsolve(Teuchos::Time::wallTime() - time), dtnonlinsolve(0.);
    discret_->Comm().MaxAll(&mydtnonlinsolve, &dtnonlinsolve, 1);

    // output performance statistics associated with nonlinear solver into *.csv file if applicable
    if (Core::UTILS::IntegralValue<int>(*params_, "OUTPUTNONLINSOLVERSTATS"))
      output_nonlin_solver_stats(iternum_, dtnonlinsolve, Step(), discret_->Comm());

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    Update();

    // -------------------------------------------------------------------
    // evaluate error for problems with analytical solution
    // -------------------------------------------------------------------
    evaluate_error_compared_to_analytical_sol();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    check_and_write_output_and_restart();

  }  // while

  // print the results of time measurements
  Teuchos::TimeMonitor::summarize();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::Solve()
{
  // safety checks
  check_is_init();
  check_is_setup();

  // -----------------------------------------------------------------
  // intermediate solution step for homogeneous isotropic turbulence
  // -----------------------------------------------------------------
  if (solvtype_ == Inpar::ScaTra::solvertype_nonlinear) calc_intermediate_solution();

  // -----------------------------------------------------------------
  //                     solve (non-)linear equation
  // -----------------------------------------------------------------
  switch (solvtype_)
  {
    case Inpar::ScaTra::solvertype_linear_incremental:
    case Inpar::ScaTra::solvertype_linear_full:
    {
      linear_solve();
      break;
    }

    case Inpar::ScaTra::solvertype_nonlinear:
    {
      nonlinear_solve();
      break;
    }

    case Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro:
    case Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken:
    case Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit:
    case Inpar::ScaTra::solvertype_nonlinear_multiscale_microtomacro:
    {
      nonlinear_multi_scale_solve();
      break;
    }

    default:
    {
      FOUR_C_THROW("Unknown solver type!");
      break;
    }
  }
  // that's all
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::Update()
{
  // update quantities associated with meshtying strategy
  strategy_->Update();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::ApplyMeshMovement(Teuchos::RCP<const Epetra_Vector> dispnp)
{
  //---------------------------------------------------------------------------
  // only required in ALE case
  //---------------------------------------------------------------------------
  if (isale_)
  {
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA: apply mesh movement");

    // check existence of displacement vector
    if (dispnp == Teuchos::null) FOUR_C_THROW("Got null pointer for displacements!");

    // provide scatra discretization with displacement field
    discret_->set_state(NdsDisp(), "dispnp", dispnp);
  }  // if (isale_)
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline void ScaTra::ScaTraTimIntImpl::print_time_step_info()
{
  if (myrank_ == 0)
  {
    std::cout << std::endl
              << "TIME: " << std::setw(11) << std::setprecision(4) << time_ << "/" << maxtime_
              << "  DT = " << dta_ << "  " << MethodTitle() << std::setw(4) << "  STEP = " << step_
              << "/" << stepmax_ << std::endl;
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseMatrix> ScaTra::ScaTraTimIntImpl::SystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(sysmat_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase> ScaTra::ScaTraTimIntImpl::BlockSystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::check_and_write_output_and_restart()
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + output of solution");

  // write result and potentially flux data
  if (IsResultStep()) WriteResult();

  // add restart data
  if (IsRestartStep()) write_restart();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::WriteResult()
{
  // step number and time (only after that data output is possible)
  output_->NewStep(step_, time_);

  // write domain decomposition for visualization (only once at the first time step!)
  if (step_ == 0) output_->WriteElementData(true);

  // write state vectors
  output_state();

  // write output to Gmsh postprocessing files
  if (outputgmsh_) output_to_gmsh(step_, time_);

  // write flux vector field (only writing, calculation was done during Update() call)
  if (calcflux_domain_ != Inpar::ScaTra::flux_none or
      calcflux_boundary_ != Inpar::ScaTra::flux_none)
  {
    // for flux output of initial field (before first solve) do:
    // flux_domain_ and flux_boundary_ vectors are initialized when CalcFlux() is called
    if (step_ == 0 or
        (calcflux_domain_ != Inpar::ScaTra::flux_none and flux_domain_ == Teuchos::null) or
        (calcflux_boundary_ != Inpar::ScaTra::flux_none and flux_boundary_ == Teuchos::null))
      CalcFlux(true);

    if (calcflux_domain_ != Inpar::ScaTra::flux_none) output_flux(flux_domain_, "domain");
    if (calcflux_boundary_ != Inpar::ScaTra::flux_none) output_flux(flux_boundary_, "boundary");
  }

  // write mean values of scalar(s)
  output_total_and_mean_scalars();

  // write domain and boundary integrals, i.e., surface areas and volumes of specified nodesets
  output_domain_or_boundary_integrals("DomainIntegral");
  output_domain_or_boundary_integrals("BoundaryIntegral");

  // write integral values of reaction(s)
  OutputIntegrReac();

  // problem-specific outputs
  output_problem_specific();

  // biofilm growth
  if (scfldgrdisp_ != Teuchos::null)
  {
    output_->WriteVector("scfld_growth_displ", scfldgrdisp_);
  }

  // biofilm growth
  if (scstrgrdisp_ != Teuchos::null)
  {
    output_->WriteVector("scstr_growth_displ", scstrgrdisp_);
  }

  // generate output associated with meshtying strategy
  strategy_->Output();

  // generate output on micro scale if necessary
  if (macro_scale_)
  {
    // create parameter list for macro elements
    Teuchos::ParameterList eleparams;

    // set action
    Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
        "action", ScaTra::Action::micro_scale_output, eleparams);

    // loop over macro-scale elements
    discret_->Evaluate(
        eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  }

  if ((step_ != 0) and (output_state_matlab_))
  {
    std::ostringstream filename;
    filename << problem_->OutputControlFile()->FileName() << "-Result_Step" << step_ << ".m";
    Core::LinAlg::PrintVectorInMatlabFormat(filename.str(), *phinp_);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::SetInitialField(
    const Inpar::ScaTra::InitialField init, const int startfuncno)
{
  switch (init)
  {
    case Inpar::ScaTra::initfield_zero_field:
    {
      phin_->PutScalar(0.0);
      phinp_->PutScalar(0.0);
      break;
    }
    case Inpar::ScaTra::initfield_field_by_function:
    case Inpar::ScaTra::initfield_disturbed_field_by_function:
    {
      const Epetra_Map* dofrowmap = discret_->dof_row_map();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        Core::Nodes::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        int numdofs = static_cast<int>(nodedofset.size());
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);
          // evaluate component k of spatial function
          double initialval =
              problem_->FunctionById<Core::UTILS::FunctionOfSpaceTime>(startfuncno - 1)
                  .Evaluate(lnode->X().data(), time_, k);
          int err = phin_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }
      }

      // for NURBS discretizations we have to solve a least squares problem,
      // with high accuracy! (do nothing for Lagrangian polynomials)
      const Teuchos::ParameterList& scatradyn = problem_->scalar_transport_dynamic_params();
      const int lstsolver = scatradyn.get<int>("LINEAR_SOLVER");

      auto* nurbsdis = dynamic_cast<Discret::Nurbs::NurbsDiscretization*>(&(*discret_));
      if (nurbsdis != nullptr)
      {
        if (lstsolver == (-1))
        {
          FOUR_C_THROW(
              "no linear solver defined for least square NURBS problem. Please set LINEAR_SOLVER "
              "in SCALAR TRANSPORT DYNAMIC to a valid number! Note: this solver block is misused "
              "for the least square problem. Maybe one should add a separate parameter for this.");
        }

        Discret::Nurbs::apply_nurbs_initial_condition(*discret_, problem_->SolverParams(lstsolver),
            problem_->FunctionById<Core::UTILS::FunctionOfSpaceTime>(startfuncno - 1), phin_);
      }

      // initialize also the solution vector. These values are a pretty good guess for the
      // solution after the first time step (much better than starting with a zero vector)
      phinp_->Update(1.0, *phin_, 0.0);

      // add random perturbation for initial field of turbulent flows
      if (init == Inpar::ScaTra::initfield_disturbed_field_by_function)
      {
        int err = 0;

        // random noise is relative to difference of max-min values of initial profile
        double perc =
            extraparams_->sublist("TURBULENCE MODEL").get<double>("CHAN_AMPL_INIT_DIST", 0.1);

        // out to screen
        if (myrank_ == 0)
        {
          std::cout << "Disturbed initial scalar profile:   max. " << perc * 100
                    << "% random perturbation\n";
          std::cout << "\n\n";
        }

        // get overall max and min values and range between min and max
        double maxphi(0.0);
        double minphi(0.0);
        err = phinp_->MaxValue(&maxphi);
        if (err > 0) FOUR_C_THROW("Error during evaluation of maximum value.");
        err = phinp_->MinValue(&minphi);
        if (err > 0) FOUR_C_THROW("Error during evaluation of minimum value.");
        double range = abs(maxphi - minphi);

        // disturb initial field for all degrees of freedom
        for (int k = 0; k < phinp_->MyLength(); ++k)
        {
          double randomnumber = problem_->Random()->Uni();
          double noise = perc * range * randomnumber;
          err += phinp_->SumIntoMyValues(1, &noise, &k);
          err += phin_->SumIntoMyValues(1, &noise, &k);
          if (err != 0) FOUR_C_THROW("Error while disturbing initial field.");
        }
      }
      break;
    }
    case Inpar::ScaTra::initfield_field_by_condition:
    {
      // set initial field for ALL existing scatra fields in condition
      const std::string field = "ScaTra";

      // get initial field conditions
      std::vector<Core::Conditions::Condition*> initfieldconditions(0);
      discret_->GetCondition("Initfield", initfieldconditions);

      if (not initfieldconditions.size())
      {
        FOUR_C_THROW(
            "Tried to evaluate initial field by condition without a corresponding condition "
            "defined on the ScaTra discretization!");
      }
      if (scalarhandler_ == Teuchos::null) FOUR_C_THROW("scalarhandler_ is null pointer!");

      std::set<int> numdofpernode;
      for (auto& initfieldcondition : initfieldconditions)
      {
        const int condmaxnumdofpernode =
            scalarhandler_->num_dof_per_node_in_condition(*initfieldcondition, discret_);

        if (condmaxnumdofpernode != 0) numdofpernode.insert(condmaxnumdofpernode);
      }

      if (numdofpernode.empty()) FOUR_C_THROW("No DOFs defined on initial field condition!");

      const int maxnumdofpernode = *(numdofpernode.rbegin());

      std::vector<int> localdofs(maxnumdofpernode);
      for (int i = 0; i < maxnumdofpernode; i++)
      {
        localdofs[i] = i;
      }
      discret_->evaluate_initial_field(
          Global::Problem::Instance()->FunctionManager(), field, phin_, localdofs);

      // initialize also the solution vector. These values are a pretty good guess for the
      // solution after the first time step (much better than starting with a zero vector)
      phinp_->Update(1.0, *phin_, 0.0);

      break;
    }
    // discontinuous 0-1 field for progress variable in 1-D
    case Inpar::ScaTra::initfield_discontprogvar_1D:
    {
      const Epetra_Map* dofrowmap = discret_->dof_row_map();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        Core::Nodes::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        // get coordinate
        const double x = lnode->X()[0];

        int numdofs = static_cast<int>(nodedofset.size());
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);

          double initialval = 0.0;
          if (x > -1e-10) initialval = 1.0;

          int err = 0;
          err += phin_->ReplaceMyValues(1, &initialval, &doflid);
          // initialize also the solution vector. These values are a pretty good guess for the
          // solution after the first time step (much better than starting with a zero vector)
          err += phinp_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }
      }
      break;
    }
    // reconstructed initial profile for progress variable in x2-direction from
    // Lessani and Papalexandris (2006), also used in Moureau et al. (2007, 2009),
    // for two-dimensional flame-vortex interaction problem (x2=0-200)
    case Inpar::ScaTra::initfield_flame_vortex_interaction:
    {
      // locations separating region 1 from region 2 and region 2 from region 3
      const double loc12 = 98.5;
      const double loc23 = 103.0;

      // define parameters for region 1 (exponential function for curve fitting)
      const double beta1 = 1.65;
      const double delta1 = 1.0;
      const double trans1 = 100.0;

      // define parameters for region 2 (linear function for curve fitting)
      const double abs2 = 0.0879;
      const double fac2 = 0.139309333;
      const double trans2 = 98.5;

      // define parameters for region 3 (exponential function for curve fitting)
      const double beta3 = 3.506209;
      const double delta3 = 4.28875;
      const double trans3 = 103.0;

      const Epetra_Map* dofrowmap = discret_->dof_row_map();

      // define variable
      double initialval = 0.0;

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        Core::Nodes::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        // get x2-coordinate
        const double x2 = lnode->X()[1];

        int numdofs = static_cast<int>(nodedofset.size());
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);

          if (x2 < loc12 - 1e-10)
            initialval = (1.0 - (1.0 / beta1)) * exp((x2 - trans1) / delta1);
          else if (x2 > loc23 + 1e-10)
            initialval = 1.0 - (exp((1.0 - beta3) * (x2 - trans3) / delta3) / beta3);
          else
            initialval = fac2 * (x2 - trans2) + abs2;

          int err = 0;
          err += phin_->ReplaceMyValues(1, &initialval, &doflid);
          // initialize also the solution vector. These values are a pretty good guess for the
          // solution after the first time step (much better than starting with a zero vector)
          err += phinp_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }
      }
      break;
    }
    // initial mixture-fraction profile for Rayleigh-Taylor instability
    case Inpar::ScaTra::initfield_raytaymixfrac:
    {
      // define interface thickness, sinusoidal disturbance wave amplitude and pi
      const double delta = 0.002;
      const double alpha = 0.001;

      const Epetra_Map* dofrowmap = discret_->dof_row_map();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        Core::Nodes::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        // get x1- and x2-coordinate
        const double x1 = lnode->X()[0];
        const double x2 = lnode->X()[1];

        // interface disturbance
        double x2_int = 0.0;
        x2_int -= std::cos(4 * M_PI * x1);
        x2_int -= std::cos(14 * M_PI * x1);
        x2_int -= std::cos(23 * M_PI * x1);
        x2_int -= std::cos(28 * M_PI * x1);
        x2_int -= std::cos(33 * M_PI * x1);
        x2_int -= std::cos(42 * M_PI * x1);
        x2_int -= std::cos(51 * M_PI * x1);
        x2_int -= std::cos(59 * M_PI * x1);
        x2_int *= alpha;

        const double value = (x2_int - x2) / (2.0 * delta);

        // values required for tanh-distribution
        const double vp = exp(value);
        const double vm = exp(-value);

        int numdofs = static_cast<int>(nodedofset.size());
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);

          // compute tanh-distribution
          double initialval = 0.0;
          initialval = 0.5 * (1.0 + (vp - vm) / (vp + vm));

          int err = 0;
          err += phin_->ReplaceMyValues(1, &initialval, &doflid);
          // initialize also the solution vector. These values are a pretty good guess for the
          // solution after the first time step (much better than starting with a zero vector)
          err += phinp_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }
      }
      break;
    }
    // initial field for skew convection of L-shaped domain
    case Inpar::ScaTra::initfield_Lshapeddomain:
    {
      const Epetra_Map* dofrowmap = discret_->dof_row_map();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        Core::Nodes::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        // get x1- and x2-coordinate
        const double x1 = lnode->X()[0];
        const double x2 = lnode->X()[1];

        int numdofs = static_cast<int>(nodedofset.size());
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);

          // compute initial values 0.0 or 1.0 depending on geometrical location
          double initialval = 0.0;
          if ((x1 <= 0.25 and x2 <= 0.5) or (x1 <= 0.5 and x2 <= 0.25)) initialval = 1.0;

          int err = 0;
          err += phin_->ReplaceMyValues(1, &initialval, &doflid);
          // initialize also the solution vector. These values are a pretty good
          // guess for the solution after the first time step (much better than
          // starting with a zero vector)
          err += phinp_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }
      }
      break;
    }
    case Inpar::ScaTra::initfield_facing_flame_fronts:
    {
      const Epetra_Map* dofrowmap = discret_->dof_row_map();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        Core::Nodes::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        // get x1- and x2-coordinate
        const double x1 = lnode->X()[0];

        int numdofs = static_cast<int>(nodedofset.size());
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);
          // evaluate component k of spatial function

          double initialval;
          if (x1 < 0.0)
            initialval = -(x1 + 0.75);
          else
            initialval = x1 - 0.75;

          int err = 0;
          err += phin_->ReplaceMyValues(1, &initialval, &doflid);
          // initialize also the solution vector. These values are a pretty good guess for the
          // solution after the first time step (much better than starting with a zero vector)
          err += phinp_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }
      }
      break;
    }
    case Inpar::ScaTra::initfield_oracles_flame:
    {
      const Epetra_Map* dofrowmap = discret_->dof_row_map();

      const double eps = 0.00152;

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        Core::Nodes::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        // get x2-coordinate
        const double x2 = lnode->X()[1];

        int numdofs = static_cast<int>(nodedofset.size());
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);
          // evaluate component k of spatial function

          double initval = 0.0;

          // initial plane implementation for periodic spanwise boundary
          if (x2 >= 0.0)
            initval = (x2 - 0.0354) - eps;
          else
            initval = (-0.0354 - x2) - eps;
          int err = 0;
          err += phin_->ReplaceMyValues(1, &initval, &doflid);
          // initialize also the solution vector. These values are a pretty good guess for the
          // solution after the first time step (much better than starting with a zero vector)
          err += phinp_->ReplaceMyValues(1, &initval, &doflid);
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }
      }
      break;
    }
    case Inpar::ScaTra::initialfield_forced_hit_high_Sc:
    case Inpar::ScaTra::initialfield_forced_hit_low_Sc:
    {
      // initialize calculation of initial field based on fast Fourier transformation
      Teuchos::RCP<HomIsoTurbInitialScalarField> HitInitialScalarField =
          Teuchos::rcp(new ScaTra::HomIsoTurbInitialScalarField(*this, init));
      // calculate initial field
      HitInitialScalarField->calculate_initial_field();

      break;
    }
    default:
      FOUR_C_THROW("Unknown option for initial field: %d", init);
      break;
  }  // switch(init)
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc)
{
  // store incremental vector to be available for convergence check
  // if incremental vector is received from outside for coupled problem
  increment_->Update(1.0, *inc, 0.0);

  // update scalar values by adding increments
  phinp_->Update(1.0, *inc, 1.0);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::setup_krylov_space_projection(Core::Conditions::Condition* kspcond)
{
  // previously, scatra was able to define actual modes that formed a
  // nullspace. factors when assigned to scalars in the dat file. it could
  // take several scatra Krylov conditions each forming one mode like:
  // 3.0*c_1 + 2.0*c_2 = const
  // since this was never used, not even in ELCH-problems, and for the sake of
  // consistency, now only a singel scatra Krylov condition can be given, with
  // flags that triggers the scalars that are to be levelled by a projection
  // (like the pressure in a pure Dirichlet fluid problem).
  // furthermore, this is a step towards the ability to have projection on
  // more than one field.
  // to see the handling of different modes, the ability that is now lost, see
  // revision 17615.

  // confirm that mode flags are number of nodal dofs/scalars
  const int nummodes = kspcond->parameters().Get<int>("NUMMODES");
  if (nummodes != NumDofPerNode())
  {
    FOUR_C_THROW(
        "Expecting as many mode flags as nodal dofs in Krylov projection definition. Check "
        "dat-file!");
  }

  // get vector of mode flags as given in dat-file
  const auto* modeflags = &kspcond->parameters().Get<std::vector<int>>("ONOFF");

  // count actual active modes selected in dat-file
  std::vector<int> activemodeids;
  for (int rr = 0; rr < NumDofPerNode(); ++rr)
  {
    if (((*modeflags)[rr]) != 0)
    {
      activemodeids.push_back(rr);
    }
  }

  // get from dat-file definition how weights are to be computed
  const auto* weighttype = &kspcond->parameters().Get<std::string>("weight vector definition");

  // set flag for projection update true only if ALE and integral weights
  if (isale_ and (*weighttype == "integration")) updateprojection_ = true;

  // create the projector
  projector_ = Teuchos::rcp(
      new Core::LinAlg::KrylovProjector(activemodeids, weighttype, discret_->dof_row_map()));

  // update the projector
  update_krylov_space_projection();
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::update_krylov_space_projection()
{
  // loop over modes to create vectors within multi-vector
  // one could tailor a MapExtractor to extract necessary dofs, however this
  // seems to be an overkill, since normally only a single scalar is
  // projected. for only a second projected scalar it seems worthwhile. feel
  // free! :)

  // get Teuchos::RCP to kernel vector of projector
  Teuchos::RCP<Epetra_MultiVector> c = projector_->GetNonConstKernel();
  c->PutScalar(0.0);

  const std::string* weighttype = projector_->WeightType();
  // compute w_ as defined in dat-file
  if (*weighttype == "pointvalues")
  {
    FOUR_C_THROW("option pointvalues not implemented");
  }
  else if (*weighttype == "integration")
  {
    // get Teuchos::RCP to weight vector of projector
    Teuchos::RCP<Epetra_MultiVector> w = projector_->GetNonConstWeights();
    w->PutScalar(0.0);

    // get number of modes and their ids
    int nummodes = projector_->Nsdim();
    std::vector<int> modeids = projector_->Modes();

    // initialize dofid vector to -1
    Teuchos::RCP<Core::LinAlg::IntSerialDenseVector> dofids =
        Teuchos::rcp(new Core::LinAlg::IntSerialDenseVector(NumDofPerNode()));
    for (int rr = 0; rr < NumDofPerNode(); ++rr)
    {
      (*dofids)[rr] = -1;
    }

    Teuchos::ParameterList mode_params;

    // set parameters for elements that do not change over mode
    Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
        "action", ScaTra::Action::integrate_shape_functions, mode_params);

    // loop over all activemodes
    for (int imode = 0; imode < nummodes; ++imode)
    {
      // activate dof of current mode and add dofids to parameter list
      (*dofids)[modeids[imode]] = 1;
      mode_params.set("dofids", dofids);

      /*
      // evaluate KrylovSpaceProjection condition in order to get
      // integrated nodal basis functions w_
      // Note that in the case of definition integration based, the average
      // increment of the scalar quantity c will vanish in an integral sense
      //
      //                    /              /                      /
      //   /    \          |              |  /          \        |  /    \
      //  | w_*c | = c_i * | N_i(x) dx =  | | N_i(x)*c_i | dx =  | | c(x) | dx = 0
      //   \    /          |              |  \          /        |  \    /
      //                   /              /                      /
      */

      // get an Teuchos::RCP of the current column Epetra_Vector of the MultiVector
      Teuchos::RCP<Epetra_Vector> wi = Teuchos::rcp((*w)(imode), false);

      // compute integral of shape functions
      discret_->evaluate_condition(mode_params, Teuchos::null, Teuchos::null, wi, Teuchos::null,
          Teuchos::null, "KrylovSpaceProjection");

      // deactivate dof of current mode
      (*dofids)[modeids[imode]] = -1;

      // set the current kernel basis vector - not very nice
      for (int inode = 0; inode < discret_->NumMyRowNodes(); inode++)
      {
        Core::Nodes::Node* node = discret_->lRowNode(inode);
        std::vector<int> gdof = discret_->Dof(0, node);
        int err = c->ReplaceGlobalValue(gdof[modeids[imode]], imode, 1);
        if (err != 0) FOUR_C_THROW("error while inserting value into c");
      }

    }  // loop over modes

    // adapt weight vector according to meshtying case
    if (msht_ != Inpar::FLUID::no_meshtying)
    {
      FOUR_C_THROW(
          "Since meshtying for scatra is not tested under Krylov projection FOUR_C_THROW is "
          "introduced. "
          "Remove at own responsibility.");
      // meshtying_->adapt_krylov_projector(w);
    }

  }  // endif integration
  else
  {
    FOUR_C_THROW("unknown definition of weight vector w for restriction of Krylov space");
  }

  // adapt kernel vector according to meshtying case
  if (msht_ != Inpar::FLUID::no_meshtying)
  {
    FOUR_C_THROW(
        "Since meshtying for scatra is not tested under Krylov projection FOUR_C_THROW is "
        "introduced. "
        "Remove at own responsibility.");
    // meshtying_->adapt_krylov_projector(c);
  }

  // fillcomplete the projector to compute (w^T c)^(-1)
  projector_->fill_complete();
}

/*----------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::create_meshtying_strategy()
{
  if (msht_ != Inpar::FLUID::no_meshtying)  // fluid meshtying
  {
    strategy_ = Teuchos::rcp(new MeshtyingStrategyFluid(this));
  }
  else if (S2IMeshtying())  // scatra-scatra interface mesh tying
  {
    strategy_ = Teuchos::rcp(new MeshtyingStrategyS2I(this, *params_));
  }
  else if (heteroreaccoupling_)  // scatra-scatra interface coupling
  {
    strategy_ = Teuchos::rcp(new HeterogeneousReactionStrategy(this));
  }
  else if (arterycoupling_)
  {
    strategy_ = Teuchos::rcp(new MeshtyingStrategyArtery(this));
  }
  else  // standard case without meshtying
  {
    strategy_ = Teuchos::rcp(new MeshtyingStrategyStd(this));
  }
}  // ScaTraTimIntImpl::create_meshtying_strategy

/*----------------------------------------------------------------------------------------*
 *----------------------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::create_scalar_handler()
{
  scalarhandler_ = Teuchos::rcp(new ScalarHandler());
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::apply_dirichlet_to_system()
{
  // -------------------------------------------------------------------
  // Apply Dirichlet boundary conditions to system matrix
  // -------------------------------------------------------------------
  if (incremental_)
  {
    // blank residual DOFs which are on Dirichlet BC
    // We can do this because the values at the Dirichlet positions
    // are not used anyway.
    // We could avoid this though, if the dofrowmap would not include
    // the Dirichlet values as well. But it is expensive to avoid that.
    dbcmaps_->InsertCondVector(dbcmaps_->ExtractCondVector(zeros_), residual_);

    //--------- Apply Dirichlet boundary conditions to system of equations
    // residual values are supposed to be zero at Dirichlet boundaries
    increment_->PutScalar(0.0);

    {
      // time measurement: application of DBC to system
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + apply DBC to system");

      Core::LinAlg::apply_dirichlet_to_system(
          *sysmat_, *increment_, *residual_, *zeros_, *(dbcmaps_->CondMap()));
    }
  }
  else
  {
    // time measurement: application of DBC to system
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + apply DBC to system");

    Core::LinAlg::apply_dirichlet_to_system(
        *sysmat_, *phinp_, *residual_, *phinp_, *(dbcmaps_->CondMap()));
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::apply_dirichlet_bc(
    const double time, Teuchos::RCP<Epetra_Vector> phinp, Teuchos::RCP<Epetra_Vector> phidt)
{
  // time measurement: apply Dirichlet conditions
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:      + apply dirich cond.");

  // Todo: what happens in  the case of generalized alpha
  // needed parameters
  Teuchos::ParameterList p;
  p.set("total time", time);  // actual time t_{n+1}
  p.set<const Core::UTILS::FunctionManager*>(
      "function_manager", &Global::Problem::Instance()->FunctionManager());
  const Core::ProblemType problem_type = Core::ProblemType::scatra;
  p.set<const Core::ProblemType*>("problem_type", &problem_type);

  // predicted Dirichlet values
  // \c  phinp then also holds prescribed new Dirichlet values
  discret_->evaluate_dirichlet(p, phinp, phidt, Teuchos::null, Teuchos::null, dbcmaps_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::scaling_and_neumann()
{
  // scaling to get true residual vector for all time integration schemes
  // in incremental case: boundary flux values can be computed from trueresidual
  if (incremental_) trueresidual_->Update(residual_scaling(), *residual_, 0.0);

  // add Neumann b.c. scaled with a factor due to time discretization
  add_neumann_to_residual();

  // add potential Neumann inflow or convective heat transfer boundary
  // conditions (simultaneous evaluation of both conditions not allowed!)
  if (neumanninflow_)
    compute_neumann_inflow(sysmat_, residual_);
  else if (convheatrans_)
    evaluate_convective_heat_transfer(sysmat_, residual_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::apply_neumann_bc(const Teuchos::RCP<Epetra_Vector>& neumann_loads)
{
  // prepare load vector
  neumann_loads->PutScalar(0.0);

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
      "action", ScaTra::BoundaryAction::calc_Neumann, condparams);

  condparams.set<const Core::UTILS::FunctionManager*>(
      "function_manager", &Global::Problem::Instance()->FunctionManager());

  // specific parameters
  add_problem_specific_parameters_and_vectors(condparams);

  // set time for evaluation of point Neumann conditions as parameter depending on time integration
  // scheme line/surface/volume Neumann conditions use the time stored in the time parameter class
  set_time_for_neumann_evaluation(condparams);

  // evaluate Neumann boundary conditions at time t_{n+alpha_F} (generalized alpha) or time t_{n+1}
  // (otherwise)
  discret_->evaluate_neumann(condparams, *neumann_loads);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::evaluate_solution_depending_conditions(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix, Teuchos::RCP<Epetra_Vector> rhs)
{
  // evaluate Robin type boundary condition
  evaluate_robin_boundary_conditions(systemmatrix, rhs);

  // evaluate meshtying
  // this needs to be done as final step for consistency
  strategy_->EvaluateMeshtying();

  //----------------------------------------------------------------------
  // apply contact terms...
  // account for partitioning algorithm. The dofs in contact discretization must be frozen
  // before calling this function
  //----------------------------------------------------------------------
  if (contact_strategy_nitsche_ != Teuchos::null)
  {
    const auto fint_scatra =
        contact_strategy_nitsche_->GetRhsBlockPtr(CONTACT::VecBlockType::scatra);
    if (residual_->Update(1.0, *fint_scatra, 1.0)) FOUR_C_THROW("update failed");
  }

  // evaluate macro-micro coupling on micro scale in multi-scale scalar transport problems
  if (micro_scale_) evaluate_macro_micro_coupling();
  strategy_->evaluate_point_coupling();
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int ScaTra::ScaTraTimIntImpl::GetMaxDofSetNumber() const
{
  return std::max({nds_disp_, nds_growth_, nds_micro_, nds_pres_, nds_scatra_, nds_thermo_,
      nds_two_tensor_quantitiy_, nds_vel_, nds_wss_});
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::evaluate_additional_solution_depending_models(
    Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix, Teuchos::RCP<Epetra_Vector> rhs)
{
  // evaluate solution depending additional models
  // this point is unequal nullptr only if a scatra
  // adapter has been constructed.
  if (additional_model_evaluator_ != nullptr)
    additional_model_evaluator_->evaluate_additional_solution_depending_models(systemmatrix, rhs);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::evaluate_robin_boundary_conditions(
    Teuchos::RCP<Core::LinAlg::SparseOperator> matrix, Teuchos::RCP<Epetra_Vector> rhs)
{
  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
      "action", ScaTra::BoundaryAction::calc_Robin, condparams);

  // add element parameters and set state vectors according to time-integration scheme
  add_time_integration_specific_vectors();

  // evaluate ElchBoundaryKinetics conditions at time t_{n+1} or t_{n+alpha_F}
  discret_->evaluate_condition(
      condparams, matrix, Teuchos::null, rhs, Teuchos::null, Teuchos::null, "TransportRobin");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::assemble_mat_and_rhs()
{
  // safety check
  check_is_init();
  check_is_setup();

  // time measurement: element calls
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + element calls");

  // get cpu time
  const double tcpuele = Teuchos::Time::wallTime();

  // zero out matrix entries
  sysmat_->Zero();

  // reset the residual vector
  residual_->PutScalar(0.0);

  // create parameter list for elements
  Teuchos::ParameterList eleparams;

  // action for elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::calc_mat_and_rhs, eleparams);

  // DO THIS AT VERY FIRST!!!
  // compute reconstructed diffusive fluxes for better consistency
  const auto consistency = Core::UTILS::IntegralValue<Inpar::ScaTra::Consistency>(
      params_->sublist("STABILIZATION"), "CONSISTENCY");
  if (consistency == Inpar::ScaTra::consistency_l2_projection_lumped)
  {
    // compute flux approximation and add it to the parameter list
    add_flux_approx_to_parameter_list(eleparams);
  }

  // prepare dynamic Smagorinsky model if required,
  // i.e. calculate turbulent Prandtl number
  if (timealgo_ != Inpar::ScaTra::timeint_stationary)
  {
    dynamic_computation_of_cs();
    dynamic_computation_of_cv();
  }
  // this parameter list is required here to get the element-based filtered constants
  eleparams.sublist("TURBULENCE MODEL") = extraparams_->sublist("TURBULENCE MODEL");

  // AVM3 separation for incremental solver: get fine-scale part of scalar
  if (incremental_ and step_ > 0 and
      (fssgd_ != Inpar::ScaTra::fssugrdiff_no or
          turbmodel_ == Inpar::FLUID::multifractal_subgrid_scales))
    av_m3_separation();

  // add state vectors according to time-integration scheme
  add_time_integration_specific_vectors();

  // set external volume force (required, e.g., for forced homogeneous isotropic turbulence)
  if (homisoturb_forcing_ != Teuchos::null) homisoturb_forcing_->UpdateForcing(step_);

  if (forcing_ != Teuchos::null) discret_->set_state("forcing", forcing_);

  // add problem specific time-integration parameters
  add_problem_specific_parameters_and_vectors(eleparams);

  // call loop over elements (with or without subgrid-diffusivity(-scaling) vector)
  if (fssgd_ != Inpar::ScaTra::fssugrdiff_no)
    discret_->Evaluate(eleparams, sysmat_, Teuchos::null, residual_, subgrdiff_, Teuchos::null);
  else
    discret_->Evaluate(eleparams, sysmat_, residual_);

  //----------------------------------------------------------------------
  // apply weak Dirichlet boundary conditions
  //----------------------------------------------------------------------
  {
    Teuchos::ParameterList mhdbcparams;

    // set action for elements
    Core::UTILS::AddEnumClassToParameterList<ScaTra::BoundaryAction>(
        "action", ScaTra::BoundaryAction::calc_weak_Dirichlet, mhdbcparams);

    add_time_integration_specific_vectors();

    // evaluate all mixed hybrid Dirichlet boundary conditions
    discret_->evaluate_condition(mhdbcparams, sysmat_, Teuchos::null, residual_, Teuchos::null,
        Teuchos::null, "LineWeakDirichlet");

    discret_->evaluate_condition(mhdbcparams, sysmat_, Teuchos::null, residual_, Teuchos::null,
        Teuchos::null, "SurfaceWeakDirichlet");
  }

  // AVM3 scaling for non-incremental solver: scaling of normalized AVM3-based
  // fine-scale subgrid-diffusivity matrix by subgrid diffusivity
  if (not incremental_ and fssgd_ != Inpar::ScaTra::fssugrdiff_no) av_m3_scaling(eleparams);

  // potential residual scaling and potential addition of Neumann terms
  scaling_and_neumann();  // TODO: do we have to call this function twice??

  // evaluate solution-depending additional models
  evaluate_additional_solution_depending_models(sysmat_, residual_);

  // evaluate solution-depending boundary and interface conditions
  evaluate_solution_depending_conditions(sysmat_, residual_);

  // finalize assembly of system matrix
  sysmat_->Complete();

  // end time measurement for element and take average over all processors via communication
  double mydtele = Teuchos::Time::wallTime() - tcpuele;
  discret_->Comm().MaxAll(&mydtele, &dtele_, 1);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::linear_solve()
{
  // safety checks
  check_is_init();
  check_is_setup();

  // output to screen
  print_time_step_info();

  // clear state vectors
  discret_->ClearState();

  // preparations for solve
  PrepareLinearSolve();

  // Solve system in incremental or non-incremental case
  if (incremental_)
  {
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + solver calls");

    // get cpu time
    const double tcpusolve = Teuchos::Time::wallTime();

    Core::LinAlg::SolverParams solver_params;

    strategy_->Solve(solver_, sysmat_, increment_, residual_, phinp_, 1, solver_params);

    // end time measurement for solver
    dtsolve_ = Teuchos::Time::wallTime() - tcpusolve;

    //------------------------------------------------ update solution vector
    UpdateIter(increment_);

    //--------------------------------------------- compute norm of increment
    double incnorm_L2(0.0);
    double scalnorm_L2(0.0);
    increment_->Norm2(&incnorm_L2);
    phinp_->Norm2(&scalnorm_L2);

    if (myrank_ == 0)
    {
      printf("+-------------------------------+-------------+\n");
      {
        if (scalnorm_L2 > 1e-10)
          printf("|  relative increment (L2 norm) | %10.3E  |", incnorm_L2 / scalnorm_L2);
        else  // prevent division by an almost zero value
          printf("|  absolute increment (L2 norm) | %10.3E  |\n", incnorm_L2);
      }
      printf(" (ts=%10.3E,te=%10.3E)\n", dtsolve_, dtele_);
      printf("+-------------------------------+-------------+\n");
    }
  }
  else
  {
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + solver calls");

    // get cpu time
    const double tcpusolve = Teuchos::Time::wallTime();

    Core::LinAlg::SolverParams solver_params;

    strategy_->Solve(solver_, sysmat_, phinp_, residual_, phinp_, 1, solver_params);

    // end time measurement for solver
    dtsolve_ = Teuchos::Time::wallTime() - tcpusolve;

    if (myrank_ == 0) printf("Solvertype linear_full (ts=%10.3E,te=%10.3E)\n", dtsolve_, dtele_);
  }

  // -------------------------------------------------------------------
  // compute values at intermediate time steps (only for gen.-alpha)
  // -------------------------------------------------------------------
  compute_intermediate_values();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::nonlinear_solve()
{
  // safety checks
  check_is_init();
  check_is_setup();

  // time measurement: nonlinear iteration
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:   + nonlin. iteration/lin. solve");

  // out to screen
  print_time_step_info();

  // special preparations for multifractal subgrid-scale model
  if (turbmodel_ == Inpar::FLUID::multifractal_subgrid_scales) recompute_mean_csgs_b();

  //------------------------------ turn adaptive solver tolerance on/off
  const double ittol = params_->sublist("NONLINEAR").get<double>("CONVTOL");
  const bool isadapttol =
      (Core::UTILS::IntegralValue<int>(params_->sublist("NONLINEAR"), "ADAPTCONV"));
  const double adaptolbetter = params_->sublist("NONLINEAR").get<double>("ADAPTCONV_BETTER");
  double actresidual(0.0);

  // prepare Newton-Raphson iteration
  iternum_ = 0;

  // perform explicit predictor step (-> better starting point for nonlinear solver)
  const bool explpredictor =
      (Core::UTILS::IntegralValue<int>(params_->sublist("NONLINEAR"), "EXPLPREDICT") == 1);
  if (explpredictor)
  {
    // explicit predictor + recovery of DBC values
    auto phinp_dirich = dbcmaps_->ExtractCondVector(phinp_);
    explicit_predictor();
    dbcmaps_->InsertCondVector(phinp_dirich, phinp_);
  }

  // start Newton-Raphson iteration
  while (true)
  {
    // clear states
    discret_->ClearState();

    iternum_++;

    // call elements to calculate system matrix and rhs and assemble
    assemble_mat_and_rhs();

    // perform finite difference check on time integrator level
    if (fdcheck_ == Inpar::ScaTra::fdcheck_global) fd_check();

    // project residual such that only part orthogonal to nullspace is considered
    if (projector_ != Teuchos::null) projector_->ApplyPT(*residual_);

    // Apply Dirichlet boundary conditions to system of equations
    // residual values are supposed to be zero at Dirichlet boundaries
    {
      // time measurement: application of DBC to system
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + apply DBC to system");

      Core::LinAlg::apply_dirichlet_to_system(
          *sysmat_, *increment_, *residual_, *zeros_, *(dbcmaps_->CondMap()));
    }

    // abort nonlinear iteration if desired
    if (strategy_->abort_nonlin_iter(*this, actresidual)) break;

    // initialize increment vector
    increment_->PutScalar(0.0);

    {
      // get cpu time
      const double tcpusolve = Teuchos::Time::wallTime();

      // time measurement: call linear solver
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + call linear solver");

      // do adaptive linear solver tolerance (not in first solve)
      Core::LinAlg::SolverParams solver_params;
      if (isadapttol && iternum_ > 1)
      {
        solver_params.nonlin_tolerance = ittol;
        solver_params.nonlin_residual = actresidual;
        solver_params.lin_tol_better = adaptolbetter;
      }

      // reprepare Krylov projection only if ale and projection required
      if (updateprojection_) update_krylov_space_projection();

      solver_params.projector = projector_;

      strategy_->Solve(solver_, sysmat_, increment_, residual_, phinp_, iternum_, solver_params);

      solver_->ResetTolerance();

      // end time measurement for solver and take average over all processors via communication
      double mydtsolve = Teuchos::Time::wallTime() - tcpusolve;
      discret_->Comm().MaxAll(&mydtsolve, &dtsolve_, 1);

      // output performance statistics associated with linear solver into text file if applicable
      if (Core::UTILS::IntegralValue<int>(*params_, "OUTPUTLINSOLVERSTATS"))
        output_lin_solver_stats(strategy_->Solver(), dtsolve_, Step(), iternum_,
            strategy_->dof_row_map().NumGlobalElements());
    }

    //------------------------------------------------ update solution vector
    phinp_->Update(1.0, *increment_, 1.0);

    //-------- update values at intermediate time steps (only for gen.-alpha)
    compute_intermediate_values();

    // compute values at the interior of the elements (required for hdg)
    compute_interior_values();

    compute_time_derivative();
  }  // nonlinear iteration

  // calculate mean concentration of micro discretization and set state to nds_micro_
  if (macro_scale_ and NdsMicro() != -1) calc_mean_micro_concentration();
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::nonlinear_multi_scale_solve()
{
  // reset number of outer iterations
  iternum_outer_ = 0;

  // initialize relaxed macro-scale state vector
  Teuchos::RCP<Epetra_Vector> phinp_relaxed = Teuchos::rcp(new Epetra_Vector(*phinp_));

  // begin outer iteration loop
  while (true)
  {
    // increment iteration number
    iternum_outer_++;

    // store current state vector on macro scale
    phinp_inc_->Update(1., *phinp_relaxed, 0.);

    // solve micro scale first and macro scale second
    if (solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro or
        solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken or
        solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit)
    {
      // backup macro-scale state vector
      const Teuchos::RCP<Epetra_Vector> phinp = phinp_;

      // replace macro-scale state vector by relaxed macro-scale state vector as input for micro
      // scale
      phinp_ = phinp_relaxed;

      // solve micro-scale problems
      nonlinear_micro_scale_solve();

      // undo state vector replacement
      phinp_ = phinp;

      // solve macro-scale problem
      nonlinear_solve();
    }

    // solve macro scale first and micro scale second
    else if (solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_microtomacro)
    {
      // solve macro-scale problem
      nonlinear_solve();

      // solve micro-scale problems
      nonlinear_micro_scale_solve();
    }

    // compute increment of macro-scale state vector
    phinp_inc_->Update(1., *phinp_, -1.);

    // convergence check
    if (strategy_->AbortOuterIter(*this)) break;

    if (solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken or
        solvtype_ == Inpar::ScaTra::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit)
    {
      // compute difference between current and previous increments of macro-scale state vector
      Epetra_Vector phinp_inc_diff(*phinp_inc_);
      phinp_inc_diff.Update(-1., *phinp_inc_old_, 1.);

      // perform Aitken relaxation
      perform_aitken_relaxation(*phinp_relaxed, phinp_inc_diff);

      // update increment of macro-scale state vector
      phinp_inc_old_->Update(1., *phinp_inc_, 0.);
    }

    else
      // no relaxation
      phinp_relaxed = phinp_;
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::nonlinear_micro_scale_solve()
{
  // initialize parameter list for evaluation of macro-scale elements
  Teuchos::ParameterList eleparams;

  // set action for macro-scale elements
  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::micro_scale_solve, eleparams);

  // clear state vectors
  discret_->ClearState();

  // set state vectors
  add_time_integration_specific_vectors();

  // evaluate macro-scale elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}

/*--------------------------------------------------------------------------*
 *--------------------------------------------------------------------------*/
std::string ScaTra::ScaTraTimIntImpl::map_tim_int_enum_to_string(
    const enum Inpar::ScaTra::TimeIntegrationScheme term)
{
  // length of return std::string is 14 due to usage in formated screen output
  switch (term)
  {
    case Inpar::ScaTra::timeint_one_step_theta:
      return "One-Step-Theta";
      break;
    case Inpar::ScaTra::timeint_bdf2:
      return "     BDF2     ";
      break;
    case Inpar::ScaTra::timeint_stationary:
      return "  Stationary  ";
      break;
    case Inpar::ScaTra::timeint_gen_alpha:
      return "  Gen. Alpha  ";
      break;
    default:
      FOUR_C_THROW("Cannot cope with name enum %d", term);
      return "";
      break;
  }

  return "";
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::output_state()
{
  // solution
  output_->WriteVector("phinp", phinp_);

  // convective velocity (written in case of coupled simulations since volmortar is now possible)
  if (velocity_field_type_ == Inpar::ScaTra::velocity_function or
      velocity_field_type_ == Inpar::ScaTra::velocity_Navier_Stokes)
  {
    Teuchos::RCP<const Epetra_Vector> convel =
        discret_->GetState(NdsVel(), "convective velocity field");
    if (convel == Teuchos::null) FOUR_C_THROW("Cannot get state vector convective velocity");

    Teuchos::RCP<Epetra_MultiVector> convel_multi =
        convert_dof_vector_to_componentwise_node_vector(convel, NdsVel());

    output_->WriteVector("convec_velocity", convel_multi, Core::IO::nodevector);
  }

  // displacement field
  if (isale_)
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = discret_->GetState(NdsDisp(), "dispnp");
    if (dispnp == Teuchos::null)
      FOUR_C_THROW("Cannot extract displacement field from discretization");

    Teuchos::RCP<Epetra_MultiVector> dispnp_multi =
        convert_dof_vector_to_componentwise_node_vector(dispnp, NdsDisp());

    output_->WriteVector("dispnp", dispnp_multi, Core::IO::nodevector);
  }

  if (NdsMicro() != -1)
  {
    // convert vector to multi vector
    auto micro_conc_multi = Teuchos::rcp(new Epetra_MultiVector(*discret_->NodeRowMap(), 1, true));

    for (int inode = 0; inode < discret_->NumMyRowNodes(); ++inode)
      (*micro_conc_multi)[0][inode] = (*phinp_micro_)[inode];

    output_->WriteVector("micro_conc", micro_conc_multi, Core::IO::nodevector);
  }

  if (has_external_force_)
  {
    Teuchos::RCP<const Epetra_Vector> external_force =
        discret_->GetState(nds_vel_, "external_force");
    Teuchos::RCP<Epetra_MultiVector> output_external_force =
        convert_dof_vector_to_componentwise_node_vector(external_force, NdsVel());
    output_->WriteVector("external_force", output_external_force, Core::IO::nodevector);

    Teuchos::RCP<const Epetra_Vector> mobility = discret_->GetState(nds_vel_, "intrinsic_mobility");
    Teuchos::RCP<Epetra_MultiVector> output_intrinsic_mobility =
        convert_dof_vector_to_componentwise_node_vector(mobility, NdsVel());
    output_->WriteVector("intrinsic_mobility", output_intrinsic_mobility, Core::IO::nodevector);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector>
ScaTra::ScaTraTimIntImpl::convert_dof_vector_to_componentwise_node_vector(
    const Teuchos::RCP<const Epetra_Vector>& dof_vector, const int nds) const
{
  Teuchos::RCP<Epetra_MultiVector> componentwise_node_vector =
      Teuchos::rcp(new Epetra_MultiVector(*discret_->NodeRowMap(), nsd_, true));
  for (int inode = 0; inode < discret_->NumMyRowNodes(); ++inode)
  {
    Core::Nodes::Node* node = discret_->lRowNode(inode);
    for (int idim = 0; idim < nsd_; ++idim)
      (*componentwise_node_vector)[idim][inode] =
          (*dof_vector)[dof_vector->Map().LID(discret_->Dof(nds, node, idim))];
  }
  return componentwise_node_vector;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
inline void ScaTra::ScaTraTimIntImpl::increment_time_and_step()
{
  step_ += 1;
  time_ += dta_;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::adapt_time_step_size()
{
  timestepadapted_ = false;

  // check flag for adaptive time stepping
  if (Core::UTILS::IntegralValue<bool>(*params_, "ADAPTIVE_TIMESTEPPING"))
  {
    // initialize time step size with original value
    double dt(params_->get<double>("TIMESTEP"));

    // reduce time step size if necessary
    compute_time_step_size(dt);

    // adapt time step size if necessary
    if (dt < dta_ or dt > dta_)
    {
      // print information about adaptation of time step size to screen
      if (myrank_ == 0)
      {
        std::cout << std::scientific << std::setprecision(2) << std::endl;
        std::cout << "ADAPTIVE TIME STEPPING:" << std::endl;
        std::cout << "Time step size is " << (dt < dta_ ? "decreased" : "increased") << " from "
                  << dta_ << " to " << dt << "!" << std::endl;
      }

      // adapt time step size
      set_dt(dt);

      // time step was adapted
      timestepadapted_ = true;
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::compute_time_step_size(double& dt)
{
  strategy_->compute_time_step_size(dt);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::add_dirich_cond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
  condmaps.push_back(maptoadd);
  condmaps.push_back(dbcmaps_->CondMap());
  Teuchos::RCP<Epetra_Map> condmerged = Core::LinAlg::MultiMapExtractor::MergeMaps(condmaps);
  *dbcmaps_ = Core::LinAlg::MapExtractor(*(discret_->dof_row_map()), condmerged);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::add_time_integration_specific_vectors(bool forcedincrementalsolver)
{
  // add global state vectors associated with meshtying strategy
  strategy_->add_time_integration_specific_vectors();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::remove_dirich_cond(const Teuchos::RCP<const Epetra_Map> maptoremove)
{
  std::vector<Teuchos::RCP<const Epetra_Map>> othermaps;
  othermaps.push_back(maptoremove);
  othermaps.push_back(dbcmaps_->OtherMap());
  Teuchos::RCP<Epetra_Map> othermerged = Core::LinAlg::MultiMapExtractor::MergeMaps(othermaps);
  *dbcmaps_ = Core::LinAlg::MapExtractor(*(discret_->dof_row_map()), othermerged, false);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ScaTra::ScaTraTimIntImpl::dof_row_map() { return dof_row_map(0); }

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> ScaTra::ScaTraTimIntImpl::dof_row_map(int nds)
{
  const Epetra_Map* dofrowmap = discret_->dof_row_map(nds);
  return Teuchos::rcp(dofrowmap, false);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int ScaTra::ScaTraTimIntImpl::NumScal() const
{
  if (scalarhandler_ == Teuchos::null) FOUR_C_THROW("scalar handler was not initialized!");
  return scalarhandler_->NumScal();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int ScaTra::ScaTraTimIntImpl::NumDofPerNode() const
{
  if (scalarhandler_ == Teuchos::null) FOUR_C_THROW("scalar handler was not initialized!");
  return scalarhandler_->NumDofPerNode();
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
int ScaTra::ScaTraTimIntImpl::num_dof_per_node_in_condition(
    const Core::Conditions::Condition& condition) const
{
  if (scalarhandler_ == Teuchos::null) FOUR_C_THROW("scalar handler was not initialized!");
  return scalarhandler_->num_dof_per_node_in_condition(condition, discret_);
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
const std::map<const int, std::vector<double>>& ScaTra::ScaTraTimIntImpl::TotalScalars() const
{
  if (outputscalarstrategy_ == Teuchos::null) FOUR_C_THROW("output strategy was not initialized!");

  return outputscalarstrategy_->TotalScalars();
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
const std::map<const int, std::vector<double>>& ScaTra::ScaTraTimIntImpl::MeanScalars() const
{
  if (outputscalarstrategy_ == Teuchos::null) FOUR_C_THROW("output strategy was not initialized!");

  return outputscalarstrategy_->MeanScalars();
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
const std::vector<double>& ScaTra::ScaTraTimIntImpl::DomainIntegrals() const
{
  if (outputdomainintegralstrategy_ == Teuchos::null)
    FOUR_C_THROW("output strategy for domain integration was not initialized!");

  return outputdomainintegralstrategy_->DomainIntegrals();
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
const std::vector<double>& ScaTra::ScaTraTimIntImpl::BoundaryIntegrals() const
{
  if (outputdomainintegralstrategy_ == Teuchos::null)
    FOUR_C_THROW("output strategy for domain integration was not initialized!");

  return outputdomainintegralstrategy_->BoundaryIntegrals();
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::evaluate_macro_micro_coupling()
{
  // extract multi-scale coupling conditions
  std::vector<Teuchos::RCP<Core::Conditions::Condition>> conditions;
  discret_->GetCondition("ScatraMultiScaleCoupling", conditions);

  // loop over conditions
  for (auto& condition : conditions)
  {
    // extract nodal cloud
    const std::vector<int>* const nodeids = condition->GetNodes();
    if (nodeids == nullptr)
      FOUR_C_THROW("Multi-scale coupling condition does not have nodal cloud!");

    // loop over all nodes in nodal cloud
    for (int inode : *nodeids)
    {
      // process row nodes only
      if (discret_->NodeRowMap()->MyGID(inode))
      {
        // extract node
        Core::Nodes::Node* node = discret_->gNode(inode);
        if (node == nullptr)
          FOUR_C_THROW(
              "Cannot extract node with global ID %d from micro-scale discretization!", inode);

        // safety check
        if (node->NumElement() != 1)
          FOUR_C_THROW("Number of 1D elements adjacent to the boundary node must be 1!");

        // compute domain integration factor
        constexpr double four_pi = 4.0 * M_PI;
        const double fac = Core::UTILS::IntegralValue<bool>(*params_, "SPHERICALCOORDS")
                               ? *node->X().data() * *node->X().data() * four_pi
                               : 1.0;

        // extract degrees of freedom from node
        const std::vector<int> dofs = discret_->Dof(0, node);

        // loop over all degrees of freedom
        for (int gid : dofs)
        {
          // extract global and local IDs of degree of freedom
          const int lid = discret_->dof_row_map()->LID(gid);
          if (lid < 0) FOUR_C_THROW("Cannot extract degree of freedom with global ID %d!", gid);

          // compute matrix and vector contributions according to kinetic model for current
          // macro-micro coupling condition
          const int kinetic_model = condition->parameters().Get<int>("kinetic model");

          switch (kinetic_model)
          {
            case Inpar::S2I::kinetics_constperm:
            {
              // access real vector of constant permeabilities
              const std::vector<double>* permeabilities =
                  condition->parameters().GetIf<std::vector<double>>("permeabilities");
              if (permeabilities == nullptr)
                FOUR_C_THROW("Cannot access vector of permeabilities for macro-micro coupling!");
              if (permeabilities->size() != (unsigned)NumScal())
                FOUR_C_THROW("Number of permeabilities does not match number of scalars!");

              // compute and store micro-scale coupling flux
              q_ = (*permeabilities)[0] * ((*phinp_)[lid] - phinp_macro_[0]);

              // compute and store derivative of micro-scale coupling flux w.r.t. macro-scale state
              // variable
              dq_dphi_[0] = -(*permeabilities)[0];

              // assemble contribution from macro-micro coupling into global residual vector
              (*residual_)[lid] -=
                  Discret::ELEMENTS::ScaTraEleParameterTimInt::Instance(discret_->Name())
                      ->TimeFacRhs() *
                  q_ * fac;

              // assemble contribution from macro-micro coupling into global system matrix
              sysmat_->Assemble(
                  Discret::ELEMENTS::ScaTraEleParameterTimInt::Instance(discret_->Name())
                          ->TimeFac() *
                      (*permeabilities)[0] * fac,
                  gid, gid);

              break;
            }

            case Inpar::S2I::kinetics_butlervolmer:
            case Inpar::S2I::kinetics_butlervolmerreduced:
            {
              // access material of electrode
              Teuchos::RCP<const Mat::Electrode> matelectrode =
                  Teuchos::rcp_dynamic_cast<const Mat::Electrode>(node->Elements()[0]->Material());
              if (matelectrode == Teuchos::null)
                FOUR_C_THROW("Invalid electrode material for multi-scale coupling!");

              // access input parameters associated with current condition
              const int nume = condition->parameters().Get<int>("e-");
              if (nume != 1)
              {
                FOUR_C_THROW(
                    "Invalid number of electrons involved in charge transfer at "
                    "electrode-electrolyte interface!");
              }
              const std::vector<int>* stoichiometries =
                  condition->parameters().GetIf<std::vector<int>>("stoichiometries");
              if (stoichiometries == nullptr)
              {
                FOUR_C_THROW(
                    "Cannot access vector of stoichiometric coefficients for multi-scale "
                    "coupling!");
              }
              if (stoichiometries->size() != 1)
                FOUR_C_THROW(
                    "Number of stoichiometric coefficients does not match number of scalars!");
              if ((*stoichiometries)[0] != -1) FOUR_C_THROW("Invalid stoichiometric coefficient!");
              const double faraday =
                  Global::Problem::Instance(0)->ELCHControlParams().get<double>("FARADAY_CONSTANT");
              const double gasconstant =
                  Global::Problem::Instance(0)->ELCHControlParams().get<double>("GAS_CONSTANT");
              const double frt =
                  faraday /
                  (gasconstant * (Global::Problem::Instance(0)->ELCHControlParams().get<double>(
                                     "TEMPERATURE")));
              const double alphaa =
                  condition->parameters().Get<double>("alpha_a");  // anodic transfer coefficient
              const double alphac =
                  condition->parameters().Get<double>("alpha_c");  // cathodic transfer coefficient
              const double kr = condition->parameters().Get<double>(
                  "k_r");  // rate constant of charge transfer reaction
              if (kr < 0.) FOUR_C_THROW("Charge transfer constant k_r is negative!");

              // extract saturation value of intercalated lithium concentration from electrode
              // material
              const double cmax = matelectrode->CMax();
              if (cmax < 1.e-12)
                FOUR_C_THROW(
                    "Saturation value c_max of intercalated lithium concentration is too small!");

              // extract electrode-side and electrolyte-side concentration values at multi-scale
              // coupling point
              const double conc_ed = (*phinp_)[lid];
              const double conc_el = phinp_macro_[0];

              // evaluate overall integration factors
              const double timefacfac =
                  Discret::ELEMENTS::ScaTraEleParameterTimInt::Instance(discret_->Name())
                      ->TimeFac() *
                  fac;
              const double timefacrhsfac =
                  Discret::ELEMENTS::ScaTraEleParameterTimInt::Instance(discret_->Name())
                      ->TimeFacRhs() *
                  fac;
              if (timefacfac < 0. or timefacrhsfac < 0.)
                FOUR_C_THROW("Integration factor is negative!");

              // no deformation available
              const double dummy_detF(1.0);

              // equilibrium electric potential difference and its derivative w.r.t. concentration
              // at electrode surface
              const double epd =
                  matelectrode->compute_open_circuit_potential(conc_ed, faraday, frt, dummy_detF);
              const double epdderiv =
                  matelectrode->compute_d_open_circuit_potential_d_concentration(
                      conc_ed, faraday, frt, dummy_detF);

              const double eta = phinp_macro_[2] - phinp_macro_[1] - epd;

              // Butler-Volmer exchange mass flux density
              const double j0 = condition->parameters().Get<int>("kinetic model") ==
                                        Inpar::S2I::kinetics_butlervolmerreduced
                                    ? kr
                                    : kr * std::pow(conc_el, alphaa) *
                                          std::pow(cmax - conc_ed, alphaa) *
                                          std::pow(conc_ed, alphac);

              // exponential Butler-Volmer terms
              const double expterm1 = std::exp(alphaa * frt * eta);
              const double expterm2 = std::exp(-alphac * frt * eta);
              const double expterm = expterm1 - expterm2;

              // core residual term associated with Butler-Volmer mass flux density
              q_ = j0 * expterm;

              const double dummyresistance(0.0);
              // define flux linearization terms
              double dj_dc_ed(0.0), dj_dc_el(0.0), dj_dpot_ed(0.0), dj_dpot_el(0.0);
              // calculate flux linearizations
              Discret::ELEMENTS::CalculateButlerVolmerElchLinearizations(kinetic_model, j0, frt,
                  epdderiv, alphaa, alphac, dummyresistance, expterm1, expterm2, kr, faraday,
                  conc_el, conc_ed, cmax, eta, dj_dc_ed, dj_dc_el, dj_dpot_ed, dj_dpot_el);

              dq_dphi_[0] = dj_dc_el;
              dq_dphi_[1] = dj_dpot_el;
              dq_dphi_[2] = dj_dpot_ed;

              // assemble contribution from macro-micro coupling into global residual vector
              (*residual_)[lid] -= timefacrhsfac * q_;

              // assemble contribution from macro-micro coupling into micro global system matrix
              sysmat_->Assemble(timefacfac * dj_dc_ed, gid, gid);

              break;
            }
            case Inpar::S2I::kinetics_nointerfaceflux:
              break;

            default:
            {
              FOUR_C_THROW("Kinetic model for macro-micro coupling not yet implemented!");
              break;
            }
          }
        }
      }
    }
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::check_is_init() const
{
  if (not is_init()) FOUR_C_THROW("ScaTraTimIntImpl is not initialized. Call Init() first.");
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::check_is_setup() const
{
  if (not is_setup()) FOUR_C_THROW("ScaTraTimIntImpl is not set up. Call Setup() first.");
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::setup_matrix_block_maps()
{
  if (matrixtype_ == Core::LinAlg::MatrixType::block_condition or
      matrixtype_ == Core::LinAlg::MatrixType::block_condition_dof)
  {
    // extract domain partitioning conditions from discretization
    std::vector<Teuchos::RCP<Core::Conditions::Condition>> partitioningconditions;
    discret_->GetCondition("ScatraPartitioning", partitioningconditions);

    // safety check
    if (partitioningconditions.empty())
    {
      FOUR_C_THROW(
          "For block preconditioning based on domain partitioning, at least one associated "
          "condition needs to be specified in the input file!");
    }

    // build maps associated with blocks of global system matrix
    std::vector<Teuchos::RCP<const Epetra_Map>> blockmaps;
    BuildBlockMaps(partitioningconditions, blockmaps);

    // initialize full map extractor associated with blocks of global system matrix
    blockmaps_ =
        Teuchos::rcp(new Core::LinAlg::MultiMapExtractor(*(discret_->dof_row_map()), blockmaps));
    // safety check
    blockmaps_->check_for_valid_map_extractor();
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::BuildBlockMaps(
    const std::vector<Teuchos::RCP<Core::Conditions::Condition>>& partitioningconditions,
    std::vector<Teuchos::RCP<const Epetra_Map>>& blockmaps) const
{
  if (matrixtype_ == Core::LinAlg::MatrixType::block_condition)
  {
    for (const auto& cond : partitioningconditions)
    {
      // all dofs that form one block map
      std::vector<int> dofs;

      for (int nodegid : *cond->GetNodes())
      {
        if (discret_->HaveGlobalNode(nodegid) and
            discret_->gNode(nodegid)->Owner() == discret_->Comm().MyPID())
        {
          const std::vector<int> nodedofs = discret_->Dof(0, discret_->gNode(nodegid));
          std::copy(nodedofs.begin(), nodedofs.end(), std::inserter(dofs, dofs.end()));
        }
      }
#ifdef FOUR_C_ENABLE_ASSERTIONS
      std::unordered_set<int> dof_set(dofs.begin(), dofs.end());
      FOUR_C_ASSERT(dof_set.size() == dofs.size(), "The dofs are not unique");
#endif

      blockmaps.emplace_back(Teuchos::rcp(
          new Epetra_Map(-1, static_cast<int>(dofs.size()), dofs.data(), 0, discret_->Comm())));
    }
  }
  else
    FOUR_C_THROW("Invalid type of global system matrix!");
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::post_setup_matrix_block_maps()
{
  // now build the null spaces
  build_block_null_spaces(Solver(), 0);

  // in case of an extended solver for scatra-scatra interface meshtying including interface growth
  // we need to equip it with the null space information generated above
  if (S2IMeshtying()) strategy_->equip_extended_solver_with_null_space_info();
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::build_block_null_spaces(
    Teuchos::RCP<Core::LinAlg::Solver> solver, int init_block_number) const
{
  // loop over blocks of global system matrix
  for (int iblock = init_block_number; iblock < BlockMaps()->NumMaps() + init_block_number;
       ++iblock)
  {
    // store number of current block as string, starting from 1
    std::ostringstream iblockstr;
    iblockstr << iblock + 1;

    // equip smoother for current matrix block with empty parameter sublists to trigger null space
    // computation
    Teuchos::ParameterList& blocksmootherparams =
        solver->Params().sublist("Inverse" + iblockstr.str());
    blocksmootherparams.sublist("Belos Parameters");
    blocksmootherparams.sublist("MueLu Parameters");

    // equip smoother for current matrix block with null space associated with all degrees of
    // freedom on discretization
    discret_->compute_null_space_if_necessary(blocksmootherparams);

    // reduce full null space to match degrees of freedom associated with current matrix block
    Core::LinearSolver::Parameters::FixNullSpace("Block " + iblockstr.str(),
        *discret_->dof_row_map(), *BlockMaps()->Map(iblock - init_block_number),
        blocksmootherparams);
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::setup_matrix_block_maps_and_meshtying()
{
  switch (MatrixType())
  {
    // case Core::LinAlg::MatrixType::undefined:
    case Core::LinAlg::MatrixType::sparse:
    {
      // only setup the meshtying in this case, as matrix has no block structure
      strategy_->setup_meshtying();

      break;
    }
    case Core::LinAlg::MatrixType::block_condition:
    case Core::LinAlg::MatrixType::block_condition_dof:
    {
      // safety check
      if (!Solver()->Params().isSublist("AMGnxn Parameters"))
        FOUR_C_THROW(
            "Global system matrix with block structure requires AMGnxn block preconditioner!");

      // setup the matrix block maps
      setup_matrix_block_maps();

      // setup the meshtying
      strategy_->setup_meshtying();

      // do some post setup matrix block map operations after the call to setup_meshtying, as they
      // rely on the fact that the interface maps have already been built
      post_setup_matrix_block_maps();

      break;
    }
    default:
    {
      FOUR_C_THROW("ScaTra Matrixtype %i not recognised", static_cast<int>(MatrixType()));
      break;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::SparseOperator> ScaTra::ScaTraTimIntImpl::init_system_matrix() const
{
  Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix(Teuchos::null);

  switch (matrixtype_)
  {
    case Core::LinAlg::MatrixType::sparse:
    {
      // initialize system matrix
      systemmatrix =
          Teuchos::rcp(new Core::LinAlg::SparseMatrix(*discret_->dof_row_map(), 27, false, true));
      break;
    }

    case Core::LinAlg::MatrixType::block_condition:
    case Core::LinAlg::MatrixType::block_condition_dof:
    {
      // initialize system matrix and associated strategy
      systemmatrix = Teuchos::rcp(
          new Core::LinAlg::BlockSparseMatrix<Core::LinAlg::DefaultBlockMatrixStrategy>(
              *BlockMaps(), *BlockMaps(), 81, false, true));

      break;
    }

    default:
    {
      FOUR_C_THROW(
          "Type of global system matrix for scatra-scatra interface coupling not recognized!");
      break;
    }
  }

  return systemmatrix;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::calc_mean_micro_concentration()
{
  phinp_micro_->PutScalar(0.0);

  if (NdsMicro() < 0) FOUR_C_THROW("must set number of dofset for micro scale concentrations");

  discret_->set_state("phinp", phinp_);

  Teuchos::ParameterList eleparams;

  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::calc_elch_elctrode_mean_concentration, eleparams);

  // evaluate nodal mean concentration of micro discretizations
  Core::FE::AssembleStrategy strategy(NdsMicro(), NdsMicro(), Teuchos::null, Teuchos::null,
      phinp_micro_, Teuchos::null, Teuchos::null);
  discret_->Evaluate(eleparams, strategy);

  // copy states from first dof of MAT_Electrode
  for (int ele_lid = 0; ele_lid < discret_->ElementRowMap()->NumMyElements(); ++ele_lid)
  {
    const int ele_gid = discret_->ElementRowMap()->GID(ele_lid);
    auto* ele = discret_->gElement(ele_gid);

    if (ele->Material()->MaterialType() != Core::Materials::m_electrode) continue;

    auto* nodes = ele->Nodes();

    for (int node_lid = 0; node_lid < ele->num_node(); ++node_lid)
    {
      // micro and macro dofs at this node
      auto* node = nodes[node_lid];
      int dof_macro = discret_->Dof(0, node)[0];
      int dof_micro = discret_->Dof(NdsMicro(), node)[0];

      const int dof_lid_micro = phinp_micro_->Map().LID(dof_micro);
      const int dof_lid_macro = phinp_->Map().LID(dof_macro);

      // only if owned by this proc
      if (dof_lid_micro != -1 and dof_lid_macro != -1)
      {
        const double macro_value = (*phinp_)[dof_lid_macro];
        // Sum, because afterwards it is divided by the number of adjacent nodes
        phinp_micro_->SumIntoMyValue(dof_lid_micro, 0, macro_value);
      }
    }
  }

  // divide nodal values by number of adjacent elements (due to assembly)
  const auto* node_row_map = discret_->NodeRowMap();
  for (int node_lid = 0; node_lid < node_row_map->NumMyElements(); ++node_lid)
  {
    const int node_gid = node_row_map->GID(node_lid);
    const auto* node = discret_->gNode(node_gid);
    std::vector<int> dofs = discret_->Dof(NdsMicro(), node);

    if (dofs.size() != 1) FOUR_C_THROW("Only one dof expected.");

    const int dof_gid = dofs[0];
    const int dof_lid = phinp_micro_->Map().LID(dof_gid);

    // only if this dof is part of the phinp_micro_ vector/map
    if (dof_lid != -1)
    {
      const double old_value = (*phinp_micro_)[dof_lid];
      const int num_elements = node->NumElement();
      const double new_value = old_value / static_cast<double>(num_elements);
      phinp_micro_->ReplaceMyValue(dof_lid, 0, new_value);
    }
  }

  // nodes with 3 dofs
  std::set<int> multiscale_nodes;
  // nodes with 2 dofs
  std::set<int> other_nodes;

  // loop over all element and search for nodes that are on elements with 2 dof on one side and 3
  // dofs at the other side
  for (int ele_lid = 0; ele_lid < discretization()->ElementRowMap()->NumMyElements(); ++ele_lid)
  {
    const int ele_gid = discretization()->ElementRowMap()->GID(ele_lid);
    auto* ele = discretization()->gElement(ele_gid);

    for (auto mat_id = 0; mat_id < ele->NumMaterial(); ++mat_id)
    {
      auto ele_mat = ele->Material(mat_id);
      auto material_type = ele_mat->MaterialType();

      if (material_type == Core::Materials::m_elchmat)
      {
        const auto* elchmat = static_cast<const Mat::ElchMat*>(ele_mat.get());

        const int num_dof_element = elchmat->NumDOF();

        const Core::Nodes::Node* const* nodes = ele->Nodes();
        for (int inode = 0; inode < ele->num_node(); ++inode)
        {
          if (num_dof_element == 3)
            Core::Communication::AddOwnedNodeGID(
                *discretization(), nodes[inode]->Id(), multiscale_nodes);
          else if (num_dof_element == 2)
            Core::Communication::AddOwnedNodeGID(
                *discretization(), nodes[inode]->Id(), other_nodes);
          else
            FOUR_C_THROW("Only 2 or 3 dofs per element supported");
        }
      }
    }
  }

  // find nodes that connect elements with 2 and 3 dofs ("hybrid nodes")
  std::vector<int> hybrid_nodes;
  for (int other_node : other_nodes)
  {
    for (int multiscale_node : multiscale_nodes)
    {
      if (other_node == multiscale_node)
      {
        hybrid_nodes.emplace_back(multiscale_node);
        break;
      }
    }
  }

  // get dofs from hybrid nodes
  std::vector<int> hybrid_dofs;
  for (int hybrid_node_gid : hybrid_nodes)
  {
    auto* hybrid_node = discretization()->gNode(hybrid_node_gid);
    auto dofs = discretization()->Dof(2, hybrid_node);
    for (int dof : dofs) hybrid_dofs.emplace_back(dof);
  }

  // correct values on hybrid dofs (value on node with 2 dofs is artificially set to 0.0)
  for (int hybrid_dof : hybrid_dofs)
  {
    const int lid = phinp_micro_->Map().LID(hybrid_dof);
    if (lid != -1)
    {
      const double value = (*phinp_micro_)[lid];
      const double corrected_value = 2.0 * value;
      phinp_micro_->ReplaceMyValue(lid, 0, corrected_value);
    }
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::add_problem_specific_parameters_and_vectors(
    Teuchos::ParameterList& params)
{
  if (micro_scale_) params.set<double>("rea_coeff", macro_micro_rea_coeff_);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::set_time_stepping_to_micro_scale()
{
  Teuchos::ParameterList eleparams;

  Core::UTILS::AddEnumClassToParameterList<ScaTra::Action>(
      "action", ScaTra::Action::micro_scale_set_time, eleparams);

  eleparams.set<double>("dt", dta_);
  eleparams.set<double>("time", time_);
  eleparams.set<int>("step", step_);

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::UTILS::ResultTest> ScaTra::ScaTraTimIntImpl::create_sca_tra_field_test()
{
  return Teuchos::rcp(new ScaTra::ScaTraResultTest(Teuchos::rcp(this, false)));
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::ScaTraTimIntImpl::TestResults()
{
  Global::Problem::Instance()->AddFieldTest(create_sca_tra_field_test());
  Global::Problem::Instance()->TestAll(discret_->Comm());
}

FOUR_C_NAMESPACE_CLOSE
