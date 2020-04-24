/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for convection-diffusion (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme

     o generalized-alpha time-integration scheme

     and stationary solver.

\level 1

\maintainer Anh-Tu Vuong

*/
/*----------------------------------------------------------------------*/
#include "scatra_timint_implicit.H"

#include "scatra_timint_heterogeneous_reaction_strategy.H"
#include "scatra_timint_meshtying_strategy_fluid.H"
#include "scatra_timint_meshtying_strategy_s2i.H"
#include "scatra_timint_meshtying_strategy_std.H"
#include "scatra_timint_meshtying_strategy_artery.H"
#include "turbulence_hit_initial_scalar_field.H"
#include "turbulence_hit_scalar_forcing.H"

#include "../drt_fluid/fluid_rotsym_periodicbc_utils.H"

#include "../drt_fluid_turbulence/dyn_vreman.H"

#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_elch.H"

#include "../drt_io/io.H"
#include "../drt_io/io_pstream.H"

#include "../drt_lib/drt_condition_selector.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_periodicbc.H"

#include "../drt_mat/electrode.H"
#include "../drt_mat/matlist.H"
#include "../drt_mat/matpar_bundle.H"
#include "../drt_mat/scatra_mat.H"

#include "../drt_nurbs_discret/drt_apply_nurbs_initial_condition.H"
#include "../drt_nurbs_discret/drt_nurbs_discret.H"

#include "../drt_scatra_ele/scatra_ele_action.H"
#include "../drt_scatra_ele/scatra_ele_parameter_timint.H"

#include "../linalg/linalg_krylov_projector.H"
#include "../linalg/linalg_solver.H"

// for the condition writer output
/*
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_lib/drt_utils.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io.H"
*/

/*
// for output of intermediate states in Newton loop
#include "../drt_io/io.H"
#include "../drt_io/io_control.H"
*/

//#define VISUALIZE_ELEDATA_GMSH
// only if VISUALIZE_ELEDATA_GMSH
//#include "../drt_io/io_gmsh.H"


namespace
{
  double distance_between_two_points(std::vector<double> point1, std::vector<double> point2)
  {
    if (point1.size() != point2.size())
      dserror("cannot calculate disctance between points of differing dimension");

    double distance = (point1[0] - point2[0]) * (point1[0] - point2[0]);
    for (unsigned int i = 1; i < point1.size(); ++i)
      distance += (point1[i] - point2[i]) * (point1[i] - point2[i]);

    return std::sqrt(distance);
  }
}  // namespace

/*==========================================================================*/
// Constructors and destructors and related methods
/*==========================================================================*/

/*----------------------------------------------------------------------*
 |  Constructor                                        (public) vg 05/07|
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntImpl::ScaTraTimIntImpl(
    Teuchos::RCP<DRT::Discretization> actdis,          //!< discretization
    Teuchos::RCP<LINALG::Solver> solver,               //!< linear solver
    Teuchos::RCP<Teuchos::ParameterList> params,       //!< parameter list
    Teuchos::RCP<Teuchos::ParameterList> extraparams,  //!< supplementary parameter list
    Teuchos::RCP<IO::DiscretizationWriter> output,     //!< output writer
    const int probnum                                  //!< global problem number
    )
    :  // call constructor for "nontrivial" objects
      problem_(DRT::Problem::Instance(probnum)),
      probnum_(probnum),
      solver_(solver),
      params_(params),
      extraparams_(extraparams),
      myrank_(actdis->Comm().MyPID()),
      splitter_(Teuchos::null),
      errfile_(extraparams->get<FILE*>("err file")),
      strategy_(Teuchos::null),
      additional_model_evaluator_(NULL),
      isale_(extraparams->get<bool>("isale")),
      equilibrationmethod_(
          Teuchos::getIntegralValue<INPAR::SCATRA::EquilibrationMethod>(*params, "EQUILIBRATION")),
      matrixtype_(Teuchos::getIntegralValue<INPAR::SCATRA::MatrixType>(*params, "MATRIXTYPE")),
      solvtype_(DRT::INPUT::IntegralValue<INPAR::SCATRA::SolverType>(*params, "SOLVERTYPE")),
      incremental_(true),
      fssgd_(DRT::INPUT::IntegralValue<INPAR::SCATRA::FSSUGRDIFF>(*params, "FSSUGRDIFF")),
      turbmodel_(INPAR::FLUID::no_model),
      s2icoupling_(actdis->GetCondition("S2ICoupling") != NULL),
      arterycoupling_(DRT::INPUT::IntegralValue<int>(
                          problem_->PoroMultiPhaseScatraDynamicParams(), "ARTERY_COUPLING") &&
                      actdis->Name() == "scatra"),
      heteroreaccoupling_(actdis->GetCondition("ScatraHeteroReactionSlave") != NULL),
      macro_scale_(problem_->Materials()->FirstIdByType(INPAR::MAT::m_scatra_multiscale) != -1 or
                   problem_->Materials()->FirstIdByType(INPAR::MAT::m_newman_multiscale) != -1),
      micro_scale_(probnum != 0),
      calcflux_domain_(
          DRT::INPUT::IntegralValue<INPAR::SCATRA::FluxType>(*params, "CALCFLUX_DOMAIN")),
      calcflux_domain_lumped_(DRT::INPUT::IntegralValue<bool>(*params, "CALCFLUX_DOMAIN_LUMPED")),
      calcflux_boundary_(
          DRT::INPUT::IntegralValue<INPAR::SCATRA::FluxType>(*params, "CALCFLUX_BOUNDARY")),
      calcflux_boundary_lumped_(
          DRT::INPUT::IntegralValue<bool>(*params, "CALCFLUX_BOUNDARY_LUMPED")),
      writefluxids_(Teuchos::rcp(new std::vector<int>)),
      flux_domain_(Teuchos::null),
      flux_boundary_(Teuchos::null),
      flux_boundary_maps_(Teuchos::null),
      sumnormfluxintegral_(Teuchos::null),
      lastfluxoutputstep_(-1),
      outputscalars_(
          DRT::INPUT::IntegralValue<INPAR::SCATRA::OutputScalarType>(*params, "OUTPUTSCALARS")),
      outputgmsh_(DRT::INPUT::IntegralValue<int>(*params, "OUTPUT_GMSH")),
      output_state_matlab_(DRT::INPUT::IntegralValue<int>(*params, "MATLAB_STATE_OUTPUT")),
      fdcheck_(DRT::INPUT::IntegralValue<INPAR::SCATRA::FDCheck>(*params, "FDCHECK")),
      fdcheckeps_(params->get<double>("FDCHECKEPS")),
      fdchecktol_(params->get<double>("FDCHECKTOL")),
      computeintegrals_(
          DRT::INPUT::IntegralValue<INPAR::SCATRA::ComputeIntegrals>(*params, "COMPUTEINTEGRALS")),
      calcerror_(DRT::INPUT::IntegralValue<INPAR::SCATRA::CalcError>(*params, "CALCERROR")),
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
          DRT::INPUT::IntegralValue<INPAR::SCATRA::TimeIntegrationScheme>(*params, "TIMEINTEGR")),
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
          DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(*params, "VELOCITYFIELD")),
      mean_conc_(Teuchos::null),
      membrane_conc_(Teuchos::null),
      nds_vel_(-1),
      nds_disp_(-1),
      nds_pres_(-1),
      nds_wss_(-1),
      densific_(0, 0.0),
      c0_(0, 0.0),
      discret_(actdis),
      output_(output),
      convform_(DRT::INPUT::IntegralValue<INPAR::SCATRA::ConvForm>(*params, "CONVFORM")),
      sysmat_(Teuchos::null),
      zeros_(Teuchos::null),
      dbcmaps_(Teuchos::null),
      neumann_loads_(Teuchos::null),
      normals_(Teuchos::null),
      residual_(Teuchos::null),
      trueresidual_(Teuchos::null),
      increment_(Teuchos::null),
      msht_(DRT::INPUT::IntegralValue<INPAR::FLUID::MeshTying>(*params, "MESHTYING")),
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
      turbinflow_(DRT::INPUT::IntegralValue<int>(
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
      neumanninflow_(DRT::INPUT::IntegralValue<int>(*params, "NEUMANNINFLOW")),
      convheatrans_(DRT::INPUT::IntegralValue<int>(*params, "CONV_HEAT_TRANS")),
      phinp_macro_(0, 0.),
      q_(0.0),
      dq_dphi_(0, 0.),
      // Initialization of Biofilm specific stuff
      scfldgrdisp_(Teuchos::null),
      scstrgrdisp_(Teuchos::null),
      outintegrreac_(DRT::INPUT::IntegralValue<int>(*params, "OUTINTEGRREAC")),
      skipinitder_(DRT::INPUT::IntegralValue<int>(*params, "SKIPINITDER")),
      issetup_(false),
      isinit_(false)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
}


/*------------------------------------------------------------------------*
 |  initialize time integration                               rauch 09/16 |
 *------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Init()
{
  SetIsSetup(false);

  // -------------------------------------------------------------------
  // safety check for spherical coordinates
  // -------------------------------------------------------------------
  if (DRT::INPUT::IntegralValue<bool>(*params_, "SPHERICALCOORDS") and nsd_ > 1)
    dserror("Spherical coordinates only available for 1D problems!");

  // -------------------------------------------------------------------
  // connect degrees of freedom for periodic boundary conditions
  // -------------------------------------------------------------------
  // note: pbcs have to be correctly set up before extended ghosting is applied
  Teuchos::RCP<PeriodicBoundaryConditions> pbc =
      Teuchos::rcp(new PeriodicBoundaryConditions(discret_, false));
  if (pbc->HasPBC() and not isinit_)
  {
    pbc->UpdateDofsForPeriodicBoundaryConditions();
  }

  // -------------------------------------------------------------------
  // determine whether linear incremental or nonlinear solver
  // -------------------------------------------------------------------
  switch (solvtype_)
  {
    case INPAR::SCATRA::solvertype_nonlinear:
    case INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro:
    case INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken:
    case INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit:
    case INPAR::SCATRA::solvertype_nonlinear_multiscale_microtomacro:
    case INPAR::SCATRA::solvertype_linear_incremental:
    {
      incremental_ = true;
    }
    break;
    case INPAR::SCATRA::solvertype_linear_full:
    {
      incremental_ = false;
    }
    break;
    default:
      dserror("Received illegal scatra solvertype enum.");
      break;
  }

  // -----------------------------------------------------------------------
  // determine number of degrees of freedom and transported scalars per node
  // -----------------------------------------------------------------------
  CreateScalarHandler();

  // -------------------------------------------------------------------
  // check compatibility of boundary conditions
  // -------------------------------------------------------------------
  if (neumanninflow_ and convheatrans_)
    dserror(
        "Neumann inflow and convective heat transfer boundary conditions must not appear "
        "simultaneously for the same problem!");

  // -----------------------------------------------------------------------------
  // initialize meshtying strategy (including standard strategy without meshtying)
  // -----------------------------------------------------------------------------
  // safety checks
  if (msht_ != INPAR::FLUID::no_meshtying and s2icoupling_)
    dserror(
        "Fluid-fluid meshtying in combination with scatra-scatra interface coupling is not "
        "implemented yet!");
  if (s2icoupling_ and !incremental_)
    dserror(
        "Scatra-scatra interface coupling only working for incremental solve so far!\n"
        "Set the parameter SOLVERTYPE in SCALAR TRANSPORT DYNAMIC section to 'nonlinear' or "
        "'linear_incremental'!");

  // create strategy
  CreateMeshtyingStrategy();

  // initialize strategy
  strategy_->InitMeshtying();

  // we have successfully initialized this class
  SetIsInit(true);
  return;
}  // ScaTraTimIntImpl::Init()


/*------------------------------------------------------------------------*
 |  initialize time integration                           rasthofer 09/13 |
 *------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Setup()
{
  // we have to call Init() first
  CheckIsInit();

  // compute Null Space
  ComputeNullSpaceIfNecessary();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices: local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->DofRowMap();

  // initialize the scalar handler
  if (scalarhandler_ == Teuchos::null)
    dserror("Make sure you construct the scalarhandler_ in initialization.");
  else
    scalarhandler_->Setup(this);

  // setup splitter (needed to solve initialization problems before SetupMeshtying())
  SetupSplitter();

  // setup the matrix block maps and the meshtying strategy
  SetupMatrixBlockMapsAndMeshtying();

  // -------------------------------------------------------------------
  // create empty system matrix (27 adjacent nodes as 'good' guess)
  // -------------------------------------------------------------------
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no and not incremental_)
    // do not save the graph if fine-scale subgrid diffusivity is used in non-incremental case (very
    // special case)
    sysmat_ = Teuchos::rcp(new LINALG::SparseMatrix(*(discret_->DofRowMap()), 27));
  else
    sysmat_ = strategy_->InitSystemMatrix();

  // -------------------------------------------------------------------
  // create vectors containing problem variables
  // -------------------------------------------------------------------
  // solutions at time n+1 and n
  phinp_ = LINALG::CreateVector(*dofrowmap, true);
  phin_ = LINALG::CreateVector(*dofrowmap, true);
  if (solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro or
      solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken or
      solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit or
      solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_microtomacro)
  {
    phinp_inc_ = LINALG::CreateVector(*dofrowmap, true);
    if (solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken or
        solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit)
    {
      phinp_inc_old_ = LINALG::CreateVector(*dofrowmap, true);
      if (solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken)
        omega_.resize(1, 1.);
      else
        omega_.resize(NumDofPerNode(), 1.);
    }
  }

  // temporal solution derivative at time n+1
  phidtnp_ = LINALG::CreateVector(*dofrowmap, true);
  // temporal solution derivative at time n
  phidtn_ = LINALG::CreateVector(*dofrowmap, true);

  // history vector (a linear combination of phinm, phin (BDF)
  // or phin, phidtn (One-Step-Theta, Generalized-alpha))
  hist_ = LINALG::CreateVector(*dofrowmap, true);

  // -------------------------------------------------------------------
  // create vectors associated to boundary conditions
  // -------------------------------------------------------------------
  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_ = LINALG::CreateVector(*dofrowmap, true);

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new LINALG::MapExtractor());
  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time", time_);
    discret_->EvaluateDirichlet(
        eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);
    zeros_->PutScalar(0.0);  // just in case of change
  }

  // -------------------------------------------------------------------
  // create vectors associated to solution process
  // -------------------------------------------------------------------
  // the vector containing body and surface forces
  neumann_loads_ = LINALG::CreateVector(*dofrowmap, true);

  // the residual vector --- more or less the rhs
  residual_ = LINALG::CreateVector(*dofrowmap, true);

  // residual vector containing the normal boundary fluxes
  trueresidual_ = LINALG::CreateVector(*dofrowmap, true);

  // incremental solution vector
  increment_ = LINALG::CreateVector(*dofrowmap, true);

  // subgrid-diffusivity(-scaling) vector
  // (used either for AVM3 approach or temperature equation
  //  with all-scale subgrid-diffusivity model)
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no) subgrdiff_ = LINALG::CreateVector(*dofrowmap, true);

  // -------------------------------------------------------------------
  // set parameters associated to potential statistical flux evaluations
  // -------------------------------------------------------------------
  // initialize vector for statistics (assume a maximum of 10 conditions)
  sumnormfluxintegral_ = Teuchos::rcp(new Epetra_SerialDenseVector(10));

  if (calcflux_domain_ != INPAR::SCATRA::flux_none or
      calcflux_boundary_ != INPAR::SCATRA::flux_none)
  {
    // safety check
    if (not scalarhandler_->EqualNumDof())
      dserror(
          "Flux output only implement for equal number of DOFs per node within ScaTra "
          "discretization!");

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
      IO::cout << "Flux output is performed for " << writefluxids_->size() << " scalars: ";
      for (unsigned int i = 0; i < writefluxids_->size(); i++)
      {
        const int id = (*writefluxids_)[i];
        IO::cout << id << " ";
        if ((id < 1) or (id > NumDofPerNode()))  // check validity of these numbers as well !
          dserror("Received illegal scalar id for flux output: %d", id);
      }
      IO::cout << IO::endl;
    }

    // initialize map extractor associated with boundary segments for flux calculation
    if (calcflux_boundary_ != INPAR::SCATRA::flux_none)
    {
      // extract conditions for boundary flux calculation
      std::vector<DRT::Condition*> conditions;
      discret_->GetCondition("ScaTraFluxCalc", conditions);

      // set up map extractor
      flux_boundary_maps_ = Teuchos::rcp(new LINALG::MultiMapExtractor());
      DRT::UTILS::MultiConditionSelector mcs;
      mcs.SetOverlapping(true);
      for (unsigned icond = 0; icond < conditions.size(); ++icond)
        mcs.AddSelector(Teuchos::rcp(new DRT::UTILS::ConditionSelector(
            *discret_, std::vector<DRT::Condition*>(1, conditions[icond]))));
      mcs.SetupExtractor(*discret_, *discret_->DofRowMap(), *flux_boundary_maps_);
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
  SetInitialField(DRT::INPUT::IntegralValue<INPAR::SCATRA::InitialField>(*params_, "INITIALFIELD"),
      params_->get<int>("INITFUNCNO"));

  // -------------------------------------------------------------------
  // preparations for natural convection
  // -------------------------------------------------------------------
  if (DRT::INPUT::IntegralValue<int>(*params_, "NATURAL_CONVECTION") == true)
  {
    // allocate global density vector and initialize
    densafnp_ = LINALG::CreateVector(*discret_->DofRowMap(), true);
    densafnp_->PutScalar(1.);
  }

  // -------------------------------------------------------------------
  // preparations for total and mean values of transported scalars
  // -------------------------------------------------------------------
  if (outputscalars_ != INPAR::SCATRA::outputscalars_none)
  {
    // input check
    if (outputscalars_ == INPAR::SCATRA::outputscalars_entiredomain)
    {
      std::vector<DRT::Condition*> conditions;
      // extract conditions for calculation of total and mean values of transported scalars
      discret_->GetCondition("TotalAndMeanScalar", conditions);
      // input check
      if (conditions.size())
        dserror(
            "Found 'DESIGN TOTAL AND MEAN SCALAR' condition on ScaTra discretization, but "
            "'OUTPUTSCALAR' \n"
            "in 'SCALAR TRANSPORT DYNAMIC' is set to 'entire domain'. Either switch on the output "
            "of mean and total scalars\n"
            "on conditions or remoove the 'DESIGN TOTAL AND MEAN SCALAR' condition from your input "
            "file!");
    }

    // build helper class for total and mean scalar output depending on input parameter
    switch (outputscalars_)
    {
      case INPAR::SCATRA::outputscalars_entiredomain:
        outputscalarstrategy_ = Teuchos::rcp(new OutputScalarsStrategyDomain);
        break;
      case INPAR::SCATRA::outputscalars_condition:
        outputscalarstrategy_ = Teuchos::rcp(new OutputScalarsStrategyCondition);
        break;
      case INPAR::SCATRA::outputscalars_entiredomain_condition:
        outputscalarstrategy_ = Teuchos::rcp(new OutputScalarsStrategyDomainAndCondition);
        break;
      default:
        dserror("Unknown option for output of total and mean scalars!");
        break;
    }

    // initialize scalar output strategy
    outputscalarstrategy_->Init(this);
  }
  else
  {
    // input check

    std::vector<DRT::Condition*> conditions;
    // extract conditions for calculation of total and mean values of transported scalars
    discret_->GetCondition("TotalAndMeanScalar", conditions);
    // input check
    if (conditions.size())
      dserror(
          "Found 'DESIGN TOTAL AND MEAN SCALAR' condition on ScaTra discretization, but "
          "'OUTPUTSCALAR' \n"
          "in 'SCALAR TRANSPORT DYNAMIC' is set to 'none'. Either switch on the output of mean and "
          "total scalars\n"
          "or remove the 'DESIGN TOTAL AND MEAN SCALAR' condition from your input file!");
  }

  // -------------------------------------------------------------------
  // preparations for domain integrals
  // -------------------------------------------------------------------
  if (computeintegrals_ != INPAR::SCATRA::computeintegrals_none)
  {
    // initialize domain integral output strategy
    outputdomainintegralstrategy_ = Teuchos::rcp(new OutputDomainIntegralStrategy);
    outputdomainintegralstrategy_->Init(this);
  }
  else
  {
    // input check
    std::vector<DRT::Condition*> conditions_boundary;
    // extract conditions for calculation of total and mean values of transported scalars
    discret_->GetCondition("BoundaryIntegral", conditions_boundary);
    std::vector<DRT::Condition*> conditions_domain;
    // extract conditions for calculation of total and mean values of transported scalars
    discret_->GetCondition("DomainIntegral", conditions_domain);
    // input check
    if (conditions_boundary.size() > 0 || conditions_domain.size() > 0)
      dserror(
          "Found 'DESIGN DOMAIN INTEGRAL SURF CONDITIONS' or 'DESIGN DOMAIN INTEGRAL VOL "
          "CONDITIONS' condition on ScaTra discretization, but COMPUTEINTEGRALS\n"
          "in 'SCALAR TRANSPORT DYNAMIC' is set to 'none'. Either switch on the output of domain "
          "integrals "
          "or remove the 'DESIGN DOMAIN INTEGRAL * CONDITIONS' condition from your input file!");
  }

  // -------------------------------------------------------------------
  // preparations for error evaluation
  // -------------------------------------------------------------------
  if (calcerror_ != INPAR::SCATRA::calcerror_no)
  {
    if (calcerror_ == INPAR::SCATRA::calcerror_bycondition)
    {
      std::vector<DRT::Condition*> relerrorconditions;
      discret_->GetCondition("ScatraRelError", relerrorconditions);
      const unsigned ncond = relerrorconditions.size();
      if (!ncond)
        dserror(
            "Calculation of relative error based on conditions desired, but no conditions "
            "specified!");
      relerrors_ =
          Teuchos::rcp(new std::vector<double>(2 * NumDofPerNode() * relerrorconditions.size()));
    }
    else if (calcerror_ == INPAR::SCATRA::calcerror_AnalyticSeries)
      relerrors_ = Teuchos::rcp(new std::vector<double>(2));  // TODO: Update two n species
    else
      relerrors_ = Teuchos::rcp(new std::vector<double>(2 * NumDofPerNode()));
  }

  // we have successfully set up this class
  SetIsSetup(true);
  return;
}  // ScaTraTimIntImpl::Setup()


/*----------------------------------------------------------------------*
 | perform setup of natural convection                       fang 08/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetupNatConv()
{
  // calculate the initial mean concentration value
  if (NumScal() < 1) dserror("Error since numscal = %d. Not allowed since < 1", NumScal());
  c0_.resize(NumScal());

  discret_->ClearState();
  discret_->SetState("phinp", phinp_);

  // set action for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", SCATRA::calc_total_and_mean_scalars);
  eleparams.set("inverting", false);

  // provide displacement field in case of ALE
  if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);

  // evaluate integrals of concentrations and domain
  Teuchos::RCP<Epetra_SerialDenseVector> scalars =
      Teuchos::rcp(new Epetra_SerialDenseVector(NumScal() + 1));
  discret_->EvaluateScalars(eleparams, scalars);
  discret_->ClearState();  // clean up

  // calculate mean concentrations
  const double domint = (*scalars)[NumScal()];
  if (std::abs(domint) < EPS15) dserror("Domain has zero volume!");
  for (int k = 0; k < NumScal(); ++k) c0_[k] = (*scalars)[k] / domint;

  // initialization of the densification coefficient vector
  densific_.resize(NumScal());
  DRT::Element* element = discret_->lRowElement(0);
  Teuchos::RCP<MAT::Material> mat = element->Material();

  if (mat->MaterialType() == INPAR::MAT::m_matlist or
      mat->MaterialType() == INPAR::MAT::m_matlist_reactions)
  {
    Teuchos::RCP<const MAT::MatList> actmat = Teuchos::rcp_static_cast<const MAT::MatList>(mat);

    for (int k = 0; k < NumScal(); ++k)
    {
      const int matid = actmat->MatID(k);
      Teuchos::RCP<const MAT::Material> singlemat = actmat->MaterialById(matid);

      if (singlemat->MaterialType() == INPAR::MAT::m_scatra)
      {
        Teuchos::RCP<const MAT::ScatraMat> actsinglemat =
            Teuchos::rcp_static_cast<const MAT::ScatraMat>(singlemat);

        densific_[k] = actsinglemat->Densification();

        if (densific_[k] < 0.0) dserror("received negative densification value");
      }
      else
        dserror("Material type is not allowed!");
    }
  }

  // for a single species calculation
  else if (mat->MaterialType() == INPAR::MAT::m_scatra)
  {
    Teuchos::RCP<const MAT::ScatraMat> actmat = Teuchos::rcp_static_cast<const MAT::ScatraMat>(mat);

    densific_[0] = actmat->Densification();

    if (densific_[0] < 0.0) dserror("received negative densification value");
    if (NumScal() > 1) dserror("Single species calculation but numscal = %d > 1", NumScal());
  }
  else
    dserror("Material type is not allowed!");

  return;
}  // ScaTraTimIntImpl::SetupNatConv


/*----------------------------------------------------------------------*
 | initialization of turbulence model                        ehrl 05/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::InitTurbulenceModel(
    const Epetra_Map* dofrowmap, const Epetra_Map* noderowmap)
{
  // get fluid turbulence sublist
  Teuchos::ParameterList* turbparams = &(extraparams_->sublist("TURBULENCE MODEL"));

  // parameters for statistical evaluation of normal fluxes
  samstart_ = turbparams->get<int>("SAMPLING_START");
  samstop_ = turbparams->get<int>("SAMPLING_STOP");
  dumperiod_ = turbparams->get<int>("DUMPING_PERIOD");
  if (dumperiod_ < 0) dserror("dumperiod_ is negative!");

  // -------------------------------------------------------------------
  // necessary only for AVM3 approach:
  // initialize subgrid-diffusivity matrix + respective output
  // -------------------------------------------------------------------
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no and
      DRT::INPUT::IntegralValue<int>(*turbparams, "TURBMODEL_LS"))
  {
    sysmat_sd_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap, 27));

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
      if (fssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_small and
          turbparams->get<std::string>("FSSUGRVISC") != "Smagorinsky_small")
        dserror("Same subgrid-viscosity approach expected!");
      if (fssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_all and
          turbparams->get<std::string>("FSSUGRVISC") != "Smagorinsky_all")
        dserror("Same subgrid-viscosity approach expected!");
    }
  }
  else
    fssgd_ = INPAR::SCATRA::fssugrdiff_no;  // in case of not "TURBMODEL_LS"

  // -------------------------------------------------------------------
  // get turbulence model and parameters
  // -------------------------------------------------------------------
  turbmodel_ = INPAR::FLUID::no_model;

  if (DRT::INPUT::IntegralValue<int>(*turbparams, "TURBMODEL_LS"))
  {
    // set turbulence model
    if (turbparams->get<std::string>("PHYSICAL_MODEL") == "Smagorinsky")
    {
      turbmodel_ = INPAR::FLUID::smagorinsky;

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
      turbmodel_ = INPAR::FLUID::dynamic_smagorinsky;
      // access to the dynamic Smagorinsky class will provided by the
      // scatra fluid couling algorithm
    }
    else if (turbparams->get<std::string>("PHYSICAL_MODEL") == "Multifractal_Subgrid_Scales")
    {
      turbmodel_ = INPAR::FLUID::multifractal_subgrid_scales;

      // initalize matrix used to build the scale separation operator
      sysmat_sd_ = Teuchos::rcp(new LINALG::SparseMatrix(*dofrowmap, 27));

      Teuchos::ParameterList* mfsparams = &(extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES"));
      if (mfsparams->get<std::string>("SCALE_SEPARATION") != "algebraic_multigrid_operator")
        dserror("Only scale separation by plain algebraic multigrid available in scatra!");

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
      turbmodel_ = INPAR::FLUID::dynamic_vreman;
    }

    // warning No. 1: if classical (all-scale) turbulence model other than
    // Smagorinsky or multifractal subrgid-scale modeling
    // is intended to be used
    if (turbparams->get<std::string>("PHYSICAL_MODEL") != "Smagorinsky" and
        turbparams->get<std::string>("PHYSICAL_MODEL") != "Dynamic_Smagorinsky" and
        turbparams->get<std::string>("PHYSICAL_MODEL") != "Multifractal_Subgrid_Scales" and
        turbparams->get<std::string>("PHYSICAL_MODEL") != "Dynamic_Vreman" and
        turbparams->get<std::string>("PHYSICAL_MODEL") != "no_model")
      dserror(
          "No classical (all-scale) turbulence model other than constant-coefficient Smagorinsky "
          "model and multifractal subrgid-scale modeling currently possible!");

    // warning No. 2: if classical (all-scale) turbulence model and fine-scale
    // subgrid-viscosity approach are intended to be used simultaneously
    if (turbmodel_ == INPAR::FLUID::smagorinsky and fssgd_ != INPAR::SCATRA::fssugrdiff_no)
      dserror(
          "No combination of classical turbulence model and fine-scale subgrid-diffusivity "
          "approach currently possible!");
  }

  if (turbmodel_ != INPAR::FLUID::no_model and NumScal() > 1)
    dserror("Turbulent passive scalar transport not supported for more than one scalar!");

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
      forcing_ = LINALG::CreateVector(*dofrowmap, true);
      forcing_->PutScalar(0.0);
    }
  }

  return;
}  // ScaTraTimIntImpl::InitTurbulenceModel()


/*----------------------------------------------------------------------*
 | create vectors for Krylov projection if necessary         ehrl 12/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::PrepareKrylovProjection()
{
  // sysmat might be singular (some modes are defined only up to a constant)
  // in this case, we need basis vectors for the nullspace/kernel

  // get condition "KrylovSpaceProjection" from discretization
  std::vector<DRT::Condition*> KSPCond;
  discret_->GetCondition("KrylovSpaceProjection", KSPCond);
  int numcond = KSPCond.size();
  int numscatra = 0;

  DRT::Condition* kspcond = NULL;
  // check if for scatra Krylov projection is required
  for (int icond = 0; icond < numcond; icond++)
  {
    const std::string* name = KSPCond[icond]->Get<std::string>("discretization");
    if (*name == "scatra")
    {
      numscatra++;
      kspcond = KSPCond[icond];
    }
  }

  // initialize variables for Krylov projection if necessary
  if (numscatra == 1)
  {
    SetupKrylovSpaceProjection(kspcond);
    if (myrank_ == 0)
      std::cout << "\nSetup of KrylovSpaceProjection in scatra field\n" << std::endl;
  }
  else if (numscatra == 0)
  {
    projector_ = Teuchos::null;
  }
  else
    dserror("Received more than one KrylovSpaceCondition for scatra field");

  return;
}


/*----------------------------------------------------------------------*
 | Destructor dtor                                   (public) gjb 04/08 |
 *----------------------------------------------------------------------*/
SCATRA::ScaTraTimIntImpl::~ScaTraTimIntImpl() { return; }


/*========================================================================*/
//! set element parameters
/*========================================================================*/

/*--------------------------------------------------------------------------------*
 | set all general parameters for element                              fang 10/14 |
 *--------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetElementGeneralParameters(bool calcinitialtimederivative) const
{
  Teuchos::ParameterList eleparams;

  // set action
  eleparams.set<int>("action", SCATRA::set_general_scatra_parameter);

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
  if (calcinitialtimederivative)  // deactivate finite difference check when calculating initial
                                  // time derivative
    eleparams.set<int>("fdcheck", INPAR::SCATRA::fdcheck_none);
  else
    eleparams.set<int>("fdcheck", fdcheck_);

  eleparams.set<double>("fdcheckeps", fdcheckeps_);
  eleparams.set<double>("fdchecktol", fdchecktol_);

  // flag for spherical coordinates
  eleparams.set<bool>(
      "sphericalcoords", DRT::INPUT::IntegralValue<bool>(*params_, "SPHERICALCOORDS"));

  // flag for truly partitioned multi-scale simulation
  eleparams.set<bool>("partitioned_multiscale",
      solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro or
          solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken or
          solvtype_ ==
              INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit or
          solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_microtomacro);

  // add parameters associated with meshtying strategy
  strategy_->SetElementGeneralParameters(eleparams);

  // additional problem-specific parameters for non-standard scalar transport problems
  // (electrochemistry etc.)
  SetElementSpecificScaTraParameters(eleparams);

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  return;
}


/*----------------------------------------------------------------------*
 | set turbulence parameters for element                rasthofer 11/13 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetElementTurbulenceParameters(bool calcinitialtimederivative) const
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", SCATRA::set_turbulence_scatra_parameter);

  eleparams.sublist("TURBULENCE MODEL") = extraparams_->sublist("TURBULENCE MODEL");
  if (calcinitialtimederivative)  // deactivate turbulence model when calculating initial time
                                  // derivative
    Teuchos::setStringToIntegralParameter<int>("PHYSICAL_MODEL", "no_model",
        "Classical LES approaches require an additional model for\nthe turbulent viscosity.",
        Teuchos::tuple<std::string>("no_model"),
        Teuchos::tuple<std::string>("If classical LES is our turbulence approach, this is a "
                                    "contradiction and should cause a dserror."),
        Teuchos::tuple<int>(INPAR::FLUID::no_model), &eleparams.sublist("TURBULENCE MODEL"));

  // set model-dependent parameters
  eleparams.sublist("SUBGRID VISCOSITY") = extraparams_->sublist("SUBGRID VISCOSITY");

  // set parameters for multifractal subgrid-scale modeling
  eleparams.sublist("MULTIFRACTAL SUBGRID SCALES") =
      extraparams_->sublist("MULTIFRACTAL SUBGRID SCALES");

  eleparams.set<bool>("turbulent inflow", turbinflow_);

  if (calcinitialtimederivative)
    eleparams.set<int>("fs subgrid diffusivity", INPAR::SCATRA::fssugrdiff_no);
  else
    eleparams.set<int>("fs subgrid diffusivity", fssgd_);

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  return;
}


/*==========================================================================*/
// general framework
/*==========================================================================*/

/*--- set, prepare, and predict --------------------------------------------*/

/*----------------------------------------------------------------------*
 | prepare time loop                                         fang 10/15 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::PrepareTimeLoop()
{
  // provide information about initial field (do not do for restarts!)
  if (step_ == 0)
  {
    // write out initial state
    Output();

    // compute error for problems with analytical solution (initial field!)
    EvaluateErrorComparedToAnalyticalSol();
  }

  return;
}  // SCATRA::ScaTraTimIntImpl::PrepareTimeLoop


/*----------------------------------------------------------------------*
 | setup the variables to do a new time step          (public)  vg 08/07|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::PrepareTimeStep()
{
  // time measurement: prepare time step
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + prepare time step");

  // -------------------------------------------------------------------
  //                       initialization
  // -------------------------------------------------------------------
  if (step_ == 0) PrepareFirstTimeStep();

  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  // adapt time step size if desired
  AdaptTimeStepSize();

  // note the order of the following three functions is important
  IncrementTimeAndStep();

  // -------------------------------------------------------------------
  // set part of the rhs vector belonging to the old timestep
  // -------------------------------------------------------------------
  SetOldPartOfRighthandside();
  // TODO (Thon): We do not really want to call SetElementTimeParameter() every time step.
  // But for now we just do it since "total time" has to be changed in the parameter class..
  SetElementTimeParameter();

  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  // TODO: Dirichlet auch im Fall von genalpha phinp
  // Neumann(n + alpha_f)
  ApplyDirichletBC(time_, phinp_, Teuchos::null);
  ApplyNeumannBC(neumann_loads_);

  // By definition: Applying DC on the slave side of an internal interface is not allowed
  //                since it leads to an over-constraint system
  // Therefore, nodes belonging to the slave side of an internal interface have to be excluded from
  // the DC. However, a velocity value (projected from the Dirichlet condition on the master side)
  // has to be assigned to the DOF's on the slave side in order to evaluate the system matrix
  // completely

  // Preparation for including DC on the master side in the condensation process
  strategy_->IncludeDirichletInCondensation();

  // -------------------------------------------------------------------
  //     update velocity field if given by function (it might depend on time)
  // -------------------------------------------------------------------
  if (velocity_field_type_ == INPAR::SCATRA::velocity_function) SetVelocityField(nds_vel_);

  // -------------------------------------------------------------------
  //           preparation of AVM3-based scale separation
  // -------------------------------------------------------------------
  if ((step_ == 1 or (turbinflow_ and step_ == numinflowsteps_ + 1)) and
      (fssgd_ != INPAR::SCATRA::fssugrdiff_no or
          turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales))
    AVM3Preparation();

  // -------------------------------------------------------------------
  // compute values at intermediate time steps (only for gen.-alpha)
  // -------------------------------------------------------------------
  ComputeIntermediateValues();

  // -------------------------------------------------------------------
  // prepare time step on micro scale if necessary
  // -------------------------------------------------------------------
  if (macro_scale_)
  {
    // create parameter list for macro-scale elements
    Teuchos::ParameterList eleparams;

    // set action
    eleparams.set<int>("action", SCATRA::micro_scale_prepare_time_step);

    // add state vectors
    AddTimeIntegrationSpecificVectors();

    // loop over macro-scale elements
    discret_->Evaluate(
        eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  }

  return;
}  // ScaTraTimIntImpl::PrepareTimeStep


/*------------------------------------------------------------------------------*
 | initialization procedure prior to evaluation of first time step   fang 09/15 |
 *------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::PrepareFirstTimeStep()
{
  // ApplyDirichletBC(time_,phin_,Teuchos::null);
  if (not skipinitder_)
  {
    if (nds_vel_ != -1)  // if some velocity field has been set
    {
      // TODO: Restructure enforcement of Dirichlet boundary conditions on phin_
      // A clean solution would incorporate ApplyDirichletBC(...) into CalcInitialTimeDerivative().
      // However, this would make a number of test cases fail. We should have a closer look at
      // this problem and fix it eventually.
      ApplyDirichletBC(time_, phin_, Teuchos::null);
      CalcInitialTimeDerivative();
    }

    // if initial velocity field has not been set here, the initial time derivative of phi will be
    // calculated wrongly for some time integration schemes
    else
      dserror("Initial velocity field has not been set!");

    // for safety; so we don't do this calculation twice
    // skipinitder_ = true;
  }

  return;
}  // SCATRA::ScaTraTimIntImpl::PrepareFirstTimeStep


/*----------------------------------------------------------------------*
 | preparations for solve                                (public) mr.x  |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::PrepareLinearSolve()
{
  // special preparations for multifractal subgrid-scale model
  if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales) RecomputeMeanCsgsB();

  // call elements to calculate system matrix and rhs and assemble
  AssembleMatAndRHS();

  // apply Dirichlet boundary conditions
  ApplyDirichletToSystem();
}


/*----------------------------------------------------------------------*
 | update the velocity field                                  gjb 04/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetVelocityField(const int nds)
{
  // safety check
  if (nds >= discret_->NumDofSets()) dserror("Too few dofsets on scatra discretization!");

  // initialize velocity vectors
  Teuchos::RCP<Epetra_Vector> convel = LINALG::CreateVector(*discret_->DofRowMap(nds), true);
  Teuchos::RCP<Epetra_Vector> vel = LINALG::CreateVector(*discret_->DofRowMap(nds), true);

  switch (velocity_field_type_)
  {
    case INPAR::SCATRA::velocity_zero:
    {
      // no action needed in case for zero velocity field
      break;
    }

    case INPAR::SCATRA::velocity_function:
    {
      int err(0);
      const int velfuncno = params_->get<int>("VELFUNCNO");

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(lnodeid);

        // get dofs associated with current node
        std::vector<int> nodedofs = discret_->Dof(nds, lnode);

        for (int index = 0; index < nsd_; ++index)
        {
          double value = problem_->Funct(velfuncno - 1).Evaluate(index, lnode->X(), time_);

          // get global and local dof IDs
          const int gid = nodedofs[index];
          const int lid = convel->Map().LID(gid);

          if (lid < 0) dserror("Local ID not found in map for given global ID!");
          err = convel->ReplaceMyValue(lid, 0, value);
          if (err != 0) dserror("error while inserting a value into convel");
          err = vel->ReplaceMyValue(lid, 0, value);
          if (err != 0) dserror("error while inserting a value into vel");
        }
      }

      break;
    }

    default:
    {
      dserror("Wrong SetVelocity() action for velocity field type %d!", velocity_field_type_);
      break;
    }
  }

  // store number of dof-set associated with velocity related dofs
  nds_vel_ = nds;

  // provide scatra discretization with convective velocity
  discret_->SetState(nds_vel_, "convective velocity field", convel);

  // provide scatra discretization with velocity
  discret_->SetState(nds_vel_, "velocity field", vel);

  return;

}  // ScaTraImplicitTimeInt::SetVelocityField


/*----------------------------------------------------------------------*
 | Set Wall Shear Stresses                                thon 11/15 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetWallShearStresses(
    Teuchos::RCP<const Epetra_Vector> wss, const int nds_wss)
{
  if (wss == Teuchos::null) dserror("WSS state is Teuchos::null");

#ifdef DEBUG
  // We rely on the fact, that the nodal distribution of both fields is the same.
  // Although Scatra discretization was constructed as a clone of the fluid or
  // structure mesh, respectively, at the beginning, the nodal distribution may
  // have changed meanwhile (e.g., due to periodic boundary conditions applied only
  // to the fluid field)!
  // We have to be sure that everything is still matching.
  if (not wss->Map().SameAs(*discret_->DofRowMap(nds_wss)))
    dserror("Maps are NOT identical. Emergency!");
#endif

  nds_wss_ = nds_wss;
  discret_->SetState(nds_wss, "WallShearStress", wss);
}


/*----------------------------------------------------------------------------------------------------------------------------------------------*
 | compute history vector, i.e., the history part of the right-hand side vector with all
 contributions from the previous time step   fang 01/17 |
 *----------------------------------------------------------------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetOldPartOfRighthandside()
{
  // compute history values associated with meshtying strategy
  strategy_->SetOldPartOfRHS();

  return;
}


/*----------------------------------------------------------------------*
 | Set Pressure Field                                        thon 11/15 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetPressureField(
    Teuchos::RCP<const Epetra_Vector> pressure, const int nds_pres)
{
  if (pressure == Teuchos::null) dserror("Pressure state is Teuchos::null");

#ifdef DEBUG
  // We rely on the fact, that the nodal distribution of both fields is the same.
  // Although Scatra discretization was constructed as a clone of the fluid or
  // structure mesh, respectively, at the beginning, the nodal distribution may
  // have changed meanwhile (e.g., due to periodic boundary conditions applied only
  // to the fluid field)!
  // We have to be sure that everything is still matching.
  if (not pressure->Map().SameAs(*discret_->DofRowMap(nds_pres)))
    dserror("Maps are NOT identical. Emergency!");
#endif

  nds_pres_ = nds_pres;
  discret_->SetState(nds_pres, "Pressure", pressure);
}

/*----------------------------------------------------------------------*
 | Set membrane concentration                                thon 08/16 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetMembraneConcentration(
    Teuchos::RCP<const Epetra_Vector> MembraneConc)
{
  if (MembraneConc == Teuchos::null) dserror("MeanConc state is Teuchos::null");

#ifdef DEBUG
  // We rely on the fact, that the nodal distribution of both fields is the same.
  // Although Scatra discretization was constructed as a clone of the fluid or
  // structure mesh, respectively, at the beginning, the nodal distribution may
  // have changed meanwhile (e.g., due to periodic boundary conditions applied only
  // to the fluid field)!
  // We have to be sure that everything is still matching.
  if (not MembraneConc->Map().SameAs(*discret_->DofRowMap(0)))
    dserror("Maps are NOT identical. Emergency!");
#endif

  // Note: we can not simply write this into the secondary discretisation here
  // since it is a variable of the primary dofset and is hence cleared
  // in between
  membrane_conc_ = MembraneConc;
  return;
}  // ScaTraTimIntImpl::SetMeanConcentration

/*----------------------------------------------------------------------*
 | Set mean concentrations                                hemmler 05/14 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetMeanConcentration(Teuchos::RCP<const Epetra_Vector> MeanConc)
{
  if (MeanConc == Teuchos::null) dserror("MeanConc state is Teuchos::null");

#ifdef DEBUG
  // We rely on the fact, that the nodal distribution of both fields is the same.
  // Although Scatra discretization was constructed as a clone of the fluid or
  // structure mesh, respectively, at the beginning, the nodal distribution may
  // have changed meanwhile (e.g., due to periodic boundary conditions applied only
  // to the fluid field)!
  // We have to be sure that everything is still matching.
  if (not MeanConc->Map().SameAs(*discret_->DofRowMap(0)))
    dserror("Maps are NOT identical. Emergency!");
#endif

  // Note: we can not simply write this into the secondary discretisation here
  // since it is a variable of the primary dofset and is hence cleared
  // in between
  mean_conc_ = MeanConc;
  return;
}  // ScaTraTimIntImpl::SetMeanConcentration


/*----------------------------------------------------------------------*
 | set convective velocity field (+ pressure and acceleration field as  |
 | well as fine-scale velocity field, if required)            gjb 05/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetVelocityField(
    Teuchos::RCP<const Epetra_Vector> convvel,  //!< convective velocity/press. vector
    Teuchos::RCP<const Epetra_Vector> acc,      //!< acceleration vector
    Teuchos::RCP<const Epetra_Vector> vel,      //!< velocity vector
    Teuchos::RCP<const Epetra_Vector> fsvel,    //!< fine-scale velocity vector
    const int nds,          //!< number of the dofset the velocity/pressure state belongs to
    const bool setpressure  //!< flag whether the fluid pressure needs to be known for the scatra
)
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA: set convective velocity field");

  //---------------------------------------------------------------------------
  // preliminaries
  //---------------------------------------------------------------------------

  if (convvel == Teuchos::null) dserror("Velocity state is Teuchos::null");

  if (velocity_field_type_ != INPAR::SCATRA::velocity_Navier_Stokes)
    dserror("Wrong SetVelocityField() called for velocity field type %d!", velocity_field_type_);

  if (nds >= discret_->NumDofSets()) dserror("Too few dofsets on scatra discretization!");

#ifdef DEBUG
  // We rely on the fact, that the nodal distribution of both fields is the same.
  // Although Scatra discretization was constructed as a clone of the fluid or
  // structure mesh, respectively, at the beginning, the nodal distribution may
  // have changed meanwhile (e.g., due to periodic boundary conditions applied only
  // to the fluid field)!
  // We have to be sure that everything is still matching.
  if (not convvel->Map().SameAs(*discret_->DofRowMap(nds)))
    dserror("Fluid/Structure and Scatra dofrowmaps are NOT identical. Emergency!");
#endif

  // boolean indicating whether fine-scale velocity vector exists
  // -> if yes, multifractal subgrid-scale modeling is applied
  bool fsvelswitch = (fsvel != Teuchos::null);

  // some thing went wrong if we want to use multifractal subgrid-scale modeling
  // and have not got the fine-scale velocity
  if (step_ >= 1 and
      (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales or
          fssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_small) and
      not fsvelswitch)
    dserror("Fine-scale velocity expected for multifractal subgrid-scale modeling!");
  // as fsvelswitch is also true for smagorinsky_all, we have to reset fsvelswitch
  // as the corresponding vector, which is not necessary, is not provided in scatra
  if (fssgd_ == INPAR::SCATRA::fssugrdiff_smagorinsky_all and fsvelswitch) fsvelswitch = false;
  // as fsvelswitch is true in case of turned-off model in scalar field,
  // we have to ensure false
  if (turbmodel_ == INPAR::FLUID::no_model and fssgd_ == INPAR::SCATRA::fssugrdiff_no)
    fsvelswitch = false;

  // store number of dofset associated with velocity related dofs
  nds_vel_ = nds;

  // provide scatra discretization with convective velocity
  discret_->SetState(nds_vel_, "convective velocity field", convvel);

  // provide scatra discretization with velocity
  if (vel != Teuchos::null)
    discret_->SetState(nds_vel_, "velocity field", vel);
  else
    // if velocity vector is not provided by the respective algorithm, we
    // assume that it equals the given convective velocity:
    discret_->SetState(nds_vel_, "velocity field", convvel);

  // provide scatra discretization with acceleration field if required
  if (acc != Teuchos::null) discret_->SetState(nds_vel_, "acceleration field", acc);

  // provide scatra discretization with fine-scale convective velocity if required
  if (fsvelswitch) discret_->SetState(nds_vel_, "fine-scale velocity field", fsvel);

  return;

}  // ScaTraTimIntImpl::SetVelocityField


/*----------------------------------------------------------------------*
 | contains the time loop                                       vg 05/07|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::TimeLoop()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:  + time loop");

  // prepare time loop
  PrepareTimeLoop();

  while (NotFinished())
  {
    // -------------------------------------------------------------------
    // prepare time step
    // -------------------------------------------------------------------
    PrepareTimeStep();

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
    if (DRT::INPUT::IntegralValue<int>(*params_, "OUTPUTNONLINSOLVERSTATS"))
      OutputNonlinSolverStats(iternum_, dtnonlinsolve, Step(), discret_->Comm());

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    Update();

    // -------------------------------------------------------------------
    // evaluate error for problems with analytical solution
    // -------------------------------------------------------------------
    EvaluateErrorComparedToAnalyticalSol();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();

  }  // while

  // print the results of time measurements
  Teuchos::TimeMonitor::summarize();

  return;
}  // ScaTraTimIntImpl::TimeLoop


/*----------------------------------------------------------------------*
 | contains the call of linear/nonlinear solver               gjb 02/10 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Solve()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // -----------------------------------------------------------------
  // intermediate solution step for homogeneous isotropic turbulence
  // -----------------------------------------------------------------
  if (solvtype_ == INPAR::SCATRA::solvertype_nonlinear) CalcIntermediateSolution();

  // -----------------------------------------------------------------
  //                     solve (non-)linear equation
  // -----------------------------------------------------------------
  switch (solvtype_)
  {
    case INPAR::SCATRA::solvertype_linear_incremental:
    case INPAR::SCATRA::solvertype_linear_full:
    {
      LinearSolve();
      break;
    }

    case INPAR::SCATRA::solvertype_nonlinear:
    {
      NonlinearSolve();
      break;
    }

    case INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro:
    case INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken:
    case INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit:
    case INPAR::SCATRA::solvertype_nonlinear_multiscale_microtomacro:
    {
      NonlinearMultiScaleSolve();
      break;
    }

    default:
    {
      dserror("Unknown solver type!");
      break;
    }
  }
  // that's all

  return;
}


/*------------------------------------------------------------------------------------------*
 | update solution after convergence of the nonlinear Newton-Raphson iteration   fang 01/17 |
 *------------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Update(const int num  //!< field number
)
{
  // update quantities associated with meshtying strategy
  strategy_->Update();

  return;
}


/*----------------------------------------------------------------------*
 | apply moving mesh data                                     gjb 05/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ApplyMeshMovement(Teuchos::RCP<const Epetra_Vector> dispnp, int nds)
{
  //---------------------------------------------------------------------------
  // only required in ALE case
  //---------------------------------------------------------------------------
  if (isale_)
  {
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA: apply mesh movement");

    // check existence of displacement vector
    if (dispnp == Teuchos::null) dserror("Got null pointer for displacements!");

    // store number of dofset associated with displacement related dofs
    nds_disp_ = nds;

    // provide scatra discretization with displacement field
    discret_->SetState(nds_disp_, "dispnp", dispnp);
  }  // if (isale_)

  return;
}  // ScaTraTimIntImpl::ApplyMeshMovement


/*----------------------------------------------------------------------*
 | print information about current time step to screen       fang 08/17 |
 *----------------------------------------------------------------------*/
inline void SCATRA::ScaTraTimIntImpl::PrintTimeStepInfo()
{
  if (myrank_ == 0)
    std::cout << std::endl
              << "TIME: " << std::setw(11) << std::setprecision(4) << time_ << "/" << maxtime_
              << "  DT = " << dta_ << "  " << MethodTitle() << std::setw(4) << "  STEP = " << step_
              << "/" << stepmax_ << std::endl;
}  // SCATRA::ScaTraTimIntImpl::PrintTimeStepInfo


/*----------------------------------------------------------------------*
 | return system matrix downcasted as sparse matrix           gjb 02/11 |
 | implemented here to be able to use forward declaration in .H         |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::SparseMatrix> SCATRA::ScaTraTimIntImpl::SystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
}


/*----------------------------------------------------------------------*
 | return system matrix downcasted as block sparse matrix     gjb 06/10 |
 | implemented here to be able to use forward declaration in .H         |
 *----------------------------------------------------------------------*/
Teuchos::RCP<LINALG::BlockSparseMatrixBase> SCATRA::ScaTraTimIntImpl::BlockSystemMatrix()
{
  return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_);
}


/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                          gjb 08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::Output(const int num)
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:    + output of solution");

  // solution output and potentially restart data and/or flux data
  if (DoOutput())
  {
    // step number and time (only after that data output is possible)
    output_->NewStep(step_, time_);

    // write domain decomposition for visualization (only once at the first time step!)
    if (step_ == 0) output_->WriteElementData(true);

    // write state vectors
    OutputState();

    // write output to Gmsh postprocessing files
    if (outputgmsh_) OutputToGmsh(step_, time_);

    // write flux vector field (only writing, calculation was done during Update() call)
    if (calcflux_domain_ != INPAR::SCATRA::flux_none or
        calcflux_boundary_ != INPAR::SCATRA::flux_none)
    {
      // for flux output of initial field (before first solve) do:
      // flux_domain_ and flux_boundary_ vectors are initialized when CalcFlux() is called
      if (step_ == 0 or
          (calcflux_domain_ != INPAR::SCATRA::flux_none and flux_domain_ == Teuchos::null) or
          (calcflux_boundary_ != INPAR::SCATRA::flux_none and flux_boundary_ == Teuchos::null))
        CalcFlux(true, num);

      if (calcflux_domain_ != INPAR::SCATRA::flux_none) OutputFlux(flux_domain_, "domain");
      if (calcflux_boundary_ != INPAR::SCATRA::flux_none) OutputFlux(flux_boundary_, "boundary");
    }

    // write mean values of scalar(s)
    OutputTotalAndMeanScalars(num);

    // write domain and boundary integrals, i.e., surface areas and volumes of specified nodesets
    OutputDomainOrBoundaryIntegrals("DomainIntegral");
    OutputDomainOrBoundaryIntegrals("BoundaryIntegral");

    // write integral values of reaction(s)
    OutputIntegrReac(num);

    // problem-specific outputs
    OutputProblemSpecific();

    // add restart data
    if (step_ % uprestart_ == 0 and step_ != 0) OutputRestart();

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
  }

  // generate output on micro scale if necessary
  if (macro_scale_)
  {
    // create parameter list for macro elements
    Teuchos::ParameterList eleparams;

    // set action
    eleparams.set<int>("action", SCATRA::micro_scale_output);

    // loop over macro-scale elements
    discret_->Evaluate(
        eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
  }

  if ((step_ != 0) and (output_state_matlab_))
  {
    std::ostringstream filename;
    filename << "Result_Step" << step_ << ".m";
    LINALG::PrintVectorInMatlabFormat(filename.str(), *phinp_);
  }
  // NOTE:
  // statistics output for normal fluxes at boundaries was already done during Update()

  return;
}  // ScaTraTimIntImpl::Output


/*==========================================================================*/
// scalar degrees of freedom and related
/*==========================================================================*/

/*----------------------------------------------------------------------*
 |  set initial field for phi                                 gjb 04/08 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetInitialField(
    const INPAR::SCATRA::InitialField init, const int startfuncno)
{
  switch (init)
  {
    case INPAR::SCATRA::initfield_zero_field:
    {
      phin_->PutScalar(0.0);
      phinp_->PutScalar(0.0);
      break;
    }
    case INPAR::SCATRA::initfield_field_by_function:
    case INPAR::SCATRA::initfield_disturbed_field_by_function:
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        int numdofs = nodedofset.size();
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);
          // evaluate component k of spatial function
          double initialval = problem_->Funct(startfuncno - 1).Evaluate(k, lnode->X(), time_);
          int err = phin_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) dserror("dof not on proc");
        }
      }

      // for NURBS discretizations we have to solve a least squares problem,
      // with high accuracy! (do nothing for Lagrangian polynomials)
      const Teuchos::ParameterList& scatradyn = problem_->ScalarTransportDynamicParams();
      const int lstsolver = scatradyn.get<int>("LINEAR_SOLVER");

      DRT::NURBS::NurbsDiscretization* nurbsdis =
          dynamic_cast<DRT::NURBS::NurbsDiscretization*>(&(*discret_));
      if (nurbsdis != NULL)
      {
        if (lstsolver == (-1))
          dserror(
              "no linear solver defined for least square NURBS problem. Please set LINEAR_SOLVER "
              "in SCALAR TRANSPORT DYNAMIC to a valid number! Note: this solver block is misused "
              "for the least square problem. Maybe one should add a separate parameter for this.");

        DRT::NURBS::apply_nurbs_initial_condition(
            *discret_, errfile_, problem_->SolverParams(lstsolver), startfuncno, phin_);
      }

      // initialize also the solution vector. These values are a pretty good guess for the
      // solution after the first time step (much better than starting with a zero vector)
      phinp_->Update(1.0, *phin_, 0.0);

      // add random perturbation for initial field of turbulent flows
      if (init == INPAR::SCATRA::initfield_disturbed_field_by_function)
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
        if (err > 0) dserror("Error during evaluation of maximum value.");
        err = phinp_->MinValue(&minphi);
        if (err > 0) dserror("Error during evaluation of minimum value.");
        double range = abs(maxphi - minphi);

        // disturb initial field for all degrees of freedom
        for (int k = 0; k < phinp_->MyLength(); ++k)
        {
          double randomnumber = problem_->Random()->Uni();
          double noise = perc * range * randomnumber;
          err += phinp_->SumIntoMyValues(1, &noise, &k);
          err += phin_->SumIntoMyValues(1, &noise, &k);
          if (err != 0) dserror("Error while disturbing initial field.");
        }
      }
      break;
    }
    case INPAR::SCATRA::initfield_field_by_condition:
    {
      // set initial field for ALL existing scatra fields in condition
      const std::string field = "ScaTra";

      // get initial field conditions
      std::vector<DRT::Condition*> initfieldconditions(0);
      discret_->GetCondition("Initfield", initfieldconditions);

      if (not initfieldconditions.size())
        dserror(
            "Tried to evaluate initial field by condition without a corresponding condition "
            "defined on the ScaTra discretization!");
      if (scalarhandler_ == Teuchos::null) dserror("scalarhandler_ is null pointer!");

      std::set<int> numdofpernode;
      for (unsigned icond = 0; icond < initfieldconditions.size(); icond++)
      {
        const int condmaxnumdofpernode =
            scalarhandler_->NumDofPerNodeInCondition(*(initfieldconditions[icond]), discret_);

        if (condmaxnumdofpernode != 0) numdofpernode.insert(condmaxnumdofpernode);
      }

      if (numdofpernode.empty()) dserror("No DOFs defined on initial field condition!");

      const int maxnumdofpernode = *(numdofpernode.rbegin());

      std::vector<int> localdofs(maxnumdofpernode);
      for (int i = 0; i < maxnumdofpernode; i++)
      {
        localdofs[i] = i;
      }
      discret_->EvaluateInitialField(field, phin_, localdofs);

      // initialize also the solution vector. These values are a pretty good guess for the
      // solution after the first time step (much better than starting with a zero vector)
      phinp_->Update(1.0, *phin_, 0.0);

      break;
    }
    // discontinuous 0-1 field for progress variable in 1-D
    case INPAR::SCATRA::initfield_discontprogvar_1D:
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        // get coordinate
        const double x = lnode->X()[0];

        int numdofs = nodedofset.size();
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);

          double initialval = 0.0;
          if (x > -EPS10) initialval = 1.0;

          int err = 0;
          err += phin_->ReplaceMyValues(1, &initialval, &doflid);
          // initialize also the solution vector. These values are a pretty good guess for the
          // solution after the first time step (much better than starting with a zero vector)
          err += phinp_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      break;
    }
    // reconstructed initial profile for progress variable in x2-direction from
    // Lessani and Papalexandris (2006), also used in Moureau et al. (2007, 2009),
    // for two-dimensional flame-vortex interaction problem (x2=0-200)
    case INPAR::SCATRA::initfield_flame_vortex_interaction:
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

      const Epetra_Map* dofrowmap = discret_->DofRowMap();

      // define variable
      double initialval = 0.0;

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        // get x2-coordinate
        const double x2 = lnode->X()[1];

        int numdofs = nodedofset.size();
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);

          if (x2 < loc12 - EPS10)
            initialval = (1.0 - (1.0 / beta1)) * exp((x2 - trans1) / delta1);
          else if (x2 > loc23 + EPS10)
            initialval = 1.0 - (exp((1.0 - beta3) * (x2 - trans3) / delta3) / beta3);
          else
            initialval = fac2 * (x2 - trans2) + abs2;

          int err = 0;
          err += phin_->ReplaceMyValues(1, &initialval, &doflid);
          // initialize also the solution vector. These values are a pretty good guess for the
          // solution after the first time step (much better than starting with a zero vector)
          err += phinp_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      break;
    }
    // initial mixture-fraction profile for Rayleigh-Taylor instability
    case INPAR::SCATRA::initfield_raytaymixfrac:
    {
      // define interface thickness, sinusoidal disturbance wave amplitude and pi
      const double delta = 0.002;
      const double alpha = 0.001;

      const Epetra_Map* dofrowmap = discret_->DofRowMap();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        // get x1- and x2-coordinate
        const double x1 = lnode->X()[0];
        const double x2 = lnode->X()[1];

        // interface disturbance
        // double x2_int = 0.05*cos(pi*(x1+0.5));
        // double x2_int = 0.05*cos(2.0*pi*x1);
        double x2_int = 0.0;
        x2_int -= cos(4 * M_PI * x1);
        x2_int -= cos(14 * M_PI * x1);
        x2_int -= cos(23 * M_PI * x1);
        x2_int -= cos(28 * M_PI * x1);
        x2_int -= cos(33 * M_PI * x1);
        x2_int -= cos(42 * M_PI * x1);
        x2_int -= cos(51 * M_PI * x1);
        x2_int -= cos(59 * M_PI * x1);
        x2_int *= alpha;

        const double value = (x2_int - x2) / (2.0 * delta);

        // values required for tanh-distribution
        const double vp = exp(value);
        const double vm = exp(-value);

        int numdofs = nodedofset.size();
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
          if (err != 0) dserror("dof not on proc");
        }
      }
      break;
    }
    // initial field for skew convection of L-shaped domain
    case INPAR::SCATRA::initfield_Lshapeddomain:
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        // get x1- and x2-coordinate
        const double x1 = lnode->X()[0];
        const double x2 = lnode->X()[1];

        int numdofs = nodedofset.size();
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
          if (err != 0) dserror("dof not on proc");
        }
      }
      break;
    }
    case INPAR::SCATRA::initfield_facing_flame_fronts:
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        // get x1- and x2-coordinate
        const double x1 = lnode->X()[0];
        // const double x2 = lnode->X()[1];

        int numdofs = nodedofset.size();
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
          if (err != 0) dserror("dof not on proc");
        }
      }
      break;
    }
    case INPAR::SCATRA::initfield_oracles_flame:
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap();

      const double eps = 0.00152;
      // const double xsing = 0.2;
      // const double zsing = 0.7525-0.05;//0.0354;

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        DRT::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        // get x1, x2 and x3-coordinate
        // const double x1 = lnode->X()[0];
        const double x2 = lnode->X()[1];
        // const double x3 = lnode->X()[2];

        int numdofs = nodedofset.size();
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

#if 0
        // initial wedge implementation for periodic spanwise boundary
        if (x1 <= 0.0)
        {
          if (x2 >= 0.0)
            initval = (x2-0.0354) - eps;
          else
            initval = (-0.0354-x2) - eps;
        }
        else if (x1 > 0.0 and x1 < xsing)
        {
          initval = abs(x2)-0.0354*(xsing-x1)/xsing - eps;
        }
        else if (x1 >= xsing)
          initval = x1 - xsing - eps;
        else
          dserror("impossible!");
#endif

#if 0
        // initial wedge implementation for spanwise walls
        if (x1 <= 0.0)
        {
          if ( x3 <= -zsing and abs(x2) <= abs(x3+zsing) )
          {
            initval = (-0.7525-x3) - eps;
          }
          else if ( x3 >= zsing and abs(x2) <= (x3-zsing) )
          {
            initval = (x3-0.7525) - eps;
          }
          else if ( x2 >= 0.0 and ( x2 > abs(x3+zsing) or x2 > (x3-zsing) ))
          {
            initval = (x2-0.0354) - eps;
          }
          else if ( x2 < 0.0 and (-x2 > abs(x3+zsing) or -x2 > (x3-zsing) ))
          {
            initval = (-0.0354-x2) - eps;
          }
          else
            dserror("coordinate out of range of ORACLES initial function");
        }
        else if (x1 > 0.0 and x1 < xsing)
        {
          if (abs(x3) <= 0.07)
            initval = abs(x2)-0.0354*(xsing-x1)/xsing - eps;
          else
          {
            initval = 0.07525-0.07;
          }
        }
        else if (x1 >= xsing)
          initval = x1 - xsing - eps;
        else
          dserror("impossible!");
#endif
          int err = 0;
          err += phin_->ReplaceMyValues(1, &initval, &doflid);
          // initialize also the solution vector. These values are a pretty good guess for the
          // solution after the first time step (much better than starting with a zero vector)
          err += phinp_->ReplaceMyValues(1, &initval, &doflid);
          if (err != 0) dserror("dof not on proc");
        }
      }
      break;
    }
    case INPAR::SCATRA::initialfield_forced_hit_high_Sc:
    case INPAR::SCATRA::initialfield_forced_hit_low_Sc:
    {
      // initialize calculation of initial field based on fast Fourier transformation
      Teuchos::RCP<HomIsoTurbInitialScalarField> HitInitialScalarField =
          Teuchos::rcp(new SCATRA::HomIsoTurbInitialScalarField(*this, init));
      // calculate initial field
      HitInitialScalarField->CalculateInitialField();

      break;
    }
    default:
      dserror("Unknown option for initial field: %d", init);
      break;
  }  // switch(init)

  return;
}  // ScaTraTimIntImpl::SetInitialField


/*----------------------------------------------------------------------*
 | iterative update of concentrations                                   |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc)
{
  // store incremental vector to be available for convergence check
  // if incremental vector is received from outside for coupled problem
  increment_->Update(1.0, *inc, 0.0);

  // update scalar values by adding increments
  phinp_->Update(1.0, *inc, 1.0);
}  // UpdateIter


/*==========================================================================*
 |                                                                          |
 | protected:                                                               |
 |                                                                          |
 *==========================================================================*/

/*==========================================================================*/
// general framework
/*==========================================================================*/

/*--------------------------------------------------------------------------*
 | setup Krylov projector including first fill                    nis Feb13 |
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetupKrylovSpaceProjection(DRT::Condition* kspcond)
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
  const int nummodes = kspcond->GetInt("NUMMODES");
  if (nummodes != NumDofPerNode())
    dserror(
        "Expecting as many mode flags as nodal dofs in Krylov projection definition. Check "
        "dat-file!");

  // get vector of mode flags as given in dat-file
  const std::vector<int>* modeflags = kspcond->Get<std::vector<int>>("ONOFF");

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
  const std::string* weighttype = kspcond->Get<std::string>("weight vector definition");

  // set flag for projection update true only if ALE and integral weights
  if (isale_ and (*weighttype == "integration")) updateprojection_ = true;

  // create the projector
  projector_ =
      Teuchos::rcp(new LINALG::KrylovProjector(activemodeids, weighttype, discret_->DofRowMap()));

  // update the projector
  UpdateKrylovSpaceProjection();
}  // SCATRA::ScaTraTimIntImpl::SetupKrylovSpaceProjection


/*--------------------------------------------------------------------------*
 | update projection vectors w_ and c_ for Krylov projection      nis Feb13 |
 *--------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::UpdateKrylovSpaceProjection()
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
    dserror("option pointvalues not implemented");
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
    Epetra_IntSerialDenseVector dofids(NumDofPerNode());
    for (int rr = 0; rr < NumDofPerNode(); ++rr)
    {
      dofids[rr] = -1;
    }

    Teuchos::ParameterList mode_params;

    // set parameters for elements that do not change over mode
    mode_params.set<int>("action", SCATRA::integrate_shape_functions);
    if (isale_) mode_params.set<int>("ndsdisp", nds_disp_);

    // loop over all activemodes
    for (int imode = 0; imode < nummodes; ++imode)
    {
      // activate dof of current mode and add dofids to parameter list
      dofids[modeids[imode]] = 1;
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
      discret_->EvaluateCondition(mode_params, Teuchos::null, Teuchos::null, wi, Teuchos::null,
          Teuchos::null, "KrylovSpaceProjection");

      // deactivate dof of current mode
      dofids[modeids[imode]] = -1;

      // set the current kernel basis vector - not very nice
      for (int inode = 0; inode < discret_->NumMyRowNodes(); inode++)
      {
        DRT::Node* node = discret_->lRowNode(inode);
        std::vector<int> gdof = discret_->Dof(0, node);
        int err = c->ReplaceGlobalValue(gdof[modeids[imode]], imode, 1);
        if (err != 0) dserror("error while inserting value into c");
      }

    }  // loop over modes

    // adapt weight vector according to meshtying case
    if (msht_ != INPAR::FLUID::no_meshtying)
    {
      dserror(
          "Since meshtying for scatra is not tested under Krylov projection dserror is introduced. "
          "Remove at own responsibility.");
      // meshtying_->AdaptKrylovProjector(w);
    }

  }  // endif integration
  else
  {
    dserror("unknown definition of weight vector w for restriction of Krylov space");
  }

  // adapt kernel vector according to meshtying case
  if (msht_ != INPAR::FLUID::no_meshtying)
  {
    dserror(
        "Since meshtying for scatra is not tested under Krylov projection dserror is introduced. "
        "Remove at own responsibility.");
    // meshtying_->AdaptKrylovProjector(c);
  }

  // fillcomplete the projector to compute (w^T c)^(-1)
  projector_->FillComplete();

  return;
}  // ScaTraTimIntImpl::UpdateKrylovSpaceProjection


/*----------------------------------------------------------------------------------------*
 | initialize meshtying strategy (including standard case without meshtying)   fang 12/14 |
 *----------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::CreateMeshtyingStrategy()
{
  // fluid meshtying
  if (msht_ != INPAR::FLUID::no_meshtying)
    strategy_ = Teuchos::rcp(new MeshtyingStrategyFluid(this));

  // scatra-scatra interface coupling
  else if (IsS2IMeshtying())
    strategy_ = Teuchos::rcp(new MeshtyingStrategyS2I(this, *params_));

  // scatra-scatra interface coupling
  else if (heteroreaccoupling_)
    strategy_ = Teuchos::rcp(new HeterogeneousReactionStrategy(this));

  else if (arterycoupling_)
    strategy_ = Teuchos::rcp(new MeshtyingStrategyArtery(this));

  // standard case without meshtying
  else
    strategy_ = Teuchos::rcp(new MeshtyingStrategyStd(this));

  return;
}  // ScaTraTimIntImpl::CreateMeshtyingStrategy

/*----------------------------------------------------------------------------------------*
 | create scalar manager                                                      fang 12/14 |
 *----------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::CreateScalarHandler()
{
  scalarhandler_ = Teuchos::rcp(new ScalarHandler());

  return;
}  // ScaTraTimIntImpl::CreateScalarHandler


/*----------------------------------------------------------------------*
 | add approximation to flux vectors to a parameter list      gjb 05/10 |
 *----------------------------------------------------------------------*/
// void SCATRA::ScaTraTimIntImpl::AddFluxApproxToParameterList
// defined in scalar_timint_implicit_service.cpp

/*----------------------------------------------------------------------*
 | application of Dirichlet boundary conditions                         |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ApplyDirichletToSystem()
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

      LINALG::ApplyDirichlettoSystem(
          sysmat_, increment_, residual_, zeros_, *(dbcmaps_->CondMap()));
    }
  }
  else
  {
    // time measurement: application of DBC to system
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + apply DBC to system");

    LINALG::ApplyDirichlettoSystem(sysmat_, phinp_, residual_, phinp_, *(dbcmaps_->CondMap()));
  }
  return;
}  // SCATRA::ScaTraTimIntImpl::ApplyDirichletToSystem


/*----------------------------------------------------------------------*
 | evaluate Dirichlet boundary conditions at t_{n+1}           gjb 07/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ApplyDirichletBC(
    const double time, Teuchos::RCP<Epetra_Vector> phinp, Teuchos::RCP<Epetra_Vector> phidt)
{
  // time measurement: apply Dirichlet conditions
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:      + apply dirich cond.");

  // Todo: what happens in  the case of generalized alpha
  // needed parameters
  Teuchos::ParameterList p;
  p.set("total time", time);  // actual time t_{n+1}

  // predicted Dirichlet values
  // \c  phinp then also holds prescribed new Dirichlet values
  discret_->ClearState();
  discret_->EvaluateDirichlet(p, phinp, phidt, Teuchos::null, Teuchos::null, dbcmaps_);
  discret_->ClearState();

  return;
}  // SCATRA::ScaTraTimIntImpl::ApplyDirichletBC


/*----------------------------------------------------------------------*
 | contains the residual scaling and addition of Neumann terms vg 08/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ScalingAndNeumann()
{
  // scaling to get true residual vector for all time integration schemes
  // in incremental case: boundary flux values can be computed from trueresidual
  if (incremental_) trueresidual_->Update(ResidualScaling(), *residual_, 0.0);

  // add Neumann b.c. scaled with a factor due to time discretization
  AddNeumannToResidual();

  // add potential Neumann inflow or convective heat transfer boundary
  // conditions (simultaneous evaluation of both conditions not allowed!)
  if (neumanninflow_)
    ComputeNeumannInflow(sysmat_, residual_);
  else if (convheatrans_)
    EvaluateConvectiveHeatTransfer(sysmat_, residual_);

  return;
}  // ScaTraTimIntImpl::ScalingAndNeumann


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary conditions                      fang 01/15 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ApplyNeumannBC(
    const Teuchos::RCP<Epetra_Vector>& neumann_loads  //!< Neumann loads
)
{
  // prepare load vector
  neumann_loads->PutScalar(0.0);

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action", SCATRA::bd_calc_Neumann);

  // specific parameters
  AddProblemSpecificParametersAndVectors(condparams);

  // set time for evaluation of point Neumann conditions as parameter depending on time integration
  // scheme line/surface/volume Neumann conditions use the time stored in the time parameter class
  SetTimeForNeumannEvaluation(condparams);

  // provide displacement field in case of ALE
  if (isale_) condparams.set<int>("ndsdisp", nds_disp_);

  // evaluate Neumann boundary conditions at time t_{n+alpha_F} (generalized alpha) or time t_{n+1}
  // (otherwise)
  discret_->EvaluateNeumann(condparams, *neumann_loads);
  discret_->ClearState();

  return;
}  // SCATRA::ScaTraTimIntImpl::ApplyNeumannBC


/*----------------------------------------------------------------------------*
 | evaluate solution-depending boundary and interface conditions   fang 10/14 |
 *----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::EvaluateSolutionDependingConditions(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,  //!< system matrix
    Teuchos::RCP<Epetra_Vector> rhs                     //!< rhs vector
)
{
  // evaluate Robin type boundary condition
  EvaluateRobinBoundaryConditions(systemmatrix, rhs);

  // evaluate meshtying
  // this needs to be done as final step for consistency
  strategy_->EvaluateMeshtying();

  // evaluate macro-micro coupling on micro scale in multi-scale scalar transport problems
  EvaluateMacroMicroCoupling();

  return;
}  // SCATRA::ScaTraTimIntImpl::EvaluateSolutionDependingConditions


/*----------------------------------------------------------------------------*
 | evaluate additional solution-depending models                  rauch 12/16 |
 *----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::EvaluateAdditionalSolutionDependingModels(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,  //!< system matrix
    Teuchos::RCP<Epetra_Vector> rhs                     //!< rhs vector
)
{
  // evaluate solution depending additional models
  // this point is unequal NULL only if a scatra
  // adapter has been constructed.
  if (additional_model_evaluator_ != NULL)
    additional_model_evaluator_->EvaluateAdditionalSolutionDependingModels(systemmatrix, rhs);

  return;
}  // SCATRA::ScaTraTimIntImpl::EvaluateAdditionalSolutionDependingModels


/*----------------------------------------------------------------------------*
 | evaluate Robin boundary conditions                          schoeder 03/15 |
 *----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::EvaluateRobinBoundaryConditions(
    Teuchos::RCP<LINALG::SparseOperator> matrix,  //!< system matrix
    Teuchos::RCP<Epetra_Vector> rhs               //!< rhs vector
)
{
  discret_->ClearState();

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action", SCATRA::bd_calc_Robin);

  // provide displacement field in case of ALE
  if (isale_) condparams.set<int>("ndsdisp", nds_disp_);

  // add element parameters and set state vectors according to time-integration scheme
  AddTimeIntegrationSpecificVectors();

  // evaluate ElchBoundaryKinetics conditions at time t_{n+1} or t_{n+alpha_F}
  discret_->EvaluateCondition(
      condparams, matrix, Teuchos::null, rhs, Teuchos::null, Teuchos::null, "TransportRobin");
  discret_->ClearState();

  return;
}  // ScaTraTimIntImpl::EvaluateRobinBoundaryConditions


/*----------------------------------------------------------------------*
 | contains the assembly process for matrix and rhs            vg 08/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AssembleMatAndRHS()
{
  // safety check
  CheckIsInit();
  CheckIsSetup();

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
  eleparams.set<int>("action", SCATRA::calc_mat_and_rhs);

  // DO THIS AT VERY FIRST!!!
  // compute reconstructed diffusive fluxes for better consistency
  const enum INPAR::SCATRA::Consistency consistency =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::Consistency>(
          params_->sublist("STABILIZATION"), "CONSISTENCY");
  if (consistency == INPAR::SCATRA::consistency_l2_projection_lumped)
  {
    // compute flux approximation and add it to the parameter list
    AddFluxApproxToParameterList(eleparams);
  }

  // prepare dynamic Smagorinsky model if required,
  // i.e. calculate turbulent Prandtl number
  if (timealgo_ != INPAR::SCATRA::timeint_stationary)
  {
    DynamicComputationOfCs();
    DynamicComputationOfCv();
  }
  // this parameter list is required here to get the element-based filtered constants
  eleparams.sublist("TURBULENCE MODEL") = extraparams_->sublist("TURBULENCE MODEL");

  // provide velocity field and potentially acceleration/pressure field
  // (export to column map necessary for parallel evaluation)
  eleparams.set<int>("ndsvel", nds_vel_);

  // provide displacement field in case of ALE
  if (isale_) eleparams.set<int>("ndsdisp", nds_disp_);

  // set vector values needed by elements
  discret_->ClearState();

  // AVM3 separation for incremental solver: get fine-scale part of scalar
  if (incremental_ and step_ > 0 and
      (fssgd_ != INPAR::SCATRA::fssugrdiff_no or
          turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales))
    AVM3Separation();

  // add state vectors according to time-integration scheme
  AddTimeIntegrationSpecificVectors();

  // set external volume force (required, e.g., for forced homogeneous isotropic turbulence)
  if (homisoturb_forcing_ != Teuchos::null) homisoturb_forcing_->UpdateForcing(step_);

  if (forcing_ != Teuchos::null) discret_->SetState("forcing", forcing_);

  // add problem specific time-integration parameters
  AddProblemSpecificParametersAndVectors(eleparams);

  // call loop over elements (with or without subgrid-diffusivity(-scaling) vector)
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no)
    discret_->Evaluate(eleparams, sysmat_, Teuchos::null, residual_, subgrdiff_, Teuchos::null);
  else
    discret_->Evaluate(eleparams, sysmat_, residual_);

  //  (SystemMatrix()->EpetraMatrix())->Print(std::cout); // kn nis

  discret_->ClearState();

  //----------------------------------------------------------------------
  // apply weak Dirichlet boundary conditions
  //----------------------------------------------------------------------
  {
    Teuchos::ParameterList mhdbcparams;

    // set action for elements
    mhdbcparams.set<int>("action", SCATRA::bd_calc_weak_Dirichlet);

    eleparams.set<int>("ndsvel", nds_vel_);
    AddTimeIntegrationSpecificVectors();

    // evaluate all mixed hybrid Dirichlet boundary conditions
    discret_->EvaluateCondition(mhdbcparams, sysmat_, Teuchos::null, residual_, Teuchos::null,
        Teuchos::null, "LineWeakDirichlet");

    discret_->EvaluateCondition(mhdbcparams, sysmat_, Teuchos::null, residual_, Teuchos::null,
        Teuchos::null, "SurfaceWeakDirichlet");

    // clear state
    discret_->ClearState();
  }

  // AVM3 scaling for non-incremental solver: scaling of normalized AVM3-based
  // fine-scale subgrid-diffusivity matrix by subgrid diffusivity
  if (not incremental_ and fssgd_ != INPAR::SCATRA::fssugrdiff_no) AVM3Scaling(eleparams);

  // potential residual scaling and potential addition of Neumann terms
  ScalingAndNeumann();  // TODO: do we have to call this function twice??

  // evaluate solution-depending additional models
  EvaluateAdditionalSolutionDependingModels(sysmat_, residual_);

  // evaluate solution-depending boundary and interface conditions
  EvaluateSolutionDependingConditions(sysmat_, residual_);

  // finalize assembly of system matrix
  sysmat_->Complete();

  // end time measurement for element and take average over all processors via communication
  double mydtele = Teuchos::Time::wallTime() - tcpuele;
  discret_->Comm().MaxAll(&mydtele, &dtele_, 1);

  return;
}  // ScaTraTimIntImpl::AssembleMatAndRHS


/*----------------------------------------------------------------------*
 | contains the linear solver                                  vg 08/09 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::LinearSolve()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // -------------------------------------------------------------------
  //                        output to screen
  // -------------------------------------------------------------------
  PrintTimeStepInfo();

  // -------------------------------------------------------------------
  //                     preparations for solve
  // -------------------------------------------------------------------
  PrepareLinearSolve();

  // -------------------------------------------------------------------
  // Solve system in incremental or non-incremental case
  // -------------------------------------------------------------------
  if (incremental_)
  {
    TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + solver calls");

    // get cpu time
    const double tcpusolve = Teuchos::Time::wallTime();

    strategy_->Solve(solver_, sysmat_, increment_, residual_, phinp_, 1, Teuchos::null);

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
        if (scalnorm_L2 > EPS10)
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

    strategy_->Solve(solver_, sysmat_, phinp_, residual_, phinp_, 1, Teuchos::null);

    // end time measurement for solver
    dtsolve_ = Teuchos::Time::wallTime() - tcpusolve;

    if (myrank_ == 0) printf("Solvertype linear_full (ts=%10.3E,te=%10.3E)\n", dtsolve_, dtele_);
  }

  // -------------------------------------------------------------------
  // compute values at intermediate time steps (only for gen.-alpha)
  // -------------------------------------------------------------------
  ComputeIntermediateValues();

  return;
}  // ScaTraTimIntImpl::LinearSolve


/*----------------------------------------------------------------------*
 | contains the nonlinear iteration loop                       gjb 09/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::NonlinearSolve()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // time measurement: nonlinear iteration
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:   + nonlin. iteration/lin. solve");

  // out to screen
  PrintTimeStepInfo();

  // special preparations for multifractal subgrid-scale model
  if (turbmodel_ == INPAR::FLUID::multifractal_subgrid_scales) RecomputeMeanCsgsB();

  //------------------------------ turn adaptive solver tolerance on/off
  const double ittol = params_->sublist("NONLINEAR").get<double>("CONVTOL");
  const bool isadapttol =
      (DRT::INPUT::IntegralValue<int>(params_->sublist("NONLINEAR"), "ADAPTCONV"));
  const double adaptolbetter = params_->sublist("NONLINEAR").get<double>("ADAPTCONV_BETTER");
  double actresidual(0.0);

  // prepare Newton-Raphson iteration
  iternum_ = 0;

  // perform explicit predictor step (-> better starting point for nonlinear solver)
  const bool explpredictor =
      (DRT::INPUT::IntegralValue<int>(params_->sublist("NONLINEAR"), "EXPLPREDICT") == 1);
  if (explpredictor) ExplicitPredictor();

  // start Newton-Raphson iteration
  while (true)
  {
    iternum_++;

    // call elements to calculate system matrix and rhs and assemble
    AssembleMatAndRHS();

    // perform finite difference check on time integrator level
    if (fdcheck_ == INPAR::SCATRA::fdcheck_global) FDCheck();

    // project residual such that only part orthogonal to nullspace is considered
    if (projector_ != Teuchos::null) projector_->ApplyPT(*residual_);

    // Apply Dirichlet boundary conditions to system of equations
    // residual values are supposed to be zero at Dirichlet boundaries
    {
      // time measurement: application of DBC to system
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + apply DBC to system");

      LINALG::ApplyDirichlettoSystem(
          sysmat_, increment_, residual_, zeros_, *(dbcmaps_->CondMap()));
    }

    // abort nonlinear iteration if desired
    if (strategy_->AbortNonlinIter(*this, actresidual)) break;

    // initialize increment vector
    increment_->PutScalar(0.0);

    {
      // get cpu time
      const double tcpusolve = Teuchos::Time::wallTime();

      // time measurement: call linear solver
      TEUCHOS_FUNC_TIME_MONITOR("SCATRA:       + call linear solver");

      // do adaptive linear solver tolerance (not in first solve)
      if (isadapttol && iternum_ > 1)
      {
        solver_->AdaptTolerance(ittol, actresidual, adaptolbetter);
      }

      // reprepare Krylov projection only if ale and projection required
      if (updateprojection_) UpdateKrylovSpaceProjection();

      strategy_->Solve(solver_, sysmat_, increment_, residual_, phinp_, iternum_, projector_);

      solver_->ResetTolerance();

      // end time measurement for solver and take average over all processors via communication
      double mydtsolve = Teuchos::Time::wallTime() - tcpusolve;
      discret_->Comm().MaxAll(&mydtsolve, &dtsolve_, 1);

      // output performance statistics associated with linear solver into text file if applicable
      if (DRT::INPUT::IntegralValue<int>(*params_, "OUTPUTLINSOLVERSTATS"))
        OutputLinSolverStats(strategy_->Solver(), dtsolve_, Step(), iternum_,
            strategy_->DofRowMap().NumGlobalElements());
    }

    //------------------------------------------------ update solution vector
    phinp_->Update(1.0, *increment_, 1.0);

    //-------- update values at intermediate time steps (only for gen.-alpha)
    ComputeIntermediateValues();

    // compute values at the interior of the elements (required for hdg)
    ComputeInteriorValues();

  }  // nonlinear iteration

  return;
}  // ScaTraTimIntImpl::NonlinearSolve


/*--------------------------------------------------------------------------------------------------*
 | contains the nonlinear iteration loop for truly partitioned multi-scale simulations   fang 08/17
 |
 *--------------------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::NonlinearMultiScaleSolve()
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
    if (solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro or
        solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken or
        solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit)
    {
      // backup macro-scale state vector
      const Teuchos::RCP<Epetra_Vector> phinp = phinp_;

      // replace macro-scale state vector by relaxed macro-scale state vector as input for micro
      // scale
      phinp_ = phinp_relaxed;

      // solve micro-scale problems
      NonlinearMicroScaleSolve();

      // undo state vector replacement
      phinp_ = phinp;

      // solve macro-scale problem
      NonlinearSolve();
    }

    // solve macro scale first and micro scale second
    else if (solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_microtomacro)
    {
      // solve macro-scale problem
      NonlinearSolve();

      // solve micro-scale problems
      NonlinearMicroScaleSolve();
    }

    // compute increment of macro-scale state vector
    phinp_inc_->Update(1., *phinp_, -1.);

    // convergence check
    if (strategy_->AbortOuterIter(*this)) break;

    if (solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken or
        solvtype_ == INPAR::SCATRA::solvertype_nonlinear_multiscale_macrotomicro_aitken_dofsplit)
    {
      // compute difference between current and previous increments of macro-scale state vector
      Epetra_Vector phinp_inc_diff(*phinp_inc_);
      phinp_inc_diff.Update(-1., *phinp_inc_old_, 1.);

      // perform Aitken relaxation
      PerformAitkenRelaxation(*phinp_relaxed, phinp_inc_diff);

      // update increment of macro-scale state vector
      phinp_inc_old_->Update(1., *phinp_inc_, 0.);
    }

    else
      // no relaxation
      phinp_relaxed = phinp_;
  }

  return;
}  // SCATRA::ScaTraTimIntImpl::NonlinearMultiScaleSolve


/*-----------------------------------------------------------------------------*
 | solve micro scale in truly partitioned multi-scale simulations   fang 08/17 |
 *-----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::NonlinearMicroScaleSolve()
{
  // clear macro-scale discretization
  discret_->ClearState();

  // initialize parameter list for evaluation of macro-scale elements
  Teuchos::ParameterList eleparams;

  // set action for macro-scale elements
  eleparams.set<int>("action", SCATRA::micro_scale_solve);

  // set state vectors
  AddTimeIntegrationSpecificVectors();

  // evaluate macro-scale elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  // clear macro-scale discretization
  discret_->ClearState();

  return;
}  // SCATRA::ScaTraTimIntImpl::NonlinearMicroScaleSolve


/*--------------------------------------------------------------------------*
| returns matching std::string for each time integration scheme   gjb 08/08 |
*---------------------------------------------------------------------------*/
std::string SCATRA::ScaTraTimIntImpl::MapTimIntEnumToString(
    const enum INPAR::SCATRA::TimeIntegrationScheme term)
{
  // length of return std::string is 14 due to usage in formated screen output
  switch (term)
  {
    case INPAR::SCATRA::timeint_one_step_theta:
      return "One-Step-Theta";
      break;
    case INPAR::SCATRA::timeint_bdf2:
      return "     BDF2     ";
      break;
    case INPAR::SCATRA::timeint_stationary:
      return "  Stationary  ";
      break;
    case INPAR::SCATRA::timeint_gen_alpha:
      return "  Gen. Alpha  ";
      break;
    default:
      dserror("Cannot cope with name enum %d", term);
      return "";
      break;
  }

  return "";
}  // ScaTraTimIntImpl::MapTimIntEnumToString


/*----------------------------------------------------------------------*
 |  write current state to BINIO                             gjb   08/08|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::OutputState()
{
  // solution
  output_->WriteVector("phinp", phinp_);

  // convective velocity (written in case of coupled simulations since volmortar is now possible)
  if (velocity_field_type_ == INPAR::SCATRA::velocity_function or
      velocity_field_type_ == INPAR::SCATRA::velocity_Navier_Stokes)
  {
    Teuchos::RCP<const Epetra_Vector> convel =
        discret_->GetState(nds_vel_, "convective velocity field");
    if (convel == Teuchos::null) dserror("Cannot get state vector convective velocity");

    // convert dof-based Epetra vector into node-based Epetra multi-vector for postprocessing
    Teuchos::RCP<Epetra_MultiVector> convel_multi =
        Teuchos::rcp(new Epetra_MultiVector(*discret_->NodeRowMap(), nsd_, true));
    for (int inode = 0; inode < discret_->NumMyRowNodes(); ++inode)
    {
      DRT::Node* node = discret_->lRowNode(inode);
      for (int idim = 0; idim < nsd_; ++idim)
        (*convel_multi)[idim][inode] =
            (*convel)[convel->Map().LID(discret_->Dof(nds_vel_, node, idim))];
    }

    output_->WriteVector("convec_velocity", convel_multi, IO::nodevector);
  }

  // displacement field
  if (isale_)
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = discret_->GetState(nds_disp_, "dispnp");
    if (dispnp == Teuchos::null) dserror("Cannot extract displacement field from discretization");

    // convert dof-based Epetra vector into node-based Epetra multi-vector for postprocessing
    Teuchos::RCP<Epetra_MultiVector> dispnp_multi =
        Teuchos::rcp(new Epetra_MultiVector(*discret_->NodeRowMap(), nsd_, true));
    for (int inode = 0; inode < discret_->NumMyRowNodes(); ++inode)
    {
      DRT::Node* node = discret_->lRowNode(inode);
      for (int idim = 0; idim < nsd_; ++idim)
        (*dispnp_multi)[idim][inode] =
            (*dispnp)[dispnp->Map().LID(discret_->Dof(nds_disp_, node, idim))];
    }

    output_->WriteVector("dispnp", dispnp_multi, IO::nodevector);
  }

  return;
}  // ScaTraTimIntImpl::OutputState


/*----------------------------------------------------------------------*
 | increment time and step for next iteration                     mr. x |
 *----------------------------------------------------------------------*/
inline void SCATRA::ScaTraTimIntImpl::IncrementTimeAndStep()
{
  step_ += 1;
  time_ += dta_;
}


/*----------------------------------------------------------------------*
 | adapt time step size if desired                           fang 02/18 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AdaptTimeStepSize()
{
  // check flag for adaptive time stepping
  if (DRT::INPUT::IntegralValue<bool>(*params_, "ADAPTIVE_TIMESTEPPING"))
  {
    // initialize time step size with original value
    double dt(params_->get<double>("TIMESTEP"));

    // reduce time step size if necessary
    ComputeTimeStepSize(dt);

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
      SetDt(dt);
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | compute time step size                                    fang 02/18 |
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::ComputeTimeStepSize(double& dt)
{
  strategy_->ComputeTimeStepSize(dt);
  return;
}


/*----------------------------------------------------------------------*
 |  add dirichlet dofs to dbcmaps_                         rauch   04/15|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd)
{
  std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
  condmaps.push_back(maptoadd);
  condmaps.push_back(dbcmaps_->CondMap());
  Teuchos::RCP<Epetra_Map> condmerged = LINALG::MultiMapExtractor::MergeMaps(condmaps);
  *dbcmaps_ = LINALG::MapExtractor(*(discret_->DofRowMap()), condmerged);
  return;
}  // ScaTraTimIntImpl::AddDirichCond


/*----------------------------------------------------------------------------*
 | add global state vectors specific for time-integration scheme   fang 01/17 |
 *----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver)
{
  // add global state vectors associated with meshtying strategy
  strategy_->AddTimeIntegrationSpecificVectors();

  return;
}


/*----------------------------------------------------------------------*
 |  remove dirichlet dofs from dbcmaps_                    rauch   04/15|
 *----------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove)
{
  std::vector<Teuchos::RCP<const Epetra_Map>> othermaps;
  othermaps.push_back(maptoremove);
  othermaps.push_back(dbcmaps_->OtherMap());
  Teuchos::RCP<Epetra_Map> othermerged = LINALG::MultiMapExtractor::MergeMaps(othermaps);
  *dbcmaps_ = LINALG::MapExtractor(*(discret_->DofRowMap()), othermerged, false);
  return;
}  // ScaTraTimIntImpl::RemoveDirichCond


/*----------------------------------------------------------------------*
 |  return pointer to const dofrowmap                      rauch   04/15|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SCATRA::ScaTraTimIntImpl::DofRowMap() { return DofRowMap(0); }


/*----------------------------------------------------------------------*
 |  return pointer to const dofrowmap of specified dofset  rauch   04/15|
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> SCATRA::ScaTraTimIntImpl::DofRowMap(int nds)
{
  const Epetra_Map* dofrowmap = discret_->DofRowMap(nds);
  return Teuchos::rcp(dofrowmap, false);
}

/*----------------------------------------------------------------------*
 |  return number of scalars in scatra discretization       vuong   04/16|
 *----------------------------------------------------------------------*/
int SCATRA::ScaTraTimIntImpl::NumScal() const
{
  if (scalarhandler_ == Teuchos::null) dserror("scalar handler was not initialized!");
  return scalarhandler_->NumScal();
}

/*----------------------------------------------------------------------*
 |  return number of DOFs in scatra discretization          vuong   04/16|
 *----------------------------------------------------------------------*/
int SCATRA::ScaTraTimIntImpl::NumDofPerNode() const
{
  if (scalarhandler_ == Teuchos::null) dserror("scalar handler was not initialized!");
  return scalarhandler_->NumDofPerNode();
}

/*----------------------------------------------------------------------*
|  return number of dofs per node in condition
*----------------------------------------------------------------------*/
int SCATRA::ScaTraTimIntImpl::NumDofPerNodeInCondition(const DRT::Condition& condition) const
{
  if (scalarhandler_ == Teuchos::null) dserror("scalar handler was not initialized!");
  return scalarhandler_->NumDofPerNodeInCondition(condition, discret_);
}

/*-----------------------------------------------------------------------------*
 |  return total values of transported scalars (for output only)   vuong   04/16|
 *-----------------------------------------------------------------------------*/
const std::map<const int, std::vector<double>>& SCATRA::ScaTraTimIntImpl::TotalScalars() const
{
  if (outputscalarstrategy_ == Teuchos::null) dserror("output strategy was not initialized!");

  return outputscalarstrategy_->TotalScalars();
}

/*-----------------------------------------------------------------------------*
 |  return mean values of transported scalars (for output only)   vuong   04/16|
 *-----------------------------------------------------------------------------*/
const std::map<const int, std::vector<double>>& SCATRA::ScaTraTimIntImpl::MeanScalars() const
{
  if (outputscalarstrategy_ == Teuchos::null) dserror("output strategy was not initialized!");

  return outputscalarstrategy_->MeanScalars();
}

/*-----------------------------------------------------------------------------*
 |  return values of domain integrals (for output only)       kremheller 11/19 |
 *-----------------------------------------------------------------------------*/
const std::vector<double>& SCATRA::ScaTraTimIntImpl::DomainIntegrals() const
{
  if (outputdomainintegralstrategy_ == Teuchos::null)
    dserror("output strategy for domain integration was not initialized!");

  return outputdomainintegralstrategy_->DomainIntegrals();
}

/*-----------------------------------------------------------------------------*
 |  return values of boundary integrals (for output only)     kremheller 11/19 |
 *-----------------------------------------------------------------------------*/
const std::vector<double>& SCATRA::ScaTraTimIntImpl::BoundaryIntegrals() const
{
  if (outputdomainintegralstrategy_ == Teuchos::null)
    dserror("output strategy for domain integration was not initialized!");

  return outputdomainintegralstrategy_->BoundaryIntegrals();
}

void SCATRA::ScaTraTimIntImpl::GetPointPhiValue(
    const std::vector<double>& point, double& value, bool evalreac, unsigned int numscal)
{
  if (int(numscal) > this->NumScal())
    dserror("you requested the point value for scalar %d but there is only %d scalar fields",
        numscal, this->NumScal());

  int dim = point.size();
  if (dim > 3 || dim < 1)
    dserror("point has less than 1 or more than 3 coordinates which is quite uncommon");

  double lvalue = 0.0;

  //{ first step: find the element in which the point is located
  double bestdistance = 1.0e-10;
  double currentdistance = std::numeric_limits<double>::max();
  int lid_closestnode = -1;
  double lh = (discret_->lRowNode(0)->X()[0] - discret_->lRowNode(1)->X()[0]) *
              (discret_->lRowNode(0)->X()[0] - discret_->lRowNode(1)->X()[0]);
  for (int d = 1; d < dim; ++d)
    lh += (discret_->lRowNode(0)->X()[d] - discret_->lRowNode(1)->X()[d]) *
          (discret_->lRowNode(0)->X()[d] - discret_->lRowNode(1)->X()[d]);
  lh = std::sqrt(lh);
  double gh = 0.0;
  discret_->Comm().MaxAll(&lh, &gh, 1);

  // find closest node
  for (int n = 0; n < discret_->NodeRowMap()->NumMyElements(); ++n)
  {
    std::vector<double> nodecoords(dim, 0.0);
    for (int i = 0; i < dim; ++i) nodecoords[i] = discret_->lRowNode(n)->X()[i];

    double distance = distance_between_two_points(point, nodecoords);
    if (distance < currentdistance)
    {
      currentdistance = distance;
      lid_closestnode = n;
    }
    if (currentdistance < bestdistance) break;
  }
  // find closest point across processors
  double globbestdistance = 0.0;
  discret_->Comm().MinAll(&currentdistance, &globbestdistance, 1);

  if (globbestdistance > 10 * gh)
  {
    // std::cout<<"warning: you called PointValue on a point which is not inside of the Scatra
    // domain"<<std::endl;
    lvalue = 0.0;
  }
  else
  {
    int owner = -1;
    if (globbestdistance <= currentdistance + 1.0e-13 &&
        globbestdistance >= currentdistance - 1.0e-13)
      owner = myrank_;

    int globowner = -1;
    discret_->Comm().MaxAll(&owner, &globowner, 1);

    discret_->SetState("phinp", phinp_);

    // check adjacent elements to node
    int anywhere_inside = 0;
    if (myrank_ == globowner)  // this process has the closest node
    {
      int numele_adjacent_to_node = discret_->lRowNode(lid_closestnode)->NumElement();
      for (int e = 0; e < numele_adjacent_to_node; ++e)
      {
        // try to calculate real coordinates to reference coordinates
        DRT::Element* actele = discret_->lColElement(discret_->ElementColMap()->LID(
            discret_->lRowNode(lid_closestnode)->Elements()[e]->Id()));

        int ndof = actele->NumNode();

        Epetra_SerialDenseMatrix elematrix1(ndof, ndof, false);
        Epetra_SerialDenseMatrix elematrix2(ndof, ndof, false);
        Epetra_SerialDenseVector elevector1(ndof);
        Epetra_SerialDenseVector elevector2(ndof);
        Epetra_SerialDenseVector elevector3(ndof);

        DRT::Element::LocationArray la(discret_->NumDofSets());
        actele->LocationVector(*discret_, la, false);

        Teuchos::ParameterList p;
        p.set<int>("action", SCATRA::transform_real_to_reference_point);
        double pointarr[dim];
        for (int d = 0; d < dim; ++d) pointarr[d] = point[d];
        p.set<double*>("point", pointarr);

        actele->Evaluate(
            p, *discret_, la, elematrix1, elematrix2, elevector1, elevector2, elevector3);
        for (int d = 0; d < dim; ++d) pointarr[d] = p.get<double*>("point")[d];

        if (p.get<bool>("inside"))  // third step: evaluate in the element
        {
          anywhere_inside = 1;
          p.set<int>("action", SCATRA::evaluate_field_in_point);
          p.set<int>("numscal", numscal);
          actele->Evaluate(
              p, *discret_, la, elematrix1, elematrix2, elevector1, elevector2, elevector3);
          lvalue = p.get<double>("value");

          if (evalreac)  // multiply with reaction coefficient
          {
            lvalue *= actele->Material()->Parameter()->GetParameter(1, actele->Id());
          }
        }
      }
    }
    int foundit = 0;
    discret_->Comm().MaxAll(&anywhere_inside, &foundit, 1);
    if (foundit < 1)
    {
      lvalue = 0.0;
      // std::cout<<"warning: you called PointValue on a point which is not inside of the Scatra
      // domain"<<std::endl;
    }
  }
  double gvalue = 0.0;
  discret_->Comm().SumAll(&lvalue, &gvalue, 1);

  value = gvalue;

  return;
}

void SCATRA::ScaTraTimIntImpl::GetPointsPhiValues(const std::vector<std::vector<double>>& points,
    std::vector<double>& values, bool evalreac, unsigned int numscal)
{
  if (int(numscal) > this->NumScal())
    dserror("you requested the point values for scalar %d but there is only %d scalar fields",
        numscal, this->NumScal());

  for (unsigned int p = 0; p < points.size(); ++p)
    GetPointPhiValue(points[p], values[p], evalreac, numscal);

  return;
}

/*----------------------------------------------------------------------------------------------------*
 | evaluate macro-micro coupling on micro scale in multi-scale scalar transport problems   fang
 01/16 |
 *----------------------------------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::EvaluateMacroMicroCoupling()
{
  // extract multi-scale coupling conditions
  std::vector<Teuchos::RCP<DRT::Condition>> conditions;
  discret_->GetCondition("ScatraMultiScaleCoupling", conditions);

  // loop over conditions
  for (unsigned icond = 0; icond < conditions.size(); ++icond)
  {
    // extract nodal cloud
    const std::vector<int>* const nodeids = conditions[icond]->Nodes();
    if (nodeids == NULL) dserror("Multi-scale coupling condition does not have nodal cloud!");

    // loop over all nodes in nodal cloud
    for (unsigned inode = 0; inode < (*nodeids).size(); ++inode)
    {
      // process row nodes only
      if (discret_->NodeRowMap()->MyGID((*nodeids)[inode]))
      {
        // extract node
        DRT::Node* node = discret_->gNode((*nodeids)[inode]);
        if (node == NULL)
          dserror("Cannot extract node with global ID %d from micro-scale discretization!",
              (*nodeids)[inode]);

        // safety check
        if (node->NumElement() != 1)
          dserror("Number of 1D elements adjacent to the boundary node must be 1!");

        // compute domain integration factor
        double fac(1.);
        if (DRT::INPUT::IntegralValue<bool>(*params_, "SPHERICALCOORDS"))
          fac *= *node->X() * *node->X();

        // extract degrees of freedom from node
        const std::vector<int> dofs = discret_->Dof(0, node);

        // loop over all degrees of freedom
        for (unsigned idof = 0; idof < dofs.size(); ++idof)
        {
          // extract global and local IDs of degree of freedom
          const int gid = dofs[idof];
          const int lid = discret_->DofRowMap()->LID(gid);
          if (lid < 0) dserror("Cannot extract degree of freedom with global ID %d!", gid);

          // compute matrix and vector contributions according to kinetic model for current
          // macro-micro coupling condition
          switch (conditions[icond]->GetInt("kinetic model"))
          {
            case INPAR::S2I::kinetics_constperm:
            {
              // access real vector of constant permeabilities
              const std::vector<double>* permeabilities =
                  conditions[icond]->GetMutable<std::vector<double>>("permeabilities");
              if (permeabilities == NULL)
                dserror("Cannot access vector of permeabilities for macro-micro coupling!");
              if (permeabilities->size() != (unsigned)NumScal())
                dserror("Number of permeabilities does not match number of scalars!");

              // compute and store micro-scale coupling flux
              q_ = (*permeabilities)[0] * ((*phinp_)[lid] - phinp_macro_[0]);

              // compute and store derivative of micro-scale coupling flux w.r.t. macro-scale state
              // variable
              dq_dphi_[0] = -(*permeabilities)[0];

              // assemble contribution from macro-micro coupling into global residual vector
              (*residual_)[lid] -=
                  DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(discret_->Name())
                      ->TimeFacRhs() *
                  q_ * fac;

              // assemble contribution from macro-micro coupling into global system matrix
              sysmat_->Assemble(
                  DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(discret_->Name())->TimeFac() *
                      (*permeabilities)[0] * fac,
                  gid, gid);

              break;
            }

            case INPAR::S2I::kinetics_butlervolmer:
            {
              // access material of electrode
              Teuchos::RCP<const MAT::Electrode> matelectrode =
                  Teuchos::rcp_dynamic_cast<const MAT::Electrode>(node->Elements()[0]->Material());
              if (matelectrode == Teuchos::null)
                dserror("Invalid electrode material for multi-scale coupling!");

              // access input parameters associated with current condition
              const int nume = conditions[icond]->GetInt("e-");
              if (nume != 1)
                dserror(
                    "Invalid number of electrons involved in charge transfer at "
                    "electrode-electrolyte interface!");
              const std::vector<int>* stoichiometries =
                  conditions[icond]->GetMutable<std::vector<int>>("stoichiometries");
              if (stoichiometries == NULL)
                dserror(
                    "Cannot access vector of stoichiometric coefficients for multi-scale "
                    "coupling!");
              if (stoichiometries->size() != 1)
                dserror("Number of stoichiometric coefficients does not match number of scalars!");
              if ((*stoichiometries)[0] != -1) dserror("Invalid stoichiometric coefficient!");
              const double faraday =
                  DRT::Problem::Instance(0)->ELCHControlParams().get<double>("FARADAY_CONSTANT");
              const double gasconstant =
                  DRT::Problem::Instance(0)->ELCHControlParams().get<double>("GAS_CONSTANT");
              const double frt =
                  faraday /
                  (gasconstant *
                      (DRT::Problem::Instance(0)->ELCHControlParams().get<double>("TEMPERATURE")));
              const double alphaa =
                  conditions[icond]->GetDouble("alpha_a");  // anodic transfer coefficient
              const double alphac =
                  conditions[icond]->GetDouble("alpha_c");  // cathodic transfer coefficient
              const double kr =
                  conditions[icond]->GetDouble("k_r");  // rate constant of charge transfer reaction
              if (kr < 0.) dserror("Charge transfer constant k_r is negative!");

              // extract saturation value of intercalated lithium concentration from electrode
              // material
              const double cmax = matelectrode->CMax();
              if (cmax < 1.e-12)
                dserror(
                    "Saturation value c_max of intercalated lithium concentration is too small!");

              // extract electrode-side and electrolyte-side concentration values at multi-scale
              // coupling point
              const double conc_ed = (*phinp_)[lid];
              const double conc_el = phinp_macro_[0];

              // evaluate overall integration factors
              const double timefacfac =
                  DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(discret_->Name())->TimeFac() *
                  fac;
              const double timefacrhsfac =
                  DRT::ELEMENTS::ScaTraEleParameterTimInt::Instance(discret_->Name())
                      ->TimeFacRhs() *
                  fac;
              if (timefacfac < 0. or timefacrhsfac < 0.) dserror("Integration factor is negative!");

              // equilibrium electric potential difference and its derivative w.r.t. concentration
              // at electrode surface
              const double epd = matelectrode->ComputeOpenCircuitPotential(conc_ed, faraday, frt);
              const double epdderiv =
                  matelectrode->ComputeFirstDerivOpenCircuitPotential(conc_ed, faraday, frt);

              // electrode-electrolyte overpotential at multi-scale coupling point
              const double eta = phinp_macro_[2] - phinp_macro_[1] - epd;

              // Butler-Volmer exchange mass flux density
              const double j0 =
                  kr * pow(conc_el, alphaa) * pow(cmax - conc_ed, alphaa) * pow(conc_ed, alphac);

              // exponential Butler-Volmer terms
              const double expterm1 = exp(alphaa * frt * eta);
              const double expterm2 = exp(-alphac * frt * eta);
              const double expterm = expterm1 - expterm2;

              // safety check
              if (abs(expterm) > 1.e5)
                dserror(
                    "Overflow of exponential term in Butler-Volmer formulation detected! Value: "
                    "%lf",
                    expterm);

              // core residual term associated with Butler-Volmer mass flux density
              q_ = j0 * expterm;

              // core linearizations associated with Butler-Volmer mass flux density
              const double dj_dc_ed =
                  kr * pow(conc_el, alphaa) * pow(cmax - conc_ed, alphaa - 1.) *
                      pow(conc_ed, alphac - 1.) * (-alphaa * conc_ed + alphac * (cmax - conc_ed)) *
                      expterm +
                  j0 * (-alphaa * frt * epdderiv * expterm1 - alphac * frt * epdderiv * expterm2);
              dq_dphi_[0] = j0 * alphaa / conc_el * expterm;  // dj_dc_el
              dq_dphi_[1] =
                  -j0 * (alphaa * frt * expterm1 + alphac * frt * expterm2);  // dj_dpot_el
              dq_dphi_[2] = -dq_dphi_[1];                                     // dj_dpot_ed

              // assemble contribution from macro-micro coupling into global residual vector
              (*residual_)[lid] -= timefacrhsfac * q_;

              // assemble contribution from macro-micro coupling into micro global system matrix
              sysmat_->Assemble(timefacfac * dj_dc_ed, gid, gid);

              break;
            }

            default:
            {
              dserror("Kinetic model for macro-micro coupling not yet implemented!");
              break;
            }
          }
        }
      }
    }
  }

  return;
}  // SCATRA::ScaTraTimIntImpl::EvaluateMacroMicroCoupling


/*-----------------------------------------------------------------------------*
 |  check if class is initialized                                  rauch 09/16 |
 *-----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::CheckIsInit() const
{
  if (not IsInit()) dserror("ScaTraTimIntImpl is not initialized. Call Init() first.");
}


/*-----------------------------------------------------------------------------*
 |  check if class is set up                                       rauch 09/16 |
 *-----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::CheckIsSetup() const
{
  if (not IsSetup()) dserror("ScaTraTimIntImpl is not set up. Call Setup() first.");
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
bool SCATRA::ScaTraTimIntImpl::IsS2IMeshtying() const
{
  auto problem = DRT::Problem::Instance();
  bool IsS2IMeshtying(false);
  // decide depending on the problem type
  switch (problem->GetProblemType())
  {
    case prb_elch:
    case prb_scatra:
    case prb_sti:
    {
      if (s2icoupling_) IsS2IMeshtying = true;
      break;
    }
    case prb_ssi:
    {
      // get structure discretization
      auto structdis = problem->GetDis("structure");

      // get ssi meshtying conditions
      std::vector<DRT::Condition*> ssiconditions;
      structdis->GetCondition("SSIInterfaceMeshtying", ssiconditions);

      // do mesh tying if there is at least one mesh tying condition
      if (!ssiconditions.empty()) IsS2IMeshtying = true;
      break;
    }
    default:
    {
      // do nothing
      break;
    }
  }

  return IsS2IMeshtying;
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetupMatrixBlockMaps()
{
  if (matrixtype_ == INPAR::SCATRA::MatrixType::block_condition or
      matrixtype_ == INPAR::SCATRA::MatrixType::block_condition_dof)
  {
    // extract domain partitioning conditions from discretization
    std::vector<Teuchos::RCP<DRT::Condition>> partitioningconditions;
    discret_->GetCondition("ScatraPartitioning", partitioningconditions);

    // safety check
    if (partitioningconditions.empty())
      dserror(
          "For block preconditioning based on domain partitioning, at least one associated "
          "condition needs to be specified in the input file!");

    // build maps associated with blocks of global system matrix
    std::vector<Teuchos::RCP<const Epetra_Map>> blockmaps;
    BuildBlockMaps(partitioningconditions, blockmaps);

    // initialize full map extractor associated with blocks of global system matrix
    blockmaps_ = Teuchos::rcp(new LINALG::MultiMapExtractor(*(discret_->DofRowMap()), blockmaps));
    // safety check
    blockmaps_->CheckForValidMapExtractor();
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::BuildBlockMaps(
    const std::vector<Teuchos::RCP<DRT::Condition>>& partitioningconditions,
    std::vector<Teuchos::RCP<const Epetra_Map>>& blockmaps) const
{
  if (matrixtype_ == INPAR::SCATRA::MatrixType::block_condition)
  {
    // extract number of domain partitioning conditions
    const unsigned ncond = partitioningconditions.size();

    // prepare vector for maps to be built
    blockmaps.resize(ncond, Teuchos::null);

    // loop over all domain partitioning conditions
    for (unsigned icond = 0; icond < ncond; ++icond)
    {
      // initialize set for dof IDs associated with current partitioning condition
      std::set<int> dofids;

      // extract nodes associated with current domain partitioning condition
      const std::vector<int>* nodegids = partitioningconditions[icond]->Nodes();

      // loop over all nodes associated with current domain partitioning condition
      for (int nodegid : *nodegids)
      {
        // extract global ID of current node
        // consider current node only if node is owned by current processor
        // need to make sure that node is stored on current processor, otherwise cannot resolve
        // "->Owner()"
        if (discret_->HaveGlobalNode(nodegid) and
            discret_->gNode(nodegid)->Owner() == discret_->Comm().MyPID())
        {
          // add dof IDs associated with current node to corresponding set
          const std::vector<int> nodedofs = discret_->Dof(0, discret_->gNode(nodegid));
          std::copy(nodedofs.begin(), nodedofs.end(), std::inserter(dofids, dofids.end()));
        }
      }

      // transform set for dof IDs into vector and then into Epetra map
      int nummyelements(0);
      int* myglobalelements(nullptr);
      std::vector<int> dofidvec;
      if (!dofids.empty())
      {
        dofidvec.reserve(dofids.size());
        dofidvec.assign(dofids.begin(), dofids.end());
        nummyelements = dofidvec.size();
        myglobalelements = &(dofidvec[0]);
      }
      blockmaps[icond] = Teuchos::rcp(new Epetra_Map(-1, nummyelements, myglobalelements,
          discret_->DofRowMap()->IndexBase(), discret_->DofRowMap()->Comm()));
    }
  }
  // safety check
  else
    dserror("Invalid type of global system matrix!");
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::PostSetupMatrixBlockMaps()
{
  // matrix block map extractor equals interface map extractor in this case
  if (matrixtype_ == INPAR::SCATRA::MatrixType::block_geometry)
    blockmaps_ = strategy_->InterfaceMaps();

  // now build the null spaces
  BuildBlockNullSpaces();

  // in case of an extended solver for scatra-scatra interface meshtying including interface growth
  // we need to equip it with the null space information generated above
  if (IsS2IMeshtying()) strategy_->EquipExtendedSolverWithNullSpaceInfo();
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::BuildBlockNullSpaces() const
{
  // loop over blocks of global system matrix
  for (int iblock = 0; iblock < BlockMaps().NumMaps(); ++iblock)
  {
    // store number of current block as string, starting from 1
    std::ostringstream iblockstr;
    iblockstr << iblock + 1;

    // equip smoother for current matrix block with empty parameter sublists to trigger null space
    // computation
    Teuchos::ParameterList& blocksmootherparams =
        Solver()->Params().sublist("Inverse" + iblockstr.str());
    blocksmootherparams.sublist("Aztec Parameters");
    blocksmootherparams.sublist("MueLu Parameters");

    // equip smoother for current matrix block with null space associated with all degrees of
    // freedom on discretization
    discret_->ComputeNullSpaceIfNecessary(blocksmootherparams);

    // reduce full null space to match degrees of freedom associated with current matrix block
    LINALG::Solver::FixMLNullspace("Block " + iblockstr.str(), *discret_->DofRowMap(),
        *BlockMaps().Map(iblock), blocksmootherparams);
  }
}

/*-----------------------------------------------------------------------------*
 *-----------------------------------------------------------------------------*/
void SCATRA::ScaTraTimIntImpl::SetupMatrixBlockMapsAndMeshtying()
{
  switch (MatrixType())
  {
    // case INPAR::SCATRA::MatrixType::undefined:
    case INPAR::SCATRA::MatrixType::sparse:
    {
      // only setup the meshtying in this case, as matrix has no block structure
      strategy_->SetupMeshtying();

      break;
    }
    case INPAR::SCATRA::MatrixType::block_condition:
    case INPAR::SCATRA::MatrixType::block_condition_dof:
    case INPAR::SCATRA::MatrixType::block_geometry:
    {
      // safety check
      if (!Solver()->Params().isSublist("AMGnxn Parameters"))
        dserror("Global system matrix with block structure requires AMGnxn block preconditioner!");

      // setup the matrix block maps
      SetupMatrixBlockMaps();

      // setup the meshtying
      strategy_->SetupMeshtying();

      // do some post setup matrix block map operations after the call to SetupMeshtying, as they
      // rely on the fact that the interface maps have already been built
      PostSetupMatrixBlockMaps();

      break;
    }
    default:
    {
      dserror("ScaTra Matrixtype %i not recognised", static_cast<int>(MatrixType()));
      break;
    }
  }
}