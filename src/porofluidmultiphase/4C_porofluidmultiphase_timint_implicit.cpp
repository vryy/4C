/*----------------------------------------------------------------------*/
/*! \file
 \brief base class of implicit integration schemes for porous multiphase
        flow problems

   \level 3

 *----------------------------------------------------------------------*/


#include "4C_porofluidmultiphase_timint_implicit.hpp"

#include "4C_fem_general_assemblestrategy.hpp"
#include "4C_fem_general_l2_projection.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_io_gmsh.hpp"
#include "4C_linalg_utils_sparse_algebra_assemble.hpp"
#include "4C_linalg_utils_sparse_algebra_create.hpp"
#include "4C_linalg_utils_sparse_algebra_print.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_mat_fluidporo_multiphase.hpp"
#include "4C_porofluidmultiphase_ele.hpp"
#include "4C_porofluidmultiphase_ele_action.hpp"
#include "4C_porofluidmultiphase_meshtying_strategy_artery.hpp"
#include "4C_porofluidmultiphase_meshtying_strategy_std.hpp"
#include "4C_porofluidmultiphase_resulttest.hpp"
#include "4C_porofluidmultiphase_utils.hpp"
#include "4C_utils_function.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*==========================================================================*/
// Constructors and destructors and related methods
/*==========================================================================*/

/*----------------------------------------------------------------------*
 | constructor                                     (public) vuong 08/16 |
 *----------------------------------------------------------------------*/
POROFLUIDMULTIPHASE::TimIntImpl::TimIntImpl(Teuchos::RCP<Core::FE::Discretization> actdis,
    const int linsolvernumber, const Teuchos::ParameterList& probparams,
    const Teuchos::ParameterList& poroparams,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output)
    :  // call constructor for "nontrivial" objects
      solver_(Teuchos::null),
      linsolvernumber_(linsolvernumber),
      params_(probparams),
      poroparams_(poroparams),
      myrank_(actdis->Comm().MyPID()),
      nsd_(Global::Problem::Instance()->NDim()),
      isale_(false),
      skipinitder_(Core::UTILS::IntegralValue<int>(poroparams_, "SKIPINITDER")),
      output_satpress_(Core::UTILS::IntegralValue<int>(poroparams_, "OUTPUT_SATANDPRESS")),
      output_solidpress_(Core::UTILS::IntegralValue<int>(poroparams_, "OUTPUT_SOLIDPRESS")),
      output_porosity_(Core::UTILS::IntegralValue<int>(poroparams_, "OUTPUT_POROSITY")),
      output_phase_velocities_(
          Core::UTILS::IntegralValue<int>(poroparams_, "OUTPUT_PHASE_VELOCITIES")),
      output_bloodvesselvolfrac_(Core::UTILS::IntegralValue<int>(
          poroparams_.sublist("ARTERY COUPLING"), "OUTPUT_BLOODVESSELVOLFRAC")),
      stab_biot_(Core::UTILS::IntegralValue<int>(poroparams_, "STAB_BIOT")),
      domainint_funct_(std::vector<int>()),
      num_domainint_funct_(0),
      calcerr_(Core::UTILS::IntegralValue<Inpar::POROFLUIDMULTIPHASE::CalcError>(
          poroparams_, "CALCERROR")),
      fluxrecon_(Core::UTILS::IntegralValue<Inpar::POROFLUIDMULTIPHASE::FluxReconstructionMethod>(
          poroparams_, "FLUX_PROJ_METHOD")),
      fluxreconsolvernum_(poroparams_.get<int>("FLUX_PROJ_SOLVER")),
      divcontype_(Core::UTILS::IntegralValue<Inpar::POROFLUIDMULTIPHASE::DivContAct>(
          poroparams_, "DIVERCONT")),
      fdcheck_(
          Core::UTILS::IntegralValue<Inpar::POROFLUIDMULTIPHASE::FdCheck>(poroparams_, "FDCHECK")),
      fdcheckeps_(poroparams_.get<double>("FDCHECKEPS")),
      fdchecktol_(poroparams_.get<double>("FDCHECKTOL")),
      stab_biot_scaling_(poroparams_.get<double>("STAB_BIOT_SCALING")),
      time_(0.0),
      maxtime_(params_.get<double>("MAXTIME")),
      step_(0),
      stepmax_(params_.get<int>("NUMSTEP")),
      dt_(params_.get<double>("TIMESTEP")),
      dtele_(0.0),
      dtsolve_(0.0),
      iternum_(0),
      itemax_(poroparams_.get<int>("ITEMAX")),
      upres_(params_.get<int>("RESULTSEVRY")),
      uprestart_(params_.get<int>("RESTARTEVRY")),
      vectornormfres_(Core::UTILS::IntegralValue<Inpar::POROFLUIDMULTIPHASE::VectorNorm>(
          poroparams_, "VECTORNORM_RESF")),
      vectornorminc_(Core::UTILS::IntegralValue<Inpar::POROFLUIDMULTIPHASE::VectorNorm>(
          poroparams_, "VECTORNORM_INC")),
      ittolres_(poroparams_.get<double>("TOLRES")),
      ittolinc_(poroparams_.get<double>("TOLINC")),
      artery_coupling_active_(Core::UTILS::IntegralValue<int>(params_, "ARTERY_COUPLING")),
      // Initialization of degrees of freedom variables
      phin_(Teuchos::null),
      phinp_(Teuchos::null),
      phidtn_(Teuchos::null),
      phidtnp_(Teuchos::null),
      hist_(Teuchos::null),
      pressure_(Teuchos::null),
      saturation_(Teuchos::null),
      solidpressure_(Teuchos::null),
      valid_volfracpress_dofs_(Teuchos::null),
      valid_volfracspec_dofs_(Teuchos::null),
      flux_(Teuchos::null),
      nds_disp_(-1),
      nds_vel_(-1),
      nds_solidpressure_(-1),
      nds_scatra_(-1),
      discret_(actdis),
      output_(output),
      sysmat_(Teuchos::null),
      zeros_(Teuchos::null),
      dbcmaps_(Teuchos::null),
      dbcmaps_with_volfracpress_(Teuchos::null),
      dbcmaps_starting_condition_(Teuchos::null),
      neumann_loads_(Teuchos::null),
      residual_(Teuchos::null),
      trueresidual_(Teuchos::null),
      increment_(Teuchos::null),
      starting_dbc_time_end_(poroparams_.get<double>("STARTING_DBC_TIME_END")),
      starting_dbc_onoff_(std::vector<bool>()),
      starting_dbc_funct_(std::vector<int>())
{
  return;
}


/*------------------------------------------------------------------------*
 | initialize time integration                                vuong 08/16 |
 *------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::Init(bool isale, int nds_disp, int nds_vel,
    int nds_solidpressure, int nds_scalar, const std::map<int, std::set<int>>* nearbyelepairs)
{
  // set flags
  isale_ = isale;
  nds_disp_ = nds_disp;
  nds_vel_ = nds_vel;
  nds_solidpressure_ = nds_solidpressure;
  nds_scatra_ = nds_scalar;

  // make sure the values make sense
  // -1 is the default value, meaning that there is no coupling
  if (nds_disp_ != -1)
    if (nds_disp_ < 0 or nds_disp_ > discret_->NumDofSets() - 1)
      FOUR_C_THROW("invalid number of dofset for mesh displacements!");

  // make sure the values make sense
  // -1 is the default value, meaning that there is no coupling
  if (nds_vel_ != -1)
    if (nds_vel_ < 0 or nds_vel_ > discret_->NumDofSets() - 1)
      FOUR_C_THROW("invalid number of dofset for mesh velocities!");

  // make sure the values make sense
  // there has to be a valid number for the solid pressure in all cases
  if (nds_solidpressure_ < 0 or nds_solidpressure_ > discret_->NumDofSets() - 1)
    FOUR_C_THROW("invalid number of dofset for solid pressure!");

  // ensure that degrees of freedom in the discretization have been set
  if ((not discret_->Filled()) or (not discret_->HaveDofs())) discret_->fill_complete();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices: local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  // -------------------------------------------------------------------
  // create empty system matrix (27 adjacent nodes as 'good' guess)
  // -------------------------------------------------------------------
  sysmat_ =
      Teuchos::rcp(new Core::LinAlg::SparseMatrix(*(discret_->dof_row_map()), 27, false, true));

  // -------------------------------------------------------------------
  // create vectors containing problem variables
  // -------------------------------------------------------------------
  // solutions at time n+1
  phinp_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  // solutions at time n
  phin_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  // time derivative of solutions at time n
  phidtn_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  // time derivative of solutions at time n+1
  phidtnp_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // history vector
  hist_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  // valid (physically meaningful) volume fraction dofs
  valid_volfracpress_dofs_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  valid_volfracspec_dofs_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  if (output_satpress_)
  {
    // pressure at time n+1
    pressure_ = Core::LinAlg::CreateVector(*dofrowmap, true);
    // saturation at time n+1
    saturation_ = Core::LinAlg::CreateVector(*dofrowmap, true);
  }
  // solid pressure at time n+1
  if (output_solidpress_)
    solidpressure_ = Core::LinAlg::CreateVector(*discret_->dof_row_map(nds_solidpressure_), true);
  // porosity at time n+1 (lives on same dofset as solid pressure)
  if (output_porosity_)
    porosity_ = Core::LinAlg::CreateVector(*discret_->dof_row_map(nds_solidpressure_), true);

  if (output_phase_velocities_)
  {
    const int num_poro_dof = discret_->NumDof(0, discret_->lRowNode(0));
    const int num_rows = num_poro_dof * nsd_;
    phase_velocities_ = Core::LinAlg::CreateMultiVector(*discret_->ElementRowMap(), num_rows, true);
  }

  // -------------------------------------------------------------------
  // create vectors associated to boundary conditions
  // -------------------------------------------------------------------
  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  zeros_ = Core::LinAlg::CreateVector(*dofrowmap, true);

  int stream;
  std::istringstream stream_dbc_onoff(
      Teuchos::getNumericStringParameter(poroparams_, "STARTING_DBC_ONOFF"));
  while (stream_dbc_onoff >> stream) starting_dbc_onoff_.push_back(static_cast<bool>(stream));

  std::istringstream stream_dbc_funct(
      Teuchos::getNumericStringParameter(poroparams_, "STARTING_DBC_FUNCT"));
  while (stream_dbc_funct >> stream) starting_dbc_funct_.push_back(static_cast<int>(stream));

  // object holds maps/subsets for DOFs subjected to Dirichlet BCs and otherwise
  dbcmaps_ = Teuchos::rcp(new Core::LinAlg::MapExtractor());
  dbcmaps_with_volfracpress_ = Teuchos::rcp(new Core::LinAlg::MapExtractor());
  dbcmaps_starting_condition_ = Teuchos::rcp(new Core::LinAlg::MapExtractor());
  {
    Teuchos::ParameterList eleparams;
    // other parameters needed by the elements
    eleparams.set("total time", time_);
    eleparams.set<const Core::UTILS::FunctionManager*>(
        "function_manager", &Global::Problem::Instance()->FunctionManager());
    discret_->evaluate_dirichlet(
        eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_);
    discret_->evaluate_dirichlet(
        eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null, dbcmaps_with_volfracpress_);
    discret_->evaluate_dirichlet(eleparams, zeros_, Teuchos::null, Teuchos::null, Teuchos::null,
        dbcmaps_starting_condition_);
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

  // -------------------------------------------------------------------
  // set initial field
  // -------------------------------------------------------------------
  SetInitialField(Core::UTILS::IntegralValue<Inpar::POROFLUIDMULTIPHASE::InitialField>(
                      poroparams_, "INITIALFIELD"),
      poroparams_.get<int>("INITFUNCNO"));

  // -------------------------------------------------------------------
  // domain integration functions for output
  // -------------------------------------------------------------------
  int word1;
  std::istringstream coupled_art_dof_stream(
      Teuchos::getNumericStringParameter(poroparams_, "DOMAININT_FUNCT"));
  while (coupled_art_dof_stream >> word1) domainint_funct_.push_back((int)(word1));
  // no domain integration function selected by user
  if (domainint_funct_.size() == 1 and domainint_funct_[0] < 0) domainint_funct_.resize(0);
  num_domainint_funct_ = domainint_funct_.size();

  // the values of the integrals
  domain_integrals_ = Teuchos::rcp(new Core::LinAlg::SerialDenseVector(num_domainint_funct_));

  // -------------------------------------------------------------------
  // set element parameters
  // -------------------------------------------------------------------
  set_element_general_parameters();

  // -------------------------------------------------------------------
  // build mesh tying strategy
  // -------------------------------------------------------------------
  if (artery_coupling_active_)
    strategy_ =
        Teuchos::rcp(new POROFLUIDMULTIPHASE::MeshtyingStrategyArtery(this, params_, poroparams_));
  else
    strategy_ =
        Teuchos::rcp(new POROFLUIDMULTIPHASE::MeshtyingStrategyStd(this, params_, poroparams_));
  // check if initial fields match
  strategy_->CheckInitialFields(phinp_);
  // set the nearby ele pairs
  strategy_->SetNearbyElePairs(nearbyelepairs);
  // setup the strategy
  strategy_->Setup();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  solver_ = Teuchos::rcp(
      new Core::LinAlg::Solver(Global::Problem::Instance()->SolverParams(linsolvernumber_),
          discret_->Comm(), Global::Problem::Instance()->solver_params_callback(),
          Core::UTILS::IntegralValue<Core::IO::Verbositylevel>(
              Global::Problem::Instance()->IOParams(), "VERBOSITY")));
  strategy_->initialize_linear_solver(solver_);

  return;
}  // TimIntImpl::Init()



/*========================================================================*/
//! set element parameters
/*========================================================================*/

/*----------------------------------------------------------------------*
 | set all general parameters for element                   vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::set_element_general_parameters() const
{
  Teuchos::ParameterList eleparams;

  eleparams.set<int>("action", POROFLUIDMULTIPHASE::set_general_parameter);

  eleparams.set<bool>("isale", isale_);
  eleparams.set<int>("nds_disp", nds_disp_);
  eleparams.set<int>("nds_vel", nds_vel_);
  eleparams.set<int>("nds_solidpressure", nds_solidpressure_);
  eleparams.set<int>("nds_scalar", nds_scatra_);
  eleparams.set<bool>("stab_biot", stab_biot_);

  eleparams.set<bool>("using generalized-alpha time integration", false);
  eleparams.set<bool>("using stationary formulation", false);
  eleparams.set<double>("alpha_F", 1.0);

  eleparams.set<int>("num_domainint_funct", num_domainint_funct_);
  for (int ifunct = 0; ifunct < num_domainint_funct_; ifunct++)
    eleparams.set<int>("domainint_funct_" + std::to_string(ifunct), domainint_funct_[ifunct]);

  // call standard loop over elements
  discret_->evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);

  return;
}


/*==========================================================================*/
// general framework
/*==========================================================================*/

/*--- set, prepare, and predict --------------------------------------------*/

/*----------------------------------------------------------------------*
 | prepare time loop                                        vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::prepare_time_loop()
{
  // compute pressure and saturations
  reconstruct_pressures_and_saturations();

  // compute velocities
  ReconstructFlux();

  // provide information about initial field (do not do for restarts!)
  if (step_ == 0)
  {
    // write out initial state
    Output();

    // compute error for problems with analytical solution (initial field!)
    evaluate_error_compared_to_analytical_sol();
  }

  // do the same also for meshtying
  strategy_->prepare_time_loop();

  return;
}  // POROFLUIDMULTIPHASE::TimIntImpl::prepare_time_loop


/*----------------------------------------------------------------------*
 | setup the variables to do a new time step       (public) vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::prepare_time_step()
{
  // time measurement: prepare time step
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:    + prepare time step");

  // -------------------------------------------------------------------
  //                       initialization
  // -------------------------------------------------------------------
  if (step_ == 0) prepare_first_time_step();

  // -------------------------------------------------------------------
  //              set time dependent parameters
  // -------------------------------------------------------------------
  // note the order of the following three functions is important
  increment_time_and_step();

  // -------------------------------------------------------------------
  // set part of the rhs vector belonging to the old timestep
  // -------------------------------------------------------------------
  set_old_part_of_righthandside();
  // reset every parameter that potentially changes for every time step
  set_element_time_step_parameter();

  // -------------------------------------------------------------------
  //         evaluate Dirichlet and Neumann boundary conditions
  // -------------------------------------------------------------------
  // TODO: Dirichlet auch im Fall von genalpha prenp
  // Neumann(n + alpha_f)
  apply_dirichlet_bc(time_, phinp_, Teuchos::null);
  apply_neumann_bc(neumann_loads_);

  // volume fraction pressure specific stuff
  evaluate_valid_volume_frac_press_and_spec();
  apply_additional_dbc_for_vol_frac_press();

  if (time_ <= starting_dbc_time_end_)
  {
    apply_starting_dbc();
  }

  // do the same also for meshtying
  strategy_->prepare_time_step();

  return;
}  // TimIntImpl::prepare_time_step


/*------------------------------------------------------------------------------*
 | initialization procedure prior to evaluation of first time step  vuong 08/16 |
 *------------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::prepare_first_time_step()
{
  if (not skipinitder_)
  {
    if (nds_vel_ != -1 || !isale_)  // if some velocity field has been set
    {
      // TODO: Restructure enforcement of Dirichlet boundary conditions on phin_
      apply_dirichlet_bc(time_, phin_, Teuchos::null);
      calc_initial_time_derivative();
    }
    // if initial velocity field has not been set here, the initial time derivative of phi will be
    // calculated wrongly for some time integration schemes
    else
      FOUR_C_THROW("Initial velocity field has not been set!");
  }
  return;
}  // POROFLUIDMULTIPHASE::TimIntImpl::prepare_first_time_step


/*----------------------------------------------------------------------*
 | contains the time loop                                   vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::TimeLoop()
{
  // time measurement: time loop
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:  + time loop");

  // prepare time loop
  prepare_time_loop();

  while ((step_ < stepmax_) and ((time_ + 1e-12) < maxtime_))
  {
    // -------------------------------------------------------------------
    // prepare time step
    // -------------------------------------------------------------------
    prepare_time_step();

    // -------------------------------------------------------------------
    //                  solve nonlinear / linear equation
    // -------------------------------------------------------------------
    Solve();

    // -------------------------------------------------------------------
    //                         update solution
    //        current solution becomes old solution of next timestep
    // -------------------------------------------------------------------
    Update();

    // -------------------------------------------------------------------
    // evaluate error for problems with analytical solution
    // -------------------------------------------------------------------
    if (calcerr_ != Inpar::POROFLUIDMULTIPHASE::calcerror_no)
      evaluate_error_compared_to_analytical_sol();

    // -------------------------------------------------------------------
    //                         output of solution
    // -------------------------------------------------------------------
    Output();

  }  // while

  // print the results of time measurements
  Teuchos::TimeMonitor::summarize();

  return;
}  // TimIntImpl::TimeLoop


/*----------------------------------------------------------------------*
 | contains the call of linear/nonlinear solver             vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::Solve()
{
  // -----------------------------------------------------------------
  //                    always solve nonlinear equation
  // -----------------------------------------------------------------
  nonlinear_solve();

  // reconstruct pressures and saturations
  reconstruct_pressures_and_saturations();

  // reconstruct velocities
  ReconstructFlux();

  return;
}

/*----------------------------------------------------------------------*
 | contains the call of linear/nonlinear solver             vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::Update() { strategy_->Update(); }


/*----------------------------------------------------------------------*
 | apply moving mesh data                                   vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::ApplyMeshMovement(Teuchos::RCP<const Epetra_Vector> dispnp)
{
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE: apply mesh movement");

  if (nds_disp_ == -1)
    FOUR_C_THROW(
        "Dof set number of displacment related dofs"
        " has not been set!");

  // check existence of displacement vector
  if (dispnp == Teuchos::null) FOUR_C_THROW("Got null pointer for displacements!");

  // provide POROFLUIDMULTIPHASE discretization with displacement field
  set_state(nds_disp_, "dispnp", dispnp);

  // apply mesh movement also on the strategy
  strategy_->ApplyMeshMovement();

  return;
}  // TimIntImpl::ApplyMeshMovement


/*----------------------------------------------------------------------*
 |  print information about current time step to screen     vuong 08/16 |
 *----------------------------------------------------------------------*/
inline void POROFLUIDMULTIPHASE::TimIntImpl::print_time_step_info()
{
  if (myrank_ == 0)
    printf("TIME: %11.4E/%11.4E  DT = %11.4E  Stationary  STEP = %4d/%4d \n", time_, maxtime_, dt_,
        step_, stepmax_);
}  // POROFLUIDMULTIPHASE::TimIntImpl::print_time_step_info


/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                       vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::Output()
{
  // time measurement: output of solution
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:    + output of solution");

  // solution output and potentially restart data and/or flux data
  if (do_output())
  {
    // do the same for the strategy
    strategy_->Output();

    // step number and time (only after that data output is possible)
    output_->new_step(step_, time_);

    // write domain decomposition for visualization (only once at step 0!)
    if (step_ == 0)
    {
      output_->write_element_data(true);
      // write output of blood vessel volume fraction
      if (output_bloodvesselvolfrac_)
        output_->write_vector("bloodvesselvolfrac", strategy_->blood_vessel_volume_fraction(),
            Core::IO::elementvector);
    }

    // reconstruct porosity for output; porosity is only needed for output and does not have to be
    // transferred between fields
    if (output_porosity_) reconstruct_porosity();

    if (output_phase_velocities_) calculate_phase_velocities();

    // evaluate domain integrals
    if (num_domainint_funct_ > 0) evaluate_domain_integrals();

    // write state vectors
    output_state();

    // add restart data
    if (step_ % uprestart_ == 0 and step_ != 0) output_restart();
  }

  return;
}  // TimIntImpl::Output

/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                       vuong 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROFLUIDMULTIPHASE::TimIntImpl::dof_row_map(unsigned nds) const
{
  return Teuchos::rcp(discret_->dof_row_map(nds), false);
}

/*----------------------------------------------------------------------*
 | output of solution vector to BINIO                       vuong 08/16 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Map> POROFLUIDMULTIPHASE::TimIntImpl::ArteryDofRowMap() const
{
  return strategy_->ArteryDofRowMap();
}

/*-----------------------------------------------------------------------*
 | access to block system matrix of artery poro problem kremheller 05/18 |
 *-----------------------------------------------------------------------*/
Teuchos::RCP<Core::LinAlg::BlockSparseMatrixBase>
POROFLUIDMULTIPHASE::TimIntImpl::artery_porofluid_sysmat() const
{
  return strategy_->artery_porofluid_sysmat();
}

/*----------------------------------------------------------------------*
 | return artery residual for coupled system           kremheller 05/18 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> POROFLUIDMULTIPHASE::TimIntImpl::ArteryPorofluidRHS() const
{
  return strategy_->ArteryPorofluidRHS();
}


/*==========================================================================*
 |                                                                          |
 | protected:                                                               |
 |                                                                          |
 *==========================================================================*/

/*==========================================================================*/
// general framework
/*==========================================================================*/


/*----------------------------------------------------------------------*
 | evaluate Dirichlet boundary conditions at t_{n+1}        vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::apply_dirichlet_bc(
    const double time, Teuchos::RCP<Epetra_Vector> prenp, Teuchos::RCP<Epetra_Vector> predt)
{
  // time measurement: apply Dirichlet conditions
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:      + apply dirich cond.");

  // needed parameters
  Teuchos::ParameterList p;
  p.set("total time", time);  // actual time t_{n+1}
  p.set<const Core::UTILS::FunctionManager*>(
      "function_manager", &Global::Problem::Instance()->FunctionManager());

  // predicted Dirichlet values
  // \c  prenp then also holds prescribed new Dirichlet values
  discret_->ClearState();
  discret_->evaluate_dirichlet(p, prenp, predt, Teuchos::null, Teuchos::null, dbcmaps_);
  discret_->ClearState();

  return;
}  // POROFLUIDMULTIPHASE::TimIntImpl::apply_dirichlet_bc


/*----------------------------------------------------------------------*
 | contains the residual scaling and addition of Neumann terms          |
 |                                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::scaling_and_neumann()
{
  // scaling to get true residual vector for all time integration schemes
  trueresidual_->Update(residual_scaling(), *residual_, 0.0);

  // add Neumann b.c. scaled with a factor due to time discretization
  add_neumann_to_residual();

  return;
}  // TimIntImpl::scaling_and_neumann


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary conditions                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::apply_neumann_bc(
    const Teuchos::RCP<Epetra_Vector>& neumann_loads  //!< Neumann loads
)
{
  // prepare load vector
  neumann_loads->PutScalar(0.0);

  // create parameter list
  Teuchos::ParameterList condparams;

  // action for elements
  condparams.set<int>("action", POROFLUIDMULTIPHASE::bd_calc_Neumann);

  // set time for evaluation of point Neumann conditions as parameter depending on time integration
  // scheme line/surface/volume Neumann conditions use the time stored in the time parameter class
  set_time_for_neumann_evaluation(condparams);

  // evaluate Neumann boundary conditions at time t_{n+alpha_F} (generalized alpha) or time t_{n+1}
  // (otherwise)
  discret_->evaluate_neumann(condparams, *neumann_loads);
  discret_->ClearState();

  return;
}  // POROFLUIDMULTIPHASE::TimIntImpl::apply_neumann_bc


/*----------------------------------------------------------------------*
 | contains the assembly process for matrix and rhs         vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::assemble_mat_and_rhs()
{
  // time measurement: element calls
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:       + element calls");

  // get cpu time
  const double tcpuele = Teuchos::Time::wallTime();

  // zero out matrix entries
  sysmat_->Zero();

  // reset the residual vector
  residual_->PutScalar(0.0);

  // create parameter list for elements
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", POROFLUIDMULTIPHASE::calc_mat_and_rhs);

  // clean up, just in case ...
  discret_->ClearState();

  // set vector values needed by elements
  // add state vectors according to time-integration scheme
  add_time_integration_specific_vectors();

  // call loop over elements
  discret_->evaluate(eleparams, sysmat_, Teuchos::null, residual_, Teuchos::null, Teuchos::null);

  // clean up
  discret_->ClearState();

  // potential residual scaling and potential addition of Neumann terms
  scaling_and_neumann();

  // finalize assembly of system matrix
  sysmat_->Complete();

  // end time measurement for element
  double mydtele = Teuchos::Time::wallTime() - tcpuele;
  discret_->Comm().MaxAll(&mydtele, &dtele_, 1);

  return;
}  // TimIntImpl::assemble_mat_and_rhs

/*------------------------------------------------------------------------*
 | contains the assembly process for 'valid...-vectors'  kremheller 09/18 |
 *------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::evaluate_valid_volume_frac_press_and_spec()
{
  // time measurement: element calls
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:       + element calls");

  // reset valid_volfracpress_dofs and _spec_dofs-vector
  valid_volfracpress_dofs_->PutScalar(0.0);
  valid_volfracspec_dofs_->PutScalar(0.0);

  // create parameter list for elements
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", POROFLUIDMULTIPHASE::calc_valid_dofs);

  // clean up, just in case ...
  discret_->ClearState();

  // set vector values needed by elements
  // add state vectors according to time-integration scheme
  add_time_integration_specific_vectors();

  // call loop over elements (with valid volume fraction pressure DOFs)
  discret_->evaluate(eleparams, Teuchos::null, Teuchos::null, Teuchos::null,
      valid_volfracpress_dofs_, valid_volfracspec_dofs_);

  // clean up
  discret_->ClearState();

  return;
}

/*------------------------------------------------------------------------------*
 | apply the additional DBC for the volume fraction pressures  kremheller 09/18 |
 *------------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::apply_additional_dbc_for_vol_frac_press()
{
  const Epetra_Map* elecolmap = discret_->ElementColMap();
  std::vector<int> mydirichdofs(0);

  // we identify the volume fraction pressure dofs which do not have a physical meaning and set
  // a DBC on them
  for (int i = 0; i < elecolmap->NumMyElements(); ++i)
  {
    // dynamic_cast necessary because virtual inheritance needs runtime information
    Discret::ELEMENTS::PoroFluidMultiPhase* myele =
        dynamic_cast<Discret::ELEMENTS::PoroFluidMultiPhase*>(
            discret_->gElement(elecolmap->GID(i)));

    const Core::Mat::Material& material = *(myele->Material());

    // check the material
    if (material.MaterialType() != Core::Materials::m_fluidporo_multiphase and
        material.MaterialType() != Core::Materials::m_fluidporo_multiphase_reactions)
      FOUR_C_THROW("only poro multiphase and poro multiphase reactions material valid");

    // cast
    const Mat::FluidPoroMultiPhase& multiphasemat =
        static_cast<const Mat::FluidPoroMultiPhase&>(material);

    const int numfluidphases = multiphasemat.NumFluidPhases();
    const int numvolfrac = multiphasemat.NumVolFrac();
    const int nummat = multiphasemat.NumMat();

    // this is only necessary if we have volume fractions present
    // TODO: this works only if we have the same number of phases in every element
    if (nummat == numfluidphases) return;

    Core::Nodes::Node** nodes = myele->Nodes();
    for (int inode = 0; inode < (myele->num_node()); inode++)
    {
      if (nodes[inode]->Owner() == myrank_)
      {
        std::vector<int> dofs = discret_->Dof(0, nodes[inode]);

        for (int idof = numfluidphases + numvolfrac; idof < nummat; ++idof)
        {
          // if not already in original dirich map     &&   if it is not a valid volume fraction
          // pressure dof identified with < 1
          if (dbcmaps_->CondMap()->LID(dofs[idof]) == -1 &&
              (int)(*valid_volfracpress_dofs_)[discret_->dof_row_map()->LID(dofs[idof])] < 1)
            if (not(std::find(mydirichdofs.begin(), mydirichdofs.end(), dofs[idof]) !=
                    mydirichdofs.end()))
            {
              mydirichdofs.push_back(dofs[idof]);
              phinp_->ReplaceGlobalValue(dofs[idof], 0, 0.0);
            }
        }
      }
    }
  }

  // build map
  int nummydirichvals = mydirichdofs.size();
  Teuchos::RCP<Epetra_Map> dirichmap =
      Teuchos::rcp(new Epetra_Map(-1, nummydirichvals, mydirichdofs.data(), 0, discret_->Comm()));

  // build vector of maps
  std::vector<Teuchos::RCP<const Epetra_Map>> condmaps;
  condmaps.push_back(dirichmap);
  condmaps.push_back(dbcmaps_->CondMap());

  // combined map
  Teuchos::RCP<Epetra_Map> condmerged = Core::LinAlg::MultiMapExtractor::MergeMaps(condmaps);
  *dbcmaps_with_volfracpress_ = Core::LinAlg::MapExtractor(*(discret_->dof_row_map()), condmerged);

  return;
}

void POROFLUIDMULTIPHASE::TimIntImpl::apply_starting_dbc()
{
  const auto& elecolmap = *discret_->ElementColMap();
  std::vector<int> dirichlet_dofs(0);
  const int num_poro_dofs = discret_->NumDof(0, discret_->lRowNode(0));

  for (int ele_idx = 0; ele_idx < elecolmap.NumMyElements(); ++ele_idx)
  {
    const auto& current_element = *discret_->gElement(elecolmap.GID(ele_idx));
    const auto* const nodes = current_element.Nodes();

    for (int node_idx = 0; node_idx < (current_element.num_node()); node_idx++)
    {
      const auto* const current_node = nodes[node_idx];
      if (current_node->Owner() == myrank_)
      {
        const std::vector<int> gid_node_dofs = discret_->Dof(0, current_node);

        for (int dof_idx = 0; dof_idx < num_poro_dofs; ++dof_idx)
        {
          if (starting_dbc_onoff_[dof_idx])
          {
            auto const gid = gid_node_dofs[dof_idx];
            if (std::find(dirichlet_dofs.begin(), dirichlet_dofs.end(), gid) ==
                dirichlet_dofs.end())
            {
              // LID returns -1 if not found in this map/on this processor
              if (dbcmaps_with_volfracpress_->CondMap()->LID(gid) == -1)
              {
                dirichlet_dofs.push_back(gid);
              }
              const double dbc_value = Global::Problem::Instance()
                                           ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(
                                               starting_dbc_funct_[dof_idx] - 1)
                                           .evaluate(current_node->X().data(), time_, 0);
              phinp_->ReplaceGlobalValue(gid, 0, dbc_value);
            }
          }
        }
      }
    }
  }

  // build combined DBC map
  Teuchos::RCP<Epetra_Map> additional_map = Teuchos::rcp(
      new Epetra_Map(-1, dirichlet_dofs.size(), dirichlet_dofs.data(), 0, discret_->Comm()));

  std::vector<Teuchos::RCP<const Epetra_Map>> condition_maps;
  condition_maps.emplace_back(additional_map);
  condition_maps.push_back(dbcmaps_with_volfracpress_->CondMap());

  Teuchos::RCP<Epetra_Map> combined_map =
      Core::LinAlg::MultiMapExtractor::MergeMaps(condition_maps);
  *dbcmaps_starting_condition_ =
      Core::LinAlg::MapExtractor(*(discret_->dof_row_map()), combined_map);
}

/*----------------------------------------------------------------------*
 | assembly process for fluid-structure-coupling matrix kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::assemble_fluid_struct_coupling_mat(
    Teuchos::RCP<Core::LinAlg::SparseOperator> k_fs)
{
  // time measurement: element calls
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:       + element calls");

  // get cpu time
  const double tcpuele = Teuchos::Time::wallTime();

  // create parameter list for elements
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", POROFLUIDMULTIPHASE::calc_fluid_struct_coupl_mat);

  // clean up, just in case ...
  discret_->ClearState();

  // set vector values needed by elements
  // add state vectors according to time-integration scheme
  add_time_integration_specific_vectors();

  // create strategy for assembly of solid pressure
  Core::FE::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
      1,                                       // structdofset for column
      k_fs,                                    // fluid-mechanical matrix
      Teuchos::null,                           // no other matrix or vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // Evaluate coupling matrix
  discret_->evaluate(eleparams, fluidstrategy);

  // clean up
  discret_->ClearState();

  // end time measurement for element
  dtele_ = Teuchos::Time::wallTime() - tcpuele;

  return;
}

/*----------------------------------------------------------------------*
 | assembly process for fluid-scatra-coupling matrix   kremheller 06/17 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::assemble_fluid_scatra_coupling_mat(
    Teuchos::RCP<Core::LinAlg::SparseOperator> k_pfs)
{
  // time measurement: element calls
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:       + element calls");

  // get cpu time
  const double tcpuele = Teuchos::Time::wallTime();

  // create parameter list for elements
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", POROFLUIDMULTIPHASE::calc_fluid_scatra_coupl_mat);

  // clean up, just in case ...
  discret_->ClearState();

  // set vector values needed by elements
  // add state vectors according to time-integration scheme
  add_time_integration_specific_vectors();

  // create strategy for assembly of solid pressure
  Core::FE::AssembleStrategy fluidstrategy(0,  // fluiddofset for row
      3,                                       // scatradofset for column
      k_pfs,                                   // fluid-scatra matrix
      Teuchos::null,                           // no other matrix or vectors
      Teuchos::null, Teuchos::null, Teuchos::null);

  // Evaluate coupling matrix
  discret_->evaluate(eleparams, fluidstrategy);

  // clean up
  discret_->ClearState();

  // end time measurement for element
  dtele_ = Teuchos::Time::wallTime() - tcpuele;

  return;
}

/*----------------------------------------------------------------------*
 | contains the nonlinear iteration loop                    vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::nonlinear_solve()
{
  // time measurement: nonlinear iteration
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:   + nonlin. iteration/lin. solve");

  // out to screen
  print_header();
  print_time_step_info();

  // print header of convergence table to screen
  print_convergence_header();

  //------------------------------ turn adaptive solver tolerance on/off
  const bool isadapttol = (Core::UTILS::IntegralValue<int>(poroparams_, "ADAPTCONV"));
  const double adaptolbetter = poroparams_.get<double>("ADAPTCONV_BETTER");
  const double abstolres = poroparams_.get<double>("ABSTOLRES");
  double actresidual(0.0);

  // prepare Newton-Raphson iteration
  iternum_ = 0;

  // start Newton-Raphson iteration
  while (true)
  {
    iternum_++;

    // call elements to calculate system matrix and rhs and assemble
    // note: DBC is applied herein
    evaluate();

    // abort nonlinear iteration if desired
    if (abort_nonlin_iter(iternum_, itemax_, abstolres, actresidual)) break;

    // initialize increment vector
    increment_->PutScalar(0.0);

    linear_solve(isadapttol, actresidual, adaptolbetter);

    //------------------------------------------------ update solution vector
    UpdateIter(strategy_->CombinedIncrement(increment_));

  }  // nonlinear iteration

  return;
}  // TimIntImpl::nonlinear_solve

/*--------------------------------------------------------------------------*
 | solve linear system of equations                        kremheller 04/18 |
 *--------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::linear_solve(
    bool isadapttol, double actresidual, double adaptolbetter)
{
  // get cpu time
  const double tcpusolve = Teuchos::Time::wallTime();

  // time measurement: call linear solver
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:       + call linear solver");

  // do adaptive linear solver tolerance (not in first solve)
  Core::LinAlg::SolverParams solver_params;
  if (isadapttol && iternum_ > 1)
  {
    solver_params.nonlin_tolerance = ittolres_;
    solver_params.nonlin_residual = actresidual;
    solver_params.lin_tol_better = adaptolbetter;
  }

  strategy_->linear_solve(solver_, sysmat_, increment_, residual_, solver_params);
  // solver_->Solve(sysmat_->EpetraOperator(),increment_,residual_,true,1,Teuchos::null);

  solver_->ResetTolerance();

  // end time measurement for solver
  double mydtsolve = Teuchos::Time::wallTime() - tcpusolve;
  discret_->Comm().MaxAll(&mydtsolve, &dtsolve_, 1);

  return;
}

/*----------------------------------------------------------------------*
 | check if to stop the nonlinear iteration                 vuong 08/16 |
 *----------------------------------------------------------------------*/
bool POROFLUIDMULTIPHASE::TimIntImpl::abort_nonlin_iter(
    const int itnum, const int itemax, const double abstolres, double& actresidual)
{
  //----------------------------------------------------- compute norms
  std::vector<double> preresnorm;
  std::vector<double> incprenorm;
  std::vector<double> prenorm;
  strategy_->CalculateNorms(preresnorm, incprenorm, prenorm, increment_);

  std::vector<double> relinc(prenorm.size());

  for (std::size_t i = 0; i < prenorm.size(); ++i)
  {
    // care for the case that nothing really happens in the pressure
    if (prenorm[i] < 1.0e-6) prenorm[i] = 1.0;
    relinc[i] = incprenorm[i] / prenorm[i];
  }

  const double maxres = *std::max_element(preresnorm.begin(), preresnorm.end());
  const double maxrelinc = *std::max_element(relinc.begin(), relinc.end());

  //-------------------------------------------------- output to screen
  // special case of very first iteration step: solution increment is not yet available
  if (itnum == 1)
  {
    // print first line of convergence table to screen
    print_convergence_values_first_iter(itnum, itemax, ittolinc_, preresnorm);
    // we have to solve at least once --> return false
    return false;
  }

  // ordinary case later iteration steps: solution increment can be printed and convergence check
  // should be done
  else
  {
    // print current line of convergence table to screen
    print_convergence_values(itnum, itemax, ittolinc_, preresnorm, incprenorm, prenorm);

    // convergence check
    if (maxres <= ittolres_ and maxrelinc <= ittolinc_)
    {
      // print finish line of convergence table to screen
      print_convergence_finish_line();

      return true;
    }
  }

  // abort iteration, when there's nothing more to do! -> more robustness
  // absolute tolerance for deciding if residual is (already) zero
  // prevents additional solver calls that will not improve the residual anymore
  if ((maxres < abstolres))
  {
    // print finish line of convergence table to screen
    print_convergence_finish_line();

    return true;
  }


  if ((itnum == itemax))
  {
    switch (divcontype_)
    {
      case Inpar::POROFLUIDMULTIPHASE::divcont_continue:
      {
        // warn if itemax is reached without convergence, but proceed to
        // next timestep...
        if (myrank_ == 0)
        {
          std::cout << "+---------------------------------------------------------------+"
                    << std::endl;
          std::cout << "|            >>>>>> not converged in itemax steps!              |"
                    << std::endl;
          std::cout << "+---------------------------------------------------------------+"
                    << std::endl
                    << std::endl;
        }
        break;
      }
      case Inpar::POROFLUIDMULTIPHASE::divcont_stop:
      {
        FOUR_C_THROW("Porofluid multiphase solver not converged in itemax steps!");
        break;
      }
      default:
        FOUR_C_THROW("unknown divercont action!");
        break;
    }
    // yes, we stop the iteration
    return true;
  }

  // return the maximum residual value -> used for adaptivity of linear solver tolerance
  actresidual = std::max(maxres, maxrelinc);

  // check for INF's and NaN's before going on...
  for (std::size_t i = 0; i < prenorm.size(); ++i)
    if (std::isnan(incprenorm[i]) or std::isnan(prenorm[i]) or std::isnan(preresnorm[i]))
      FOUR_C_THROW("calculated vector norm is NaN.");
  for (std::size_t i = 0; i < prenorm.size(); ++i)
    if (std::isinf(incprenorm[i]) or std::isinf(prenorm[i]) or std::isinf(preresnorm[i]))
      FOUR_C_THROW("calculated vector norm is INF.");

  return false;
}  // TimIntImpl::abort_nonlin_iter

/*----------------------------------------------------------------------*
 | Print Header to screen                              kremheller 05/17 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::print_header()
{
  if (myrank_ == 0)
  {
    std::cout << "\n";
    std::cout << "+--------------------------------------------------------------------------------"
                 "--------------------------------+\n";
    std::cout << "| PORO MULTIPHASE FLUID SOLVER                                                   "
                 "                                |\n";
  }
  return;
}

/*----------------------------------------------------------------------------*
 | reconstruct pressures and saturation from current solution     vuong 08/16 |
 *--------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::reconstruct_pressures_and_saturations()
{
  if (output_satpress_)
  {
    // reset
    pressure_->PutScalar(0.0);
    saturation_->PutScalar(0.0);

    // create parameter list for elements
    Teuchos::ParameterList eleparams;

    // action for elements
    eleparams.set<int>("action", POROFLUIDMULTIPHASE::calc_pres_and_sat);

    // set vector values needed by elements
    discret_->ClearState();

    // add state vectors according to time-integration scheme
    add_time_integration_specific_vectors();

    // initialize counter vector (will store how many times the node has been evaluated)
    Teuchos::RCP<Epetra_Vector> counter =
        Core::LinAlg::CreateVector(*discret_->dof_row_map(), true);

    // call loop over elements
    discret_->evaluate(eleparams, Teuchos::null, Teuchos::null, pressure_, saturation_, counter);

    discret_->ClearState();

    // dummy way: the values have been assembled too many times -> just divide by number of
    // evaluations
    for (int i = 0; i < discret_->dof_row_map()->NumMyElements(); i++)
    {
      (*pressure_)[i] *= 1.0 / (*counter)[i];
      (*saturation_)[i] *= 1.0 / (*counter)[i];
    }
  }

  // reconstruct also the solid pressures
  if (output_solidpress_) reconstruct_solid_pressures();

  return;
}

/*----------------------------------------------------------------------------*
 | reconstruct pressures from current solution                    vuong 08/16 |
 *---------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::reconstruct_solid_pressures()
{
  // reset
  solidpressure_->PutScalar(0.0);

  // create parameter list for elements
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", POROFLUIDMULTIPHASE::calc_solidpressure);

  // set vector values needed by elements
  discret_->ClearState();

  // add state vectors according to time-integration scheme
  add_time_integration_specific_vectors();

  // initialize counter vector (will store how many times the node has been evaluated)
  Teuchos::RCP<Epetra_Vector> counter =
      Core::LinAlg::CreateVector(*discret_->dof_row_map(nds_solidpressure_), true);

  // create strategy for assembly of solid pressure
  Core::FE::AssembleStrategy strategysolidpressure(
      nds_solidpressure_, 0, Teuchos::null, Teuchos::null, solidpressure_, counter, Teuchos::null);

  // call loop over elements
  discret_->evaluate(eleparams, strategysolidpressure);

  discret_->ClearState();

  // dummy way: the values have been assembled too many times -> just divide by number of
  // evaluations
  for (int i = 0; i < discret_->dof_row_map(nds_solidpressure_)->NumMyElements(); i++)
  {
    (*solidpressure_)[i] *= 1.0 / (*counter)[i];
  }
}

/*----------------------------------------------------------------------------*
 | reconstruct velocicities from current solution                vuong 08/16 |
 *--------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::ReconstructFlux()
{
  if (fluxrecon_ == Inpar::POROFLUIDMULTIPHASE::gradreco_none) return;

  // create the parameters for the discretization
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", POROFLUIDMULTIPHASE::recon_flux_at_nodes);

  const int dim = Global::Problem::Instance()->NDim();
  // we assume same number of dofs per node in the whole dis here
  const int totalnumdof = discret_->NumDof(0, discret_->lRowNode(0));
  const int numvec = totalnumdof * dim;

  // add state vectors according to time-integration scheme
  add_time_integration_specific_vectors();

  switch (fluxrecon_)
  {
    case Inpar::POROFLUIDMULTIPHASE::gradreco_l2:
    {
      const auto& solverparams = Global::Problem::Instance()->SolverParams(fluxreconsolvernum_);
      flux_ = Core::FE::compute_nodal_l2_projection(discret_, "phinp_fluid", numvec, eleparams,
          solverparams, Global::Problem::Instance()->solver_params_callback());
      break;
    }
    default:
      FOUR_C_THROW("unknown method for recovery of fluxes!");
      break;
  }
}

/*----------------------------------------------------------------------------*
 *---------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::calculate_phase_velocities()
{
  phase_velocities_->PutScalar(0.0);

  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", POROFLUIDMULTIPHASE::calc_phase_velocities);

  discret_->ClearState();

  add_time_integration_specific_vectors();

  discret_->EvaluateScalars(eleparams, phase_velocities_);
}

/*----------------------------------------------------------------------------*
 | reconstruct porosity from current solution                kremheller 04/17 |
 *---------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::reconstruct_porosity()
{
  // time measurement: reconstruction of porosity
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:    + reconstruct porosity");

  // reset
  porosity_->PutScalar(0.0);

  // create parameter list for elements
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", POROFLUIDMULTIPHASE::calc_porosity);

  // set vector values needed by elements
  discret_->ClearState();

  // add state vectors according to time-integration scheme
  add_time_integration_specific_vectors();

  // initialize counter vector (will store how many times the node has been evaluated)
  Teuchos::RCP<Epetra_Vector> counter =
      Core::LinAlg::CreateVector(*discret_->dof_row_map(nds_solidpressure_), true);

  // create strategy for assembly of porosity
  Core::FE::AssembleStrategy strategyporosity(
      nds_solidpressure_, 0, Teuchos::null, Teuchos::null, porosity_, counter, Teuchos::null);

  // call loop over elements
  discret_->evaluate(eleparams, strategyporosity);

  discret_->ClearState();

  // dummy way: the values have been assembled too many times -> just divide by number of
  // evaluations
  for (int i = 0; i < discret_->dof_row_map(nds_solidpressure_)->NumMyElements(); i++)
  {
    (*porosity_)[i] *= 1.0 / (*counter)[i];
  }

  return;
}

/*----------------------------------------------------------------------------*
 | evaluate domain integrals                                 kremheller 03/19 |
 *---------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::evaluate_domain_integrals()
{
  // time measurement: evaluation of domain integrals
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:    + evaluate domain integrals");

  // create parameter list for elements
  Teuchos::ParameterList eleparams;

  // action for elements
  eleparams.set<int>("action", POROFLUIDMULTIPHASE::calc_domain_integrals);

  // set vector values needed by elements
  discret_->ClearState();

  // add state vectors according to time-integration scheme
  add_time_integration_specific_vectors();

  // evaluate
  discret_->EvaluateScalars(eleparams, domain_integrals_);

  if (myrank_ == 0)  // only one core should do output
  {
    // set filename and file
    const std::string filename(
        Global::Problem::Instance()->OutputControlFile()->file_name() + ".domain_int" + ".csv");
    std::ofstream file;

    if (step_ == 0)
    {
      // inform user that function output has been requested
      std::cout << "\nINFO for domain integrals:\nUser requested " << num_domainint_funct_
                << " function(s) which will be integrated over the entire domain" << std::endl;
      std::cout << "The result will be written into " << filename << "\n" << std::endl;

      // open file and write header
      file.open(filename, std::fstream::trunc);
      file << "Step,Time";
      for (int i = 0; i < num_domainint_funct_; i++)
      {
        file << ",Function_" << domainint_funct_[i];
      }
      file << "\n";
      file.close();
    }

    // usual output
    file.open(filename, std::fstream::app);
    // step, time and results for each function
    file << step_ << "," << time_;
    for (int i = 0; i < num_domainint_funct_; i++)
      file << "," << std::setprecision(14) << (*domain_integrals_)[i];

    // close file
    file << "\n";
    file.close();
  }  // if myrank == 0

  discret_->ClearState();

  return;
}

/*----------------------------------------------------------------------*
 | print header of convergence table to screen              vuong 08/16 |
 *----------------------------------------------------------------------*/
inline void POROFLUIDMULTIPHASE::TimIntImpl::print_convergence_header()
{
  if (myrank_ == 0)
  {
    if (artery_coupling_active_)
    {
      std::cout
          << "+------------+-------------------+--------------+--------------+-------------------+-"
             "-------------+--------------+\n"
          << "|- step/max -|- tol-res   [norm]-|-- pre-res ---|--- 1D-res ---|- "
             "tol-relinc[norm]-|-- pre-inc ---|--- 1D-inc ---|"
          << std::endl;
    }
    else
    {
      std::cout
          << "+------------+-------------------+--------------+-------------------+--------------+-"
             "--"
             "--------------------------+\n"
          << "|- step/max -|- tol-res   [norm]-|-- pre-res ---|- tol-relinc[norm]-|-- pre-inc ---|"
          << std::endl;
    }
  }

  return;
}  // POROFLUIDMULTIPHASE::TimIntImpl::print_convergence_header


/*----------------------------------------------------------------------*
 | print first line of convergence table to screen          vuong 08/16 |
 *----------------------------------------------------------------------*/
inline void POROFLUIDMULTIPHASE::TimIntImpl::print_convergence_values_first_iter(
    const int& itnum,                      //!< current Newton-Raphson iteration step
    const int& itemax,                     //!< maximum number of Newton-Raphson iteration steps
    const double& ittol,                   //!< relative tolerance for Newton-Raphson scheme
    const std::vector<double>& preresnorm  //!< L2 norm of pressure residual
)
{
  if (myrank_ == 0)
  {
    std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itemax << "   | "
              << std::setw(10) << std::setprecision(3) << std::scientific << ittolres_ << " ["
              << std::setw(3) << VectorNormString(vectornormfres_).c_str() << "]  | ";

    for (std::size_t i = 0; i < preresnorm.size(); ++i)
      std::cout << std::setw(10) << std::setprecision(3) << std::scientific << preresnorm[i]
                << "   | ";

    std::cout << std::setw(10) << std::setprecision(3) << std::scientific << ittolinc_ << " ["
              << std::setw(3) << VectorNormString(vectornorminc_).c_str() << "]  |";

    for (std::size_t i = 0; i < preresnorm.size(); ++i) std::cout << "      --      |";
    std::cout << " (    --   ,te=" << std::setw(10) << std::setprecision(3) << std::scientific
              << dtele_ << ")" << std::endl;
  }

  return;
}  // POROFLUIDMULTIPHASE::TimIntImpl::print_convergence_values_first_iter


/*----------------------------------------------------------------------*
 | print current line of convergence table to screen        vuong 08/16 |
 *----------------------------------------------------------------------*/
inline void POROFLUIDMULTIPHASE::TimIntImpl::print_convergence_values(
    const int& itnum,                       //!< current Newton-Raphson iteration step
    const int& itemax,                      //!< maximum number of Newton-Raphson iteration steps
    const double& ittol,                    //!< relative tolerance for Newton-Raphson scheme
    const std::vector<double>& preresnorm,  //!< norm of pressure residual
    const std::vector<double>& incprenorm,  //!< norm of pressure increment
    const std::vector<double>& prenorm      //!< norm of pressure state vector
)
{
  if (myrank_ == 0)
  {
    std::cout << "|  " << std::setw(3) << itnum << "/" << std::setw(3) << itemax << "   | "
              << std::setw(10) << std::setprecision(3) << std::scientific << ittolres_ << " ["
              << std::setw(3) << VectorNormString(vectornormfres_).c_str() << "]  | ";
    for (std::size_t i = 0; i < preresnorm.size(); ++i)
      std::cout << std::setw(10) << std::setprecision(3) << std::scientific << preresnorm[i]
                << "   | ";
    std::cout << std::setw(10) << std::setprecision(3) << std::scientific << ittolres_ << " ["
              << std::setw(3) << VectorNormString(vectornorminc_).c_str() << "]  | ";
    for (std::size_t i = 0; i < preresnorm.size(); ++i)
      std::cout << std::setw(10) << std::setprecision(3) << std::scientific
                << incprenorm[i] / prenorm[i] << "   | ";
    std::cout << std::setw(10) << std::setprecision(3) << std::scientific << dtsolve_
              << ",te=" << std::setw(10) << std::setprecision(3) << std::scientific << dtele_ << ")"
              << std::endl;
  }

  return;
}  // POROFLUIDMULTIPHASE::TimIntImpl::print_convergence_values


/*----------------------------------------------------------------------*
 | print finish line of convergence table to screen         vuong 08/16 |
 *----------------------------------------------------------------------*/
inline void POROFLUIDMULTIPHASE::TimIntImpl::print_convergence_finish_line()
{
  if (myrank_ == 0)
  {
    if (artery_coupling_active_)
    {
      std::cout << "+------------+-------------------+--------------+--------------+---------------"
                   "----+--------------+--------------+"
                << std::endl;
    }
    else
    {
      std::cout
          << "+------------+-------------------+--------------+-------------------+--------------+"
          << std::endl;
    }
  }

  return;
}  // POROFLUIDMULTIPHASE::TimIntImpl::print_convergence_finish_line


/*----------------------------------------------------------------------*
 |  write current state to BINIO                            vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::output_state()
{
  // solution
  output_->write_vector("phinp_fluid", phinp_);
  // time derivative of solution
  // output_->write_vector("phidtnp", phidtnp_);
  if (output_satpress_)
  {
    // pressures
    output_->write_vector("pressure", pressure_);
    // saturations
    output_->write_vector("saturation", saturation_);
  }

  // solid pressure
  if (output_solidpress_)
  {
    // convert dof-based Epetra vector into node-based Epetra multi-vector for postprocessing
    Teuchos::RCP<Epetra_MultiVector> solidpressure_multi =
        POROFLUIDMULTIPHASE::UTILS::ConvertDofVectorToNodeBasedMultiVector(
            *discret_, *solidpressure_, nds_solidpressure_, 1);

    output_->write_vector("solidpressure", solidpressure_multi, Core::IO::nodevector);
  }

  // displacement field
  if (isale_)
  {
    Teuchos::RCP<const Epetra_Vector> dispnp = discret_->GetState(nds_disp_, "dispnp");
    if (dispnp == Teuchos::null)
      FOUR_C_THROW("Cannot extract displacement field from discretization");

    // convert dof-based Epetra vector into node-based Epetra multi-vector for postprocessing
    Teuchos::RCP<Epetra_MultiVector> dispnp_multi =
        POROFLUIDMULTIPHASE::UTILS::ConvertDofVectorToNodeBasedMultiVector(
            *discret_, *dispnp, nds_disp_, nsd_);

    output_->write_vector("dispnp", dispnp_multi, Core::IO::nodevector);
  }
  // fluxes
  if (flux_ != Teuchos::null)
  {
    // post_ensight does not support multivectors based on the dofmap
    // for now, I create single vectors that can be handled by the filter

    const int dim = Global::Problem::Instance()->NDim();
    const int numdof = discret_->NumDof(0, discret_->lRowNode(0));
    // get the noderowmap
    const Epetra_Map* noderowmap = discret_->NodeRowMap();
    for (int k = 0; k < numdof; k++)
    {
      Teuchos::RCP<Epetra_MultiVector> flux_k =
          Teuchos::rcp(new Epetra_MultiVector(*noderowmap, 3, true));

      std::ostringstream temp;
      temp << k + 1;
      std::string name = "flux_" + temp.str();
      for (int i = 0; i < flux_k->MyLength(); ++i)
      {
        // get value for each component of flux vector
        for (int idim = 0; idim < dim; idim++)
        {
          double value = ((*flux_)[k * dim + idim])[i];
          int err = flux_k->ReplaceMyValue(i, idim, value);
          if (err != 0) FOUR_C_THROW("Detected error in ReplaceMyValue");
        }
      }
      output_->write_vector(name, flux_k, Core::IO::nodevector);
    }
  }

  if (output_phase_velocities_)
  {
    const int num_dim = Global::Problem::Instance()->NDim();
    const int num_poro_dof = discret_->NumDof(0, discret_->lRowNode(0));

    const Epetra_Map* element_row_map = discret_->ElementRowMap();

    for (int k = 0; k < num_poro_dof; k++)
    {
      Teuchos::RCP<Epetra_MultiVector> velocity_k =
          Teuchos::rcp(new Epetra_MultiVector(*element_row_map, num_dim, true));

      for (int i = 0; i < velocity_k->MyLength(); ++i)
      {
        for (int idim = 0; idim < num_dim; idim++)
        {
          double value = ((*phase_velocities_)[k * num_dim + idim])[i];
          int err = velocity_k->ReplaceMyValue(i, idim, value);
          if (err != 0) FOUR_C_THROW("Detected error in ReplaceMyValue");
        }
      }

      std::string output_name = "velocity_" + std::to_string(k + 1);
      output_->write_vector(output_name, velocity_k, Core::IO::elementvector);
    }
  }

  // porosity
  if (output_porosity_)
  {
    // convert dof-based Epetra vector into node-based Epetra multi-vector for postprocessing
    Teuchos::RCP<Epetra_MultiVector> porosity_multi =
        POROFLUIDMULTIPHASE::UTILS::ConvertDofVectorToNodeBasedMultiVector(
            *discret_, *porosity_, nds_solidpressure_, 1);

    output_->write_vector("porosity", porosity_multi, Core::IO::nodevector);
  }

  return;
}  // TimIntImpl::output_state


/*----------------------------------------------------------------------*
 | increment time and step for next iteration               vuong 08/16 |
 *----------------------------------------------------------------------*/
inline void POROFLUIDMULTIPHASE::TimIntImpl::increment_time_and_step()
{
  step_ += 1;
  time_ += dt_;
}


/*----------------------------------------------------------------------*
 |  calculate error compared to analytical solution         vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::evaluate_error_compared_to_analytical_sol()
{
  if (calcerr_ == Inpar::POROFLUIDMULTIPHASE::calcerror_no) return;
  FOUR_C_THROW("Error calculation not yet implemented! Element evaluation is missing.");

  // create the parameters for the error calculation
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", POROFLUIDMULTIPHASE::calc_error);
  eleparams.set("total time", time_);
  eleparams.set<int>("calcerrorflag", calcerr_);

  switch (calcerr_)
  {
    case Inpar::POROFLUIDMULTIPHASE::calcerror_byfunction:
    {
      const int errorfunctnumber = poroparams_.get<int>("CALCERRORNO");
      if (errorfunctnumber < 1)
        FOUR_C_THROW("invalid value of paramter CALCERRORNO for error function evaluation!");

      eleparams.set<int>("error function number", errorfunctnumber);
      break;
    }
    default:
      FOUR_C_THROW("Cannot calculate error. Unknown type of analytical test problem");
      break;
  }

  // set vector values needed by elements
  discret_->ClearState();
  discret_->set_state("phinp_fluid", phinp_);

  // get (squared) error values
  Teuchos::RCP<Core::LinAlg::SerialDenseVector> errors =
      Teuchos::rcp(new Core::LinAlg::SerialDenseVector(4));
  discret_->EvaluateScalars(eleparams, errors);
  discret_->ClearState();

  // std::vector containing
  // [0]: relative L2 pressure error
  // [1]: relative H1 pressure error
  Teuchos::RCP<std::vector<double>> relerror = Teuchos::rcp(new std::vector<double>(2));

  if (std::abs((*errors)[2]) > 1e-14)
    (*relerror)[0] = sqrt((*errors)[0]) / sqrt((*errors)[2]);
  else
    (*relerror)[0] = sqrt((*errors)[0]);
  if (std::abs((*errors)[2]) > 1e-14)
    (*relerror)[1] = sqrt((*errors)[1]) / sqrt((*errors)[3]);
  else
    (*relerror)[1] = sqrt((*errors)[1]);

  if (myrank_ == 0)
  {
    // print last error in a separate file

    const std::string simulation = Global::Problem::Instance()->OutputControlFile()->file_name();
    const std::string fname = simulation + "_pressure_time.relerror";

    if (step_ == 0)
    {
      std::ofstream f;
      f.open(fname.c_str());
      f << "#| Step | Time | rel. L2-error  | rel. H1-error  |\n";
      f << std::setprecision(10) << step_ << " " << std::setw(1) << std::setprecision(5) << time_
        << std::setw(1) << std::setprecision(6) << " " << (*relerror)[0] << std::setw(1)
        << std::setprecision(6) << " " << (*relerror)[1] << "\n";

      f.flush();
      f.close();
    }
    else
    {
      std::ofstream f;
      f.open(fname.c_str(), std::fstream::ate | std::fstream::app);
      f << std::setprecision(10) << step_ << " " << std::setw(3) << std::setprecision(5) << time_
        << std::setw(1) << std::setprecision(6) << " " << (*relerror)[0] << std::setw(1)
        << std::setprecision(6) << " " << (*relerror)[1] << "\n";

      f.flush();
      f.close();
    }
  }

  return;
}  // POROFLUIDMULTIPHASE::TimIntImpl::evaluate_error_compared_to_analytical_sol


/*----------------------------------------------------------------------*
 | build linear system tangent matrix, rhs/force residual   vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::evaluate()
{
  // call elements to calculate system matrix and rhs and assemble
  assemble_mat_and_rhs();

  // perform finite difference check on time integrator level
  if (fdcheck_ == Inpar::POROFLUIDMULTIPHASE::fdcheck_global) fd_check();

  // Apply Dirichlet Boundary Condition
  prepare_system_for_newton_solve();

  // evaluate mesh tying
  strategy_->evaluate();
}

/*----------------------------------------------------------------------*
 | basically apply Dirichlet BC                        kremheller 06/17 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::prepare_system_for_newton_solve()
{
  // Apply Dirichlet boundary conditions to system of equations
  // residual values are supposed to be zero at Dirichlet boundaries
  {
    // time measurement: application of DBC to system
    TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:       + apply DBC to system");

    if (time_ <= starting_dbc_time_end_)
    {
      Core::LinAlg::apply_dirichlet_to_system(
          *sysmat_, *increment_, *residual_, *zeros_, *(dbcmaps_starting_condition_->CondMap()));
    }
    else
    {
      Core::LinAlg::apply_dirichlet_to_system(
          *sysmat_, *increment_, *residual_, *zeros_, *(dbcmaps_with_volfracpress_->CondMap()));
    }
  }
}

/*----------------------------------------------------------------------*
 | iterative update of scalars                              vuong 08/16  |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc)
{
  Teuchos::RCP<const Epetra_Vector> extractedinc = strategy_->extract_and_update_iter(inc);
  // store incremental vector to be available for convergence check
  // if incremental vector is received from outside for coupled problem
  increment_->Update(1.0, *extractedinc, 0.0);

  // update scalar values by adding increments
  phinp_->Update(1.0, *extractedinc, 1.0);

  // compute time derivative at time n+1
  compute_time_derivative();

}  // UpdateIter


/*----------------------------------------------------------------------*
 | set convective velocity field                            vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::set_velocity_field(
    Teuchos::RCP<const Epetra_Vector> vel  //!< velocity vector
)
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE: set convective velocity field");

  //---------------------------------------------------------------------------
  // preliminaries
  //---------------------------------------------------------------------------
  if (nds_vel_ == -1)
    FOUR_C_THROW("Dof set for velocity degreess of freedom has not been assigned!");

  if (vel == Teuchos::null) FOUR_C_THROW("Velocity state is Teuchos::null");

  if (nds_vel_ >= discret_->NumDofSets())
    FOUR_C_THROW("Too few dofsets on poro fluid discretization!");

  if (not vel->Map().SameAs(*discret_->dof_row_map(nds_vel_)))
    FOUR_C_THROW(
        "Map of given velocity and associated dof row map in poro fluid discretization"
        " do not match!");

  // provide discretization with velocity
  set_state(nds_vel_, "velocity field", vel);

  return;

}  // TimIntImpl::set_velocity_field

/*----------------------------------------------------------------------*
 | set state on discretization                                          |
 |                                                          vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::set_state(
    unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state)
{
  // provide discretization with velocity
  discret_->set_state(nds, name, state);
}


/*----------------------------------------------------------------------*
 |  set initial field for phi                                 vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::SetInitialField(
    const Inpar::POROFLUIDMULTIPHASE::InitialField init, const int startfuncno)
{
  switch (init)
  {
    case Inpar::POROFLUIDMULTIPHASE::initfield_zero_field:
    {
      phin_->PutScalar(0.0);
      phinp_->PutScalar(0.0);
      break;
    }
    case Inpar::POROFLUIDMULTIPHASE::initfield_field_by_function:
    {
      const Epetra_Map* dofrowmap = discret_->dof_row_map();

      // loop all nodes on the processor
      for (int lnodeid = 0; lnodeid < discret_->NumMyRowNodes(); lnodeid++)
      {
        // get the processor local node
        Core::Nodes::Node* lnode = discret_->lRowNode(lnodeid);
        // the set of degrees of freedom associated with the node
        std::vector<int> nodedofset = discret_->Dof(0, lnode);

        int numdofs = nodedofset.size();
        for (int k = 0; k < numdofs; ++k)
        {
          const int dofgid = nodedofset[k];
          int doflid = dofrowmap->LID(dofgid);
          // evaluate component k of spatial function
          double initialval = Global::Problem::Instance()
                                  ->FunctionById<Core::UTILS::FunctionOfSpaceTime>(startfuncno - 1)
                                  .evaluate(lnode->X().data(), time_, k);
          int err = phin_->ReplaceMyValues(1, &initialval, &doflid);
          if (err != 0) FOUR_C_THROW("dof not on proc");
        }
      }

      // initialize also the solution vector. These values are a pretty good guess for the
      // solution after the first time step (much better than starting with a zero vector)
      phinp_->Update(1.0, *phin_, 0.0);

      break;
    }
    case Inpar::POROFLUIDMULTIPHASE::initfield_field_by_condition:
    {
      // set initial field for ALL existing scatra fields in condition
      const std::string field = "PoroMultiFluid";

      const int numdof = discret_->NumDof(0, discret_->lRowNode(0));

      // get initial field conditions
      std::vector<Core::Conditions::Condition*> initfieldconditions(0);
      discret_->GetCondition("Initfield", initfieldconditions);

      if (not initfieldconditions.size())
        FOUR_C_THROW(
            "Tried to evaluate initial field by condition without a corresponding condition "
            "defined on the PoroMultiFluid discretization!");
      std::set<int> numdofpernode;
      for (unsigned icond = 0; icond < initfieldconditions.size(); icond++)
      {
        const int condmaxnumdofpernode = numdof;

        if (condmaxnumdofpernode != 0) numdofpernode.insert(condmaxnumdofpernode);
      }

      if (numdofpernode.empty()) FOUR_C_THROW("No DOFs defined on initial field condtion!");

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
    default:
      FOUR_C_THROW("Unknown option for initial field: %d", init);
      break;
  }  // switch(init)

  return;
}  // TimIntImpl::SetInitialField

/*----------------------------------------------------------------------*
 | create result test for this field                        vuong 08/16  |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Core::UTILS::ResultTest> POROFLUIDMULTIPHASE::TimIntImpl::CreateFieldTest()
{
  strategy_->CreateFieldTest();
  return Teuchos::rcp(new POROFLUIDMULTIPHASE::ResultTest(*this));
}

/*----------------------------------------------------------------------*
 |  read restart data                                  kremheller 04/18 |
 -----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::read_restart(const int step)
{
  strategy_->read_restart(step);
  return;
}

/*--------------------------------------------------------------------*
 | calculate init time derivatives of state variables kremheller 03/17 |
 *--------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::calc_initial_time_derivative()
{
  // time measurement
  TEUCHOS_FUNC_TIME_MONITOR("POROFLUIDMULTIPHASE:       + calculate initial time derivative");

  // initial screen output
  if (myrank_ == 0)
    std::cout << "POROFLUIDMULTIPHASE: calculating initial time derivative of state variables on "
                 "discretization \""
              << discret_->Name().c_str() << "\" (step " << Step() << ", time " << Time()
              << ") ... ... ";

  // reset global system matrix
  sysmat_->Zero();

  // reset the residual vector
  residual_->PutScalar(0.0);

  // evaluate Dirichlet and Neumann boundary conditions at time t = 0 to ensure consistent
  // computation of initial time derivative vector Dirichlet boundary conditions should be
  // consistent with initial field
  apply_dirichlet_bc(time_, phinp_, Teuchos::null);
  compute_intermediate_values();
  apply_neumann_bc(neumann_loads_);

  // create and fill parameter list for elements
  Teuchos::ParameterList eleparams;
  eleparams.set<int>("action", POROFLUIDMULTIPHASE::calc_initial_time_deriv);

  // add state vectors according to time integration scheme
  discret_->ClearState();
  add_time_integration_specific_vectors();

  // We evaluate the discretization such that
  // mp * dphidt + msp * dphidt + msat * dphidt = - rhsfac (sdivvel + diff + reac)
  // (      only    matrix                    )   (          only rhs              )
  // later we will also have to scale the system matrix with rhsfac
  discret_->evaluate(eleparams, sysmat_, residual_);
  discret_->ClearState();

  // potential residual scaling and potential addition of Neumann terms
  scaling_and_neumann();

  // We have to Scale the system matrix consistently
  // TODO: this is kind of a hack, does it work for other schemes than one-step theta??
  // sysmat_->Scale(1.0/residual_scaling());
  residual_->Scale(residual_scaling());

  // finalize assembly of system matrix
  sysmat_->Complete();

  // solve global system of equations for initial time derivative of state variables
  Core::LinAlg::SolverParams solver_params;
  solver_params.refactor = true;
  solver_params.reset = true;
  solver_->Solve(sysmat_->EpetraOperator(), phidtnp_, residual_, solver_params);

  // copy solution
  phidtn_->Update(1.0, *phidtnp_, 0.0);

  // reset global system matrix and its graph, since we solved a very special problem with a special
  // sparsity pattern
  sysmat_->Reset();

  // reset solver
  solver_->Reset();

  // reset true residual vector computed during assembly of the standard global system of equations,
  // since not yet needed
  trueresidual_->PutScalar(0.0);

  double maxval = 0.0;
  phidtnp_->MaxValue(&maxval);
  // final screen output
  if (myrank_ == 0)
  {
    std::cout << "done!" << std::endl;
    std::cout << "MAX value: " << maxval << std::endl;
  }

  return;
}
/*----------------------------------------------------------------------------------------------*
 | finite difference check for scalar transport system matrix (for debugging only)   vuong 08/16 |
 *----------------------------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::TimIntImpl::fd_check()
{
  // initial screen output
  if (myrank_ == 0)
    std::cout << std::endl << "FINITE DIFFERENCE CHECK FOR SYSTEM MATRIX" << std::endl;

  // make a copy of state variables to undo perturbations later
  Teuchos::RCP<Epetra_Vector> phinp_original = Teuchos::rcp(new Epetra_Vector(*phinp_));

  // make a copy of system matrix as Epetra_CrsMatrix
  Teuchos::RCP<Epetra_CrsMatrix> sysmat_original = Teuchos::null;
  if (Teuchos::rcp_dynamic_cast<Core::LinAlg::SparseMatrix>(sysmat_) != Teuchos::null)
    sysmat_original = (new Core::LinAlg::SparseMatrix(
                           *(Teuchos::rcp_static_cast<Core::LinAlg::SparseMatrix>(sysmat_))))
                          ->EpetraMatrix();
  else if (Teuchos::rcp_dynamic_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat_) != Teuchos::null)
    sysmat_original =
        (new Core::LinAlg::SparseMatrix(
             *(Teuchos::rcp_static_cast<Core::LinAlg::BlockSparseMatrixBase>(sysmat_)->Merge())))
            ->EpetraMatrix();
  else
    FOUR_C_THROW("Type of system matrix unknown!");
  sysmat_original->FillComplete();

  // make a copy of system right-hand side vector
  Teuchos::RCP<Epetra_Vector> rhs_original = Teuchos::rcp(new Epetra_Vector(*residual_));

  // initialize counter for system matrix entries with failing finite difference check
  int counter(0);

  // initialize tracking variable for maximum absolute and relative errors
  double maxabserr(0.);
  double maxrelerr(0.);

  for (int colgid = 0; colgid <= sysmat_original->ColMap().MaxAllGID(); ++colgid)
  {
    // check whether current column index is a valid global column index and continue loop if not
    int collid(sysmat_original->ColMap().LID(colgid));
    int maxcollid(-1);
    discret_->Comm().MaxAll(&collid, &maxcollid, 1);
    if (maxcollid < 0) continue;

    // fill state vector with original state variables
    phinp_->Update(1., *phinp_original, 0.);

    // impose perturbation
    if (phinp_->Map().MyGID(colgid))
      if (phinp_->SumIntoGlobalValue(colgid, 0, fdcheckeps_))
        FOUR_C_THROW(
            "Perturbation could not be imposed on state vector for finite difference check!");

    // carry perturbation over to state vectors at intermediate time stages if necessary
    compute_intermediate_values();
    compute_time_derivative();

    // calculate element right-hand side vector for perturbed state
    assemble_mat_and_rhs();

    // Now we compare the difference between the current entries in the system matrix
    // and their finite difference approximations according to
    // entries ?= (residual_perturbed - residual_original) / epsilon

    // Note that the residual_ vector actually denotes the right-hand side of the linear
    // system of equations, i.e., the negative system residual.
    // To account for errors due to numerical cancellation, we additionally consider
    // entries + residual_original / epsilon ?= residual_perturbed / epsilon

    // Note that we still need to evaluate the first comparison as well. For small entries in the
    // system matrix, the second comparison might yield good agreement in spite of the entries being
    // wrong!
    for (int rowlid = 0; rowlid < discret_->dof_row_map()->NumMyElements(); ++rowlid)
    {
      // get global index of current matrix row
      const int rowgid = sysmat_original->RowMap().GID(rowlid);
      if (rowgid < 0) FOUR_C_THROW("Invalid global ID of matrix row!");

      // get relevant entry in current row of original system matrix
      double entry(0.);
      int length = sysmat_original->NumMyEntries(rowlid);
      int numentries;
      std::vector<double> values(length);
      std::vector<int> indices(length);
      sysmat_original->ExtractMyRowCopy(rowlid, length, numentries, values.data(), indices.data());
      for (int ientry = 0; ientry < length; ++ientry)
      {
        if (sysmat_original->ColMap().GID(indices[ientry]) == colgid)
        {
          entry = values[ientry];
          break;
        }
      }

      // finite difference suggestion (first divide by epsilon and then add for better conditioning)
      const double fdval =
          -(*residual_)[rowlid] / fdcheckeps_ + (*rhs_original)[rowlid] / fdcheckeps_;

      // confirm accuracy of first comparison
      if (abs(fdval) > 1.e-17 and abs(fdval) < 1.e-15)
        FOUR_C_THROW("Finite difference check involves values too close to numerical zero!");

      // absolute and relative errors in first comparison
      const double abserr1 = entry - fdval;
      if (abs(abserr1) > maxabserr) maxabserr = abs(abserr1);
      double relerr1(0.);
      if (abs(entry) > 1.e-17)
        relerr1 = abserr1 / abs(entry);
      else if (abs(fdval) > 1.e-17)
        relerr1 = abserr1 / abs(fdval);
      if (abs(relerr1) > maxrelerr) maxrelerr = abs(relerr1);

      // evaluate first comparison
      if (abs(relerr1) > fdchecktol_)
      {
        std::cout << "sysmat[" << rowgid << "," << colgid << "]:  " << entry << "   ";
        std::cout << "finite difference suggestion:  " << fdval << "   ";
        std::cout << "absolute error:  " << abserr1 << "   ";
        std::cout << "relative error:  " << relerr1 << std::endl;

        counter++;
      }

      // first comparison OK
      else
      {
        // left-hand side in second comparison
        const double left = entry - (*rhs_original)[rowlid] / fdcheckeps_;

        // right-hand side in second comparison
        const double right = -(*residual_)[rowlid] / fdcheckeps_;

        // confirm accuracy of second comparison
        if (abs(right) > 1.e-17 and abs(right) < 1.e-15)
          FOUR_C_THROW("Finite difference check involves values too close to numerical zero!");

        // absolute and relative errors in second comparison
        const double abserr2 = left - right;
        if (abs(abserr2) > maxabserr) maxabserr = abs(abserr2);
        double relerr2(0.);
        if (abs(left) > 1.e-17)
          relerr2 = abserr2 / abs(left);
        else if (abs(right) > 1.e-17)
          relerr2 = abserr2 / abs(right);
        if (abs(relerr2) > maxrelerr) maxrelerr = abs(relerr2);

        // evaluate second comparison
        if (abs(relerr2) > fdchecktol_)
        {
          std::cout << "sysmat[" << rowgid << "," << colgid << "]-rhs[" << rowgid
                    << "]/eps:  " << left << "   ";
          std::cout << "-rhs_perturbed[" << rowgid << "]/eps:  " << right << "   ";
          std::cout << "absolute error:  " << abserr2 << "   ";
          std::cout << "relative error:  " << relerr2 << std::endl;

          counter++;
        }
      }
    }
  }

  // communicate tracking variables
  int counterglobal(0);
  discret_->Comm().SumAll(&counter, &counterglobal, 1);
  double maxabserrglobal(0.);
  discret_->Comm().MaxAll(&maxabserr, &maxabserrglobal, 1);
  double maxrelerrglobal(0.);
  discret_->Comm().MaxAll(&maxrelerr, &maxrelerrglobal, 1);

  // final screen output
  if (myrank_ == 0)
  {
    if (counterglobal)
    {
      printf(
          "--> FAILED AS LISTED ABOVE WITH %d CRITICAL MATRIX ENTRIES IN TOTAL\n\n", counterglobal);
      FOUR_C_THROW("Finite difference check failed for scalar transport system matrix!");
    }
    else
      printf(
          "--> PASSED WITH MAXIMUM ABSOLUTE ERROR %+12.5e AND MAXIMUM RELATIVE ERROR %+12.5e\n\n",
          maxabserrglobal, maxrelerrglobal);
  }

  // undo perturbations of state variables
  phinp_->Update(1., *phinp_original, 0.);
  compute_intermediate_values();
  compute_time_derivative();

  // recompute system matrix and right-hand side vector based on original state variables
  assemble_mat_and_rhs();

  return;
}

/*----------------------------------------------------------------*
 | return arterial network time integrator       kremheller 04/18 |
 *----------------------------------------------------------------*/
Teuchos::RCP<Adapter::ArtNet> POROFLUIDMULTIPHASE::TimIntImpl::ArtNetTimInt()
{
  return strategy_->ArtNetTimInt();
}

FOUR_C_NAMESPACE_CLOSE
