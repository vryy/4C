/*----------------------------------------------------------------------*/
/*! \file
\brief solution algorithm for stationary problems

\level 1

*/
/*----------------------------------------------------------------------*/

#include "4C_scatra_timint_stat.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_lib_discret.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_timint_meshtying_strategy_base.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor (public)                                      gjb 08/08 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntStationary::TimIntStationary(Teuchos::RCP<DRT::Discretization> actdis,
    Teuchos::RCP<CORE::LINALG::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, params, extraparams, output), fsphinp_(Teuchos::null)
{
  // DO NOT DEFINE ANY STATE VECTORS HERE (i.e., vectors based on row or column maps)
  // this is important since we have problems which require an extended ghosting
  // this has to be done before all state vectors are initialized
}


/*----------------------------------------------------------------------*
 |  initialize time integration                         rasthofer 09/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::Init()
{
  // initialize base class
  ScaTraTimIntImpl::Init();

  // -------------------------------------------------------------------
  // get a vector layout from the discretization to construct matching
  // vectors and matrices
  //                 local <-> global dof numbering
  // -------------------------------------------------------------------
  const Epetra_Map* dofrowmap = discret_->dof_row_map();

  // fine-scale vector
  if (fssgd_ != INPAR::SCATRA::fssugrdiff_no)
    fsphinp_ = CORE::LINALG::CreateVector(*dofrowmap, true);
  if (turbmodel_ != INPAR::FLUID::no_model) FOUR_C_THROW("Turbulence is not stationary problem!");

  // -------------------------------------------------------------------
  // set element parameters
  // -------------------------------------------------------------------
  // note: - this has to be done before element routines are called
  //       - order is important here: for safety checks in set_element_general_parameters(),
  //         we have to know the time-integration parameters
  set_element_time_parameter();
  set_element_general_parameters();
  set_element_turbulence_parameters();

  // setup krylov
  prepare_krylov_projection();

  // safety check
  if (static_cast<bool>(CORE::UTILS::IntegralValue<int>(*params_, "NATURAL_CONVECTION")))
    FOUR_C_THROW("Natural convection for stationary time integration scheme is not implemented!");
}



/*----------------------------------------------------------------------*
 | set time parameter for element evaluation (usual call)    ehrl 11/13 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::set_element_time_parameter(bool forcedincrementalsolver) const
{
  Teuchos::ParameterList eleparams;

  CORE::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
      "action", SCATRA::Action::set_time_parameter, eleparams);
  eleparams.set<bool>("using generalized-alpha time integration", false);
  eleparams.set<bool>("using stationary formulation", true);
  if (!forcedincrementalsolver)
    eleparams.set<bool>("incremental solver", incremental_);
  else
    eleparams.set<bool>("incremental solver", true);

  eleparams.set<double>("time-step length", dta_);
  eleparams.set<double>("total time", time_);
  eleparams.set<double>("time factor", 1.0);
  eleparams.set<double>("alpha_F", 1.0);

  // call standard loop over elements
  discret_->Evaluate(
      eleparams, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::Setup()
{
  // setup base class
  ScaTraTimIntImpl::Setup();

  set_element_nodeset_parameters();
}

/*----------------------------------------------------------------------*
 | set time for evaluation of Neumann boundary conditions      vg 12/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::set_time_for_neumann_evaluation(Teuchos::ParameterList& params)
{
  params.set("total time", time_);
}


/*----------------------------------------------------------------------*
 | set part of the residual vector belonging to the old timestep        |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::set_old_part_of_righthandside()
{
  // call base class routine
  ScaTraTimIntImpl::set_old_part_of_righthandside();

  hist_->PutScalar(0.0);
}


/*----------------------------------------------------------------------*
 | add actual Neumann loads                                    vg 11/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::add_neumann_to_residual()
{
  residual_->Update(1.0, *neumann_loads_, 1.0);
}


/*----------------------------------------------------------------------*
 | AVM3-based scale separation                                 vg 03/09 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::AVM3Separation()
{
  // time measurement: avm3
  TEUCHOS_FUNC_TIME_MONITOR("SCATRA:            + avm3");

  // AVM3 separation
  Sep_->Multiply(false, *phinp_, *fsphinp_);

  // set fine-scale vector
  discret_->set_state("fsphinp", fsphinp_);
}


/*--------------------------------------------------------------------------*
 | add global state vectors specific for time-integration scheme   vg 11/08 |
 *--------------------------------------------------------------------------*/
void SCATRA::TimIntStationary::add_time_integration_specific_vectors(bool forcedincrementalsolver)
{
  // call base class routine
  ScaTraTimIntImpl::add_time_integration_specific_vectors(forcedincrementalsolver);

  discret_->set_state("hist", hist_);
  discret_->set_state("phinp", phinp_);
}


/*----------------------------------------------------------------------*
 |                                                            gjb 09/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::read_restart(const int step, Teuchos::RCP<IO::InputControl> input)
{
  // call base class routine
  ScaTraTimIntImpl::read_restart(step, input);

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
  {
    reader = Teuchos::rcp(new IO::DiscretizationReader(
        discret_, GLOBAL::Problem::Instance()->InputControlFile(), step));
  }
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  time_ = reader->ReadDouble("time");
  step_ = reader->ReadInt("step");

  if (myrank_ == 0)
    std::cout << "Reading ScaTra restart data (time=" << time_ << " ; step=" << step_ << ")"
              << '\n';

  // read state vectors that are needed for restart
  reader->ReadVector(phinp_, "phinp");

  read_restart_problem_specific(step, *reader);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::Update()
{
  // call base class routine
  ScaTraTimIntImpl::Update();

  // compute flux vector field for later output BEFORE time shift of results
  // is performed below !!
  if (calcflux_domain_ != INPAR::SCATRA::flux_none or
      calcflux_boundary_ != INPAR::SCATRA::flux_none)
  {
    if (IsResultStep() or do_boundary_flux_statistics()) CalcFlux(true);
  }
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::TimIntStationary::WriteRestart() const
{
  // call base class routine
  ScaTraTimIntImpl::WriteRestart();

  // This feature enables starting a time-dependent simulation from
  // a non-trivial steady-state solution that was calculated before.
  output_->WriteVector("phin", phinp_);    // for OST and BDF2
  output_->WriteVector("phinm", phinp_);   // for BDF2
  output_->WriteVector("phidtn", zeros_);  // for OST
}

FOUR_C_NAMESPACE_CLOSE
