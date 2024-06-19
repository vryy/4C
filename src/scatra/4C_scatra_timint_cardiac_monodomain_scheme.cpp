/*----------------------------------------------------------------------*/
/*! \file
\brief time-integration scheme with extensions for
       cardiac monodomain problems

\level 2


*/
/*----------------------------------------------------------------------*/

#include "4C_scatra_timint_cardiac_monodomain_scheme.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_scatra_ele_action.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ljag 01/14 |
 *----------------------------------------------------------------------*/
ScaTra::TimIntCardiacMonodomainOST::TimIntCardiacMonodomainOST(
    Teuchos::RCP<Core::FE::Discretization> actdis, Teuchos::RCP<Core::LinAlg::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      TimIntCardiacMonodomain(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntOneStepTheta(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                              ehrl 01/14 |
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainOST::setup()
{
  // call setup()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::setup();

  TimIntCardiacMonodomain::setup();

  return;
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainOST::update()
{
  // Standard Update
  TimIntOneStepTheta::update();

  // time update of myocard material
  TimIntCardiacMonodomain::element_material_time_update();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainOST::write_restart() const
{
  // Call function from baseclass
  TimIntOneStepTheta::write_restart();

  // Cardiac Monodomain specific
  output_->write_mesh(
      step_, time_);  // add info to control file for reading all variables in restart

  return;
}

/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainOST::read_restart(
    const int step, Teuchos::RCP<Core::IO::InputControl> input)
{
  // Call function from baseclass
  TimIntOneStepTheta::read_restart(step, input);

  Teuchos::RCP<Core::IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
    reader = Teuchos::rcp(new Core::IO::DiscretizationReader(
        discret_, Global::Problem::Instance()->InputControlFile(), step));
  else
    reader = Teuchos::rcp(new Core::IO::DiscretizationReader(discret_, input, step));

  // Cardiac Monodomain specific
  reader->read_vector(activation_time_np_, "activation_time_np");
  reader->read_history_data(step);  // Read all saved data in nodes and elements und call nodal and
                                    // element Unpacking each global variable has to be read

  return;
}

/*--------------------------------------------------------------------------*
 | add global state vectors specific for time-integration scheme  hoe 06/16 |
 *--------------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainOST::add_time_integration_specific_vectors(
    bool forcedincrementalsolver)
{
  // Call function from baseclass
  TimIntOneStepTheta::add_time_integration_specific_vectors(forcedincrementalsolver);
  discret_->set_state("phin", phin_);

  return;
}

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ljag 01/14 |
 *----------------------------------------------------------------------*/
ScaTra::TimIntCardiacMonodomainBDF2::TimIntCardiacMonodomainBDF2(
    Teuchos::RCP<Core::FE::Discretization> actdis, Teuchos::RCP<Core::LinAlg::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      TimIntCardiacMonodomain(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntBDF2(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                              ehrl 01/14 |
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainBDF2::setup()
{
  // call setup()-functions of base classes
  // note: this order is important
  TimIntBDF2::setup();

  TimIntCardiacMonodomain::setup();

  return;
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainBDF2::update()
{
  // Standard Update
  TimIntBDF2::update();

  // time update of myocard material
  TimIntCardiacMonodomain::element_material_time_update();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainBDF2::write_restart() const
{
  // Call function from baseclass
  TimIntBDF2::write_restart();

  // Cardiac Monodomain specific
  output_->write_mesh(
      step_, time_);  // add info to control file for reading all variables in restart

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainBDF2::read_restart(
    const int step, Teuchos::RCP<Core::IO::InputControl> input)
{
  // Call function from baseclass
  TimIntBDF2::read_restart(step, input);

  Teuchos::RCP<Core::IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
    reader = Teuchos::rcp(new Core::IO::DiscretizationReader(
        discret_, Global::Problem::Instance()->InputControlFile(), step));
  else
    reader = Teuchos::rcp(new Core::IO::DiscretizationReader(discret_, input, step));

  // Cardiac Monodomain specific
  reader->read_vector(activation_time_np_, "activation_time_np");
  reader->read_history_data(step);  // Read all saved data in nodes and elements und call nodal and
                                    // element Unpacking each global variable has to be read

  return;
}


/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ljag 01/14 |
 *----------------------------------------------------------------------*/
ScaTra::TimIntCardiacMonodomainGenAlpha::TimIntCardiacMonodomainGenAlpha(
    Teuchos::RCP<Core::FE::Discretization> actdis, Teuchos::RCP<Core::LinAlg::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams,
    Teuchos::RCP<Core::IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      TimIntCardiacMonodomain(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntGenAlpha(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                              ehrl 01/14 |
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainGenAlpha::setup()
{
  // call init()-functions of base classes
  // note: this order is important
  TimIntGenAlpha::setup();

  TimIntCardiacMonodomain::setup();

  return;
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainGenAlpha::update()
{
  // Standard Update
  TimIntGenAlpha::update();

  // time update of myocard material
  TimIntCardiacMonodomain::element_material_time_update();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainGenAlpha::write_restart() const
{
  // Call function from baseclass
  TimIntGenAlpha::write_restart();

  // Cardiac Monodomain specific
  output_->write_mesh(
      step_, time_);  // add info to control file for reading all variables in restart

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainGenAlpha::read_restart(
    const int step, Teuchos::RCP<Core::IO::InputControl> input)
{
  // Call function from baseclass
  TimIntGenAlpha::read_restart(step, input);

  Core::IO::DiscretizationReader reader(
      discret_, Global::Problem::Instance()->InputControlFile(), step);

  // Cardiac Monodomain specific
  reader.read_vector(activation_time_np_, "activation_time_np");
  reader.read_history_data(step);  // Read all saved data in nodes and elements und call nodal and
                                   // element Unpacking each global variable has to be read

  return;
}

/*--------------------------------------------------------------------------*
 | add global state vectors specific for time-integration scheme  hoe 12/16 |
 *--------------------------------------------------------------------------*/
void ScaTra::TimIntCardiacMonodomainGenAlpha::add_time_integration_specific_vectors(
    bool forcedincrementalsolver)
{
  // Call function from baseclass
  TimIntGenAlpha::add_time_integration_specific_vectors(forcedincrementalsolver);

  if (incremental_ or forcedincrementalsolver) discret_->set_state("phin", phin_);

  return;
}

FOUR_C_NAMESPACE_CLOSE
