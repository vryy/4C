/*----------------------------------------------------------------------*/
/*! \file
\brief time-integration scheme with extensions for
       cardiac monodomain problems

\level 2


*/
/*----------------------------------------------------------------------*/

#include "baci_scatra_timint_cardiac_monodomain_scheme.hpp"

#include "baci_global_data.hpp"
#include "baci_io.hpp"
#include "baci_lib_discret.hpp"
#include "baci_scatra_ele_action.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ljag 01/14 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntCardiacMonodomainOST::TimIntCardiacMonodomainOST(
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<CORE::LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      TimIntCardiacMonodomain(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntOneStepTheta(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                              ehrl 01/14 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainOST::Setup()
{
  // call Setup()-functions of base classes
  // note: this order is important
  TimIntOneStepTheta::Setup();

  TimIntCardiacMonodomain::Setup();

  return;
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainOST::Update()
{
  // Standard Update
  TimIntOneStepTheta::Update();

  // time update of myocard material
  TimIntCardiacMonodomain::ElementMaterialTimeUpdate();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainOST::WriteRestart() const
{
  // Call function from baseclass
  TimIntOneStepTheta::WriteRestart();

  // Cardiac Monodomain specific
  output_->WriteMesh(
      step_, time_);  // add info to control file for reading all variables in restart

  return;
}

/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainOST::ReadRestart(
    const int step, Teuchos::RCP<IO::InputControl> input)
{
  // Call function from baseclass
  TimIntOneStepTheta::ReadRestart(step, input);

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
    reader = Teuchos::rcp(new IO::DiscretizationReader(
        discret_, GLOBAL::Problem::Instance()->InputControlFile(), step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  // Cardiac Monodomain specific
  reader->ReadVector(activation_time_np_, "activation_time_np");
  reader->ReadHistoryData(step);  // Read all saved data in nodes and elements und call nodal and
                                  // element Unpacking each global variable has to be read

  return;
}

/*--------------------------------------------------------------------------*
 | add global state vectors specific for time-integration scheme  hoe 06/16 |
 *--------------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainOST::AddTimeIntegrationSpecificVectors(
    bool forcedincrementalsolver)
{
  // Call function from baseclass
  TimIntOneStepTheta::AddTimeIntegrationSpecificVectors(forcedincrementalsolver);
  discret_->SetState("phin", phin_);

  return;
}

/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ljag 01/14 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntCardiacMonodomainBDF2::TimIntCardiacMonodomainBDF2(
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<CORE::LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      TimIntCardiacMonodomain(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntBDF2(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                              ehrl 01/14 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainBDF2::Setup()
{
  // call Setup()-functions of base classes
  // note: this order is important
  TimIntBDF2::Setup();

  TimIntCardiacMonodomain::Setup();

  return;
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainBDF2::Update()
{
  // Standard Update
  TimIntBDF2::Update();

  // time update of myocard material
  TimIntCardiacMonodomain::ElementMaterialTimeUpdate();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainBDF2::WriteRestart() const
{
  // Call function from baseclass
  TimIntBDF2::WriteRestart();

  // Cardiac Monodomain specific
  output_->WriteMesh(
      step_, time_);  // add info to control file for reading all variables in restart

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainBDF2::ReadRestart(
    const int step, Teuchos::RCP<IO::InputControl> input)
{
  // Call function from baseclass
  TimIntBDF2::ReadRestart(step, input);

  Teuchos::RCP<IO::DiscretizationReader> reader(Teuchos::null);
  if (input == Teuchos::null)
    reader = Teuchos::rcp(new IO::DiscretizationReader(
        discret_, GLOBAL::Problem::Instance()->InputControlFile(), step));
  else
    reader = Teuchos::rcp(new IO::DiscretizationReader(discret_, input, step));

  // Cardiac Monodomain specific
  reader->ReadVector(activation_time_np_, "activation_time_np");
  reader->ReadHistoryData(step);  // Read all saved data in nodes and elements und call nodal and
                                  // element Unpacking each global variable has to be read

  return;
}


/*----------------------------------------------------------------------*
 |  Constructor (public)                                     ljag 01/14 |
 *----------------------------------------------------------------------*/
SCATRA::TimIntCardiacMonodomainGenAlpha::TimIntCardiacMonodomainGenAlpha(
    Teuchos::RCP<DRT::Discretization> actdis, Teuchos::RCP<CORE::LINALG::Solver> solver,
    Teuchos::RCP<Teuchos::ParameterList> params,
    Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
    Teuchos::RCP<Teuchos::ParameterList> extraparams, Teuchos::RCP<IO::DiscretizationWriter> output)
    : ScaTraTimIntImpl(actdis, solver, sctratimintparams, extraparams, output),
      TimIntCardiacMonodomain(actdis, solver, params, sctratimintparams, extraparams, output),
      TimIntGenAlpha(actdis, solver, sctratimintparams, extraparams, output)
{
  return;
}


/*----------------------------------------------------------------------*
 |  initialize time integration                              ehrl 01/14 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainGenAlpha::Setup()
{
  // call Init()-functions of base classes
  // note: this order is important
  TimIntGenAlpha::Setup();

  TimIntCardiacMonodomain::Setup();

  return;
}

/*----------------------------------------------------------------------*
 | current solution becomes most recent solution of next timestep       |
 |                                                            gjb 08/08 |
 *----------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainGenAlpha::Update()
{
  // Standard Update
  TimIntGenAlpha::Update();

  // time update of myocard material
  TimIntCardiacMonodomain::ElementMaterialTimeUpdate();

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainGenAlpha::WriteRestart() const
{
  // Call function from baseclass
  TimIntGenAlpha::WriteRestart();

  // Cardiac Monodomain specific
  output_->WriteMesh(
      step_, time_);  // add info to control file for reading all variables in restart

  return;
}


/*----------------------------------------------------------------------*
 |                                                            gjb 08/08 |
 -----------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainGenAlpha::ReadRestart(
    const int step, Teuchos::RCP<IO::InputControl> input)
{
  // Call function from baseclass
  TimIntGenAlpha::ReadRestart(step, input);

  IO::DiscretizationReader reader(discret_, GLOBAL::Problem::Instance()->InputControlFile(), step);

  // Cardiac Monodomain specific
  reader.ReadVector(activation_time_np_, "activation_time_np");
  reader.ReadHistoryData(step);  // Read all saved data in nodes and elements und call nodal and
                                 // element Unpacking each global variable has to be read

  return;
}

/*--------------------------------------------------------------------------*
 | add global state vectors specific for time-integration scheme  hoe 12/16 |
 *--------------------------------------------------------------------------*/
void SCATRA::TimIntCardiacMonodomainGenAlpha::AddTimeIntegrationSpecificVectors(
    bool forcedincrementalsolver)
{
  // Call function from baseclass
  TimIntGenAlpha::AddTimeIntegrationSpecificVectors(forcedincrementalsolver);

  if (incremental_ or forcedincrementalsolver) discret_->SetState("phin", phin_);

  return;
}

FOUR_C_NAMESPACE_CLOSE
