/*---------------------------------------------------------------------------*/
/*!
\file particle_timint.cpp

\brief time integration for particle simulations

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_timint.H"

#include "particle_dirichlet_bc.H"
#include "particle_temperature_bc.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_container_bundle.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_io/io.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::TimInt::TimInt(const Teuchos::ParameterList& params)
    : params_(params), time_(0.0), dt_(params.get<double>("TIMESTEP")), modifiedstates_(false)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | destructor                                                 sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::TimInt::~TimInt()
{
  // note: destructor declaration here since at compile-time a complete type
  // of class T as used in class member std::unique_ptr<T> ptr_T_ is required
}

/*---------------------------------------------------------------------------*
 | init particle time integration                             sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimInt::Init()
{
  // init dirichlet boundary condition handler
  InitDirichletBoundaryCondition();

  // init temperature boundary condition handler
  InitTemperatureBoundaryCondition();
}

/*---------------------------------------------------------------------------*
 | setup particle time integration                            sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimInt::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const bool modifiedstates)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // modified velocity and acceleration states
  modifiedstates_ = modifiedstates;

  // setup dirichlet boundary condition handler
  if (dirichletboundarycondition_) dirichletboundarycondition_->Setup(particleengineinterface_);

  // setup temperature boundary condition handler
  if (temperatureboundarycondition_) temperatureboundarycondition_->Setup(particleengineinterface_);

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // set of particle types being excluded from time integration
  std::set<PARTICLEENGINE::TypeEnum> typesexludedfromtimeintegration;

  // boundary particles are not integrated in time
  typesexludedfromtimeintegration.insert(PARTICLEENGINE::BoundaryPhase);

  if (dirichletboundarycondition_)
  {
    // get reference to set of particle types subjected to dirichlet boundary conditions
    const std::set<PARTICLEENGINE::TypeEnum>& typessubjectedtodirichletbc =
        dirichletboundarycondition_->GetParticleTypesSubjectedToDirichletBCSet();

    // particles subjected to dirichlet boundary conditions are not integrated in time
    for (PARTICLEENGINE::TypeEnum currtype : typessubjectedtodirichletbc)
      typesexludedfromtimeintegration.insert(currtype);
  }

  // determine set of particle types to integrate in time
  for (auto typeIt : particlecontainerbundle->GetRefToAllContainersMap())
  {
    // get enum of particle types
    PARTICLEENGINE::TypeEnum particleType = typeIt.first;

    // insert current particle type into set of particles to integrate in time
    if (typesexludedfromtimeintegration.find(particleType) == typesexludedfromtimeintegration.end())
      typestointegrate_.insert(particleType);
  }
}

/*---------------------------------------------------------------------------*
 | write restart of particle time integration                 sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimInt::WriteRestart(const int step, const double time) const
{
  // write restart of dirichlet boundary condition handler
  if (dirichletboundarycondition_) dirichletboundarycondition_->WriteRestart(step, time);

  // write restart of temperature boundary condition handler
  if (temperatureboundarycondition_) temperatureboundarycondition_->WriteRestart(step, time);
}

/*---------------------------------------------------------------------------*
 | read restart of particle time integration                  sfuchs 04/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimInt::ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // read restart of dirichlet boundary condition handler
  if (dirichletboundarycondition_) dirichletboundarycondition_->ReadRestart(reader);

  // read restart of temperature boundary condition handler
  if (temperatureboundarycondition_) temperatureboundarycondition_->ReadRestart(reader);
}

/*---------------------------------------------------------------------------*
 | insert integration dependent states of all particle types  sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimInt::InsertParticleStatesOfParticleTypes(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // insert dbc dependent states of all particle types
  if (dirichletboundarycondition_)
    dirichletboundarycondition_->InsertParticleStatesOfParticleTypes(particlestatestotypes);

  // insert tempbc dependent states of all particle types
  if (temperatureboundarycondition_)
    temperatureboundarycondition_->InsertParticleStatesOfParticleTypes(particlestatestotypes);
}

/*---------------------------------------------------------------------------*
 | time integration scheme specific initialization routine    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimInt::SetInitialStates()
{
  // set particle reference position
  if (dirichletboundarycondition_) dirichletboundarycondition_->SetParticleReferencePosition();

  // set particle reference position
  if (temperatureboundarycondition_) temperatureboundarycondition_->SetParticleReferencePosition();

  // evaluate dirichlet boundary condition
  if (dirichletboundarycondition_)
    dirichletboundarycondition_->EvaluateDirichletBoundaryCondition(0.0, true, true, true);

  // evaluate temperature boundary condition
  if (temperatureboundarycondition_)
    temperatureboundarycondition_->EvaluateTemperatureBoundaryCondition(0.0);
}

/*---------------------------------------------------------------------------*
 | set current time                                           sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimInt::SetCurrentTime(const double currenttime) { time_ = currenttime; }

/*---------------------------------------------------------------------------*
 | init dirichlet boundary condition handler                  sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimInt::InitDirichletBoundaryCondition()
{
  // create dirichlet boundary condition handler
  dirichletboundarycondition_ =
      std::unique_ptr<PARTICLEALGORITHM::DirichletBoundaryConditionHandler>(
          new PARTICLEALGORITHM::DirichletBoundaryConditionHandler(params_));

  // init dirichlet boundary condition handler
  dirichletboundarycondition_->Init();

  // get reference to set of particle types subjected to dirichlet boundary conditions
  const std::set<PARTICLEENGINE::TypeEnum>& typessubjectedtodirichletbc =
      dirichletboundarycondition_->GetParticleTypesSubjectedToDirichletBCSet();

  // no particle types are subjected to dirichlet boundary conditions
  if (typessubjectedtodirichletbc.size() == 0) dirichletboundarycondition_.release();
}

/*---------------------------------------------------------------------------*
 | init temperature boundary condition handler                 meier 09/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimInt::InitTemperatureBoundaryCondition()
{
  // create temperature boundary condition handler
  temperatureboundarycondition_ =
      std::unique_ptr<PARTICLEALGORITHM::TemperatureBoundaryConditionHandler>(
          new PARTICLEALGORITHM::TemperatureBoundaryConditionHandler(params_));

  // init temperature boundary condition handler
  temperatureboundarycondition_->Init();

  // get reference to set of particle types subjected to temperature boundary conditions
  const std::set<PARTICLEENGINE::TypeEnum>& typessubjectedtotempbc =
      temperatureboundarycondition_->GetParticleTypesSubjectedToTemperatureBCSet();

  // no particle types are subjected to temperature boundary conditions
  if (typessubjectedtotempbc.size() == 0) temperatureboundarycondition_.release();
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::TimIntSemiImplicitEuler::TimIntSemiImplicitEuler(
    const Teuchos::ParameterList& params)
    : PARTICLEALGORITHM::TimInt(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | setup particle time integration                            sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimIntSemiImplicitEuler::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const bool modifiedstates)
{
  PARTICLEALGORITHM::TimInt::Setup(particleengineinterface, modifiedstates);

  // safety check
  if (modifiedstates_)
    dserror(
        "modified velocity and acceleration states not implemented yet for semi-implicit Euler "
        "time integration scheme!");
}

/*---------------------------------------------------------------------------*
 | time integration scheme specific initialization routine    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimIntSemiImplicitEuler::SetInitialStates()
{
  // call base class method
  TimInt::SetInitialStates();
}

/*---------------------------------------------------------------------------*
 | time integration scheme specific pre-interaction routine   sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimIntSemiImplicitEuler::PreInteractionRoutine()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::TimIntSemiImplicitEuler::PreInteractionRoutine");

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // iterate over particle types
  for (auto& particleType : typestointegrate_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // update velocity of all particles
    container->UpdateState(1.0, PARTICLEENGINE::Velocity, dt_, PARTICLEENGINE::Acceleration);

    // clear acceleration of all particles
    container->ClearState(PARTICLEENGINE::Acceleration);

    // update position of all particles
    container->UpdateState(1.0, PARTICLEENGINE::Position, dt_, PARTICLEENGINE::Velocity);
  }

  // evaluate dirichlet boundary condition
  if (dirichletboundarycondition_)
    dirichletboundarycondition_->EvaluateDirichletBoundaryCondition(time_, true, true, true);

  // evaluate temperature boundary condition
  if (temperatureboundarycondition_)
    temperatureboundarycondition_->EvaluateTemperatureBoundaryCondition(time_);
}

/*---------------------------------------------------------------------------*
 | time integration scheme specific post-interaction routine  sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimIntSemiImplicitEuler::PostInteractionRoutine()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::TimIntVelocityVerlet::TimIntVelocityVerlet(const Teuchos::ParameterList& params)
    : PARTICLEALGORITHM::TimInt(params), dthalf_(0.5 * dt_)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | time integration scheme specific initialization routine    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimIntVelocityVerlet::SetInitialStates()
{
  // call base class method
  TimInt::SetInitialStates();

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // iterate over particle types
  for (auto& particleType : typestointegrate_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // modified velocity and acceleration states
    if (modifiedstates_)
    {
      // update modified velocity of all particles
      container->UpdateState(0.0, PARTICLEENGINE::ModifiedVelocity, 1.0, PARTICLEENGINE::Velocity);
    }
  }
}

/*---------------------------------------------------------------------------*
 | time integration scheme specific pre-interaction routine   sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimIntVelocityVerlet::PreInteractionRoutine()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::TimIntVelocityVerlet::PreInteractionRoutine");

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // iterate over particle types
  for (auto& particleType : typestointegrate_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // update velocity of all particles
    container->UpdateState(1.0, PARTICLEENGINE::Velocity, dthalf_, PARTICLEENGINE::Acceleration);

    // clear acceleration of all particles
    container->ClearState(PARTICLEENGINE::Acceleration);

    // modified velocity and acceleration states
    if (modifiedstates_)
    {
      // update modified velocity of all particles
      container->UpdateState(0.0, PARTICLEENGINE::ModifiedVelocity, 1.0, PARTICLEENGINE::Velocity);
      container->UpdateState(
          1.0, PARTICLEENGINE::ModifiedVelocity, dthalf_, PARTICLEENGINE::ModifiedAcceleration);

      // clear modified acceleration of all particles
      container->ClearState(PARTICLEENGINE::ModifiedAcceleration);

      // update position of all particles
      container->UpdateState(1.0, PARTICLEENGINE::Position, dt_, PARTICLEENGINE::ModifiedVelocity);
    }
    else
    {
      // update position of all particles
      container->UpdateState(1.0, PARTICLEENGINE::Position, dt_, PARTICLEENGINE::Velocity);
    }
  }

  // evaluate dirichlet boundary condition
  if (dirichletboundarycondition_)
  {
    dirichletboundarycondition_->EvaluateDirichletBoundaryCondition(time_, true, false, true);
    dirichletboundarycondition_->EvaluateDirichletBoundaryCondition(
        time_ - dthalf_, false, true, false);
  }

  // evaluate temperature boundary condition
  if (temperatureboundarycondition_)
    temperatureboundarycondition_->EvaluateTemperatureBoundaryCondition(time_);
}

/*---------------------------------------------------------------------------*
 | time integration scheme specific post-interaction routine  sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::TimIntVelocityVerlet::PostInteractionRoutine()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::TimIntVelocityVerlet::PostInteractionRoutine");

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // iterate over particle types
  for (auto& particleType : typestointegrate_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainerShrdPtr container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // update velocity of all particles
    container->UpdateState(1.0, PARTICLEENGINE::Velocity, dthalf_, PARTICLEENGINE::Acceleration);
  }

  // evaluate dirichlet boundary condition
  if (dirichletboundarycondition_)
    dirichletboundarycondition_->EvaluateDirichletBoundaryCondition(time_, false, true, false);
}
