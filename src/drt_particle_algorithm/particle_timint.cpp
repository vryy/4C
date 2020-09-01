/*---------------------------------------------------------------------------*/
/*! \file
\brief time integration for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_timint.H"

#include "particle_dirichlet_bc.H"
#include "particle_temperature_bc.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_container_bundle.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_io/io.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::TimInt::TimInt(const Teuchos::ParameterList& params)
    : params_(params), time_(0.0), dt_(params.get<double>("TIMESTEP"))
{
  // empty constructor
}

PARTICLEALGORITHM::TimInt::~TimInt() = default;

void PARTICLEALGORITHM::TimInt::Init()
{
  // init dirichlet boundary condition handler
  InitDirichletBoundaryCondition();

  // init temperature boundary condition handler
  InitTemperatureBoundaryCondition();
}

void PARTICLEALGORITHM::TimInt::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // setup dirichlet boundary condition handler
  if (dirichletboundarycondition_) dirichletboundarycondition_->Setup(particleengineinterface_);

  // setup temperature boundary condition handler
  if (temperatureboundarycondition_) temperatureboundarycondition_->Setup(particleengineinterface_);

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // set of particle types being excluded from time integration
  std::set<PARTICLEENGINE::TypeEnum> typesexludedfromtimeintegration;

  // boundary and rigid particles are not integrated in time
  typesexludedfromtimeintegration.insert(
      {PARTICLEENGINE::BoundaryPhase, PARTICLEENGINE::RigidPhase});

  if (dirichletboundarycondition_)
  {
    // get reference to set of particle types subjected to dirichlet boundary conditions
    const std::set<PARTICLEENGINE::TypeEnum>& typessubjectedtodirichletbc =
        dirichletboundarycondition_->GetParticleTypesSubjectedToDirichletBCSet();

    // particles subjected to dirichlet boundary conditions are not integrated in time
    for (PARTICLEENGINE::TypeEnum currtype : typessubjectedtodirichletbc)
      typesexludedfromtimeintegration.insert(currtype);
  }

  // determine set of particle types to be integrated in time
  for (auto& typeEnum : particlecontainerbundle->GetParticleTypes())
    if (not typesexludedfromtimeintegration.count(typeEnum)) typestointegrate_.insert(typeEnum);
}

void PARTICLEALGORITHM::TimInt::WriteRestart(const int step, const double time) const
{
  // write restart of dirichlet boundary condition handler
  if (dirichletboundarycondition_) dirichletboundarycondition_->WriteRestart(step, time);

  // write restart of temperature boundary condition handler
  if (temperatureboundarycondition_) temperatureboundarycondition_->WriteRestart(step, time);
}

void PARTICLEALGORITHM::TimInt::ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // read restart of dirichlet boundary condition handler
  if (dirichletboundarycondition_) dirichletboundarycondition_->ReadRestart(reader);

  // read restart of temperature boundary condition handler
  if (temperatureboundarycondition_) temperatureboundarycondition_->ReadRestart(reader);
}

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

void PARTICLEALGORITHM::TimInt::SetInitialStates()
{
  // add initial random noise to particle position
  AddInitialRandomNoiseToPosition();

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

void PARTICLEALGORITHM::TimInt::SetCurrentTime(const double currenttime) { time_ = currenttime; }

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
  if (typessubjectedtodirichletbc.empty()) dirichletboundarycondition_.release();
}

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
  if (typessubjectedtotempbc.empty()) temperatureboundarycondition_.release();
}

void PARTICLEALGORITHM::TimInt::AddInitialRandomNoiseToPosition()
{
  // init vector of initial position amplitude for each spatial direction
  std::vector<double> amplitude;
  double value;
  std::istringstream amplitudestream(
      Teuchos::getNumericStringParameter(params_, "INITIAL_POSITION_AMPLITUDE"));

  while (amplitudestream >> value) amplitude.push_back(value);

  // safety check
  if (static_cast<int>(amplitude.size()) != 3)
    dserror("dimension (dim = %d) of initial position amplitude vector is wrong!",
        static_cast<int>(amplitude.size()));

  // safety check
  for (double a : amplitude)
    if (a < 0.0)
      dserror("no negative initial position amplitude allowed (set a positive or zero value)!");

  // get magnitude of initial position amplitude
  double temp = 0.0;
  for (double a : amplitude) temp += a * a;
  const double amplitude_norm = std::sqrt(temp);

  // no initial position amplitude defined
  if (not(amplitude_norm > 0.0)) return;

  // safety check
  const double max_amplitude = *std::max_element(amplitude.begin(), amplitude.end());
  if (max_amplitude > particleengineinterface_->MinBinSize())
    dserror("amplitude of noise added to initial position larger than minimum relevant bin size!");

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // iterate over particle types
  for (auto& particleType : typestointegrate_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get pointer to particle state
    double* pos = container->GetPtrToParticleState(PARTICLEENGINE::Position, 0);

    // get particle state dimension
    int statedim = container->GetParticleStateDim(PARTICLEENGINE::Position);

    // iterate over owned particles of current type
    for (int i = 0; i < particlestored; ++i)
    {
      // iterate over spatial dimension
      for (int dim = 0; dim < statedim; ++dim)
      {
        // generate random value
        const double randomvalue = DRT::Problem::Instance()->Random()->Uni();

        // update position of particle
        pos[statedim * i + dim] += randomvalue * amplitude[dim];
      }
    }
  }
}

PARTICLEALGORITHM::TimIntSemiImplicitEuler::TimIntSemiImplicitEuler(
    const Teuchos::ParameterList& params)
    : PARTICLEALGORITHM::TimInt(params)
{
  // empty constructor
}

void PARTICLEALGORITHM::TimIntSemiImplicitEuler::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  PARTICLEALGORITHM::TimInt::Setup(particleengineinterface);

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // iterate over particle types
  for (auto& particleType : typestointegrate_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // safety check
    if (container->HaveStoredState(PARTICLEENGINE::ModifiedVelocity) or
        container->HaveStoredState(PARTICLEENGINE::ModifiedAcceleration))
      dserror(
          "modified velocity and acceleration states not implemented yet for semi-implicit Euler "
          "time integration scheme!");
  }
}

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
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // update velocity of all particles
    container->UpdateState(1.0, PARTICLEENGINE::Velocity, dt_, PARTICLEENGINE::Acceleration);

    // clear acceleration of all particles
    container->ClearState(PARTICLEENGINE::Acceleration);

    // angular velocity and acceleration states
    if (container->HaveStoredState(PARTICLEENGINE::AngularVelocity) and
        container->HaveStoredState(PARTICLEENGINE::AngularAcceleration))
    {
      // update angular velocity of all particles
      container->UpdateState(
          1.0, PARTICLEENGINE::AngularVelocity, dt_, PARTICLEENGINE::AngularAcceleration);

      // clear angular acceleration of all particles
      container->ClearState(PARTICLEENGINE::AngularAcceleration);
    }

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

void PARTICLEALGORITHM::TimIntSemiImplicitEuler::PostInteractionRoutine()
{
  // nothing to do
}

PARTICLEALGORITHM::TimIntVelocityVerlet::TimIntVelocityVerlet(const Teuchos::ParameterList& params)
    : PARTICLEALGORITHM::TimInt(params), dthalf_(0.5 * dt_)
{
  // empty constructor
}

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
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // modified velocity states
    if (container->HaveStoredState(PARTICLEENGINE::ModifiedVelocity))
    {
      // update modified velocity of all particles
      container->UpdateState(0.0, PARTICLEENGINE::ModifiedVelocity, 1.0, PARTICLEENGINE::Velocity);
    }
  }
}

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
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // update velocity of all particles
    container->UpdateState(1.0, PARTICLEENGINE::Velocity, dthalf_, PARTICLEENGINE::Acceleration);

    // clear acceleration of all particles
    container->ClearState(PARTICLEENGINE::Acceleration);

    // angular velocity and acceleration states
    if (container->HaveStoredState(PARTICLEENGINE::AngularVelocity) and
        container->HaveStoredState(PARTICLEENGINE::AngularAcceleration))
    {
      // update angular velocity of all particles
      container->UpdateState(
          1.0, PARTICLEENGINE::AngularVelocity, dthalf_, PARTICLEENGINE::AngularAcceleration);

      // clear angular acceleration of all particles
      container->ClearState(PARTICLEENGINE::AngularAcceleration);
    }

    // modified velocity and acceleration states
    if (container->HaveStoredState(PARTICLEENGINE::ModifiedVelocity) and
        container->HaveStoredState(PARTICLEENGINE::ModifiedAcceleration))
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
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

    // update velocity of all particles
    container->UpdateState(1.0, PARTICLEENGINE::Velocity, dthalf_, PARTICLEENGINE::Acceleration);

    // angular velocity and acceleration states
    if (container->HaveStoredState(PARTICLEENGINE::AngularVelocity) and
        container->HaveStoredState(PARTICLEENGINE::AngularAcceleration))
    {
      // update angular velocity of all particles
      container->UpdateState(
          1.0, PARTICLEENGINE::AngularVelocity, dthalf_, PARTICLEENGINE::AngularAcceleration);
    }
  }

  // evaluate dirichlet boundary condition
  if (dirichletboundarycondition_)
    dirichletboundarycondition_->EvaluateDirichletBoundaryCondition(time_, false, true, false);
}
