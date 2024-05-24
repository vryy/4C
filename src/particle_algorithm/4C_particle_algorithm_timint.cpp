/*---------------------------------------------------------------------------*/
/*! \file
\brief time integration for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_algorithm_timint.hpp"

#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_particle_algorithm_dirichlet_bc.hpp"
#include "4C_particle_algorithm_temperature_bc.hpp"
#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_container_bundle.hpp"
#include "4C_particle_engine_enums.hpp"
#include "4C_particle_engine_interface.hpp"
#include "4C_particle_rigidbody_interface.hpp"
#include "4C_utils_parameter_list.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

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
  init_dirichlet_boundary_condition();

  // init temperature boundary condition handler
  init_temperature_boundary_condition();
}

void PARTICLEALGORITHM::TimInt::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLERIGIDBODY::RigidBodyHandlerInterface> particlerigidbodyinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set interface to rigid body handler
  particlerigidbodyinterface_ = particlerigidbodyinterface;

  // setup dirichlet boundary condition handler
  if (dirichletboundarycondition_) dirichletboundarycondition_->Setup(particleengineinterface_);

  // setup temperature boundary condition handler
  if (temperatureboundarycondition_) temperatureboundarycondition_->Setup(particleengineinterface_);

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  // set of particle types being excluded from time integration
  std::set<PARTICLEENGINE::TypeEnum> typesexludedfromtimeintegration;

  // boundary and rigid particles are not integrated in time
  typesexludedfromtimeintegration.insert(
      {PARTICLEENGINE::BoundaryPhase, PARTICLEENGINE::RigidPhase});

  if (dirichletboundarycondition_)
  {
    // get reference to set of particle types subjected to dirichlet boundary conditions
    const std::set<PARTICLEENGINE::TypeEnum>& typessubjectedtodirichletbc =
        dirichletboundarycondition_->get_particle_types_subjected_to_dirichlet_bc_set();

    // particles subjected to dirichlet boundary conditions are not integrated in time
    for (PARTICLEENGINE::TypeEnum currtype : typessubjectedtodirichletbc)
      typesexludedfromtimeintegration.insert(currtype);
  }

  // determine set of particle types to be integrated in time
  for (auto& typeEnum : particlecontainerbundle->GetParticleTypes())
    if (not typesexludedfromtimeintegration.count(typeEnum)) typestointegrate_.insert(typeEnum);
}

void PARTICLEALGORITHM::TimInt::insert_particle_states_of_particle_types(
    std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>& particlestatestotypes)
    const
{
  // insert dbc dependent states of all particle types
  if (dirichletboundarycondition_)
    dirichletboundarycondition_->insert_particle_states_of_particle_types(particlestatestotypes);

  // insert tempbc dependent states of all particle types
  if (temperatureboundarycondition_)
    temperatureboundarycondition_->insert_particle_states_of_particle_types(particlestatestotypes);
}

void PARTICLEALGORITHM::TimInt::SetInitialStates()
{
  // add initial random noise to particle position
  add_initial_random_noise_to_position();

  // set particle reference position
  if (dirichletboundarycondition_) dirichletboundarycondition_->set_particle_reference_position();

  // set particle reference position
  if (temperatureboundarycondition_)
    temperatureboundarycondition_->set_particle_reference_position();

  // evaluate dirichlet boundary condition
  if (dirichletboundarycondition_)
    dirichletboundarycondition_->evaluate_dirichlet_boundary_condition(0.0, true, true, true);

  // evaluate temperature boundary condition
  if (temperatureboundarycondition_)
    temperatureboundarycondition_->evaluate_temperature_boundary_condition(0.0);
}

void PARTICLEALGORITHM::TimInt::set_current_time(const double currenttime) { time_ = currenttime; }

void PARTICLEALGORITHM::TimInt::init_dirichlet_boundary_condition()
{
  // create dirichlet boundary condition handler
  dirichletboundarycondition_ =
      std::unique_ptr<PARTICLEALGORITHM::DirichletBoundaryConditionHandler>(
          new PARTICLEALGORITHM::DirichletBoundaryConditionHandler(params_));

  // init dirichlet boundary condition handler
  dirichletboundarycondition_->Init();

  // get reference to set of particle types subjected to dirichlet boundary conditions
  const std::set<PARTICLEENGINE::TypeEnum>& typessubjectedtodirichletbc =
      dirichletboundarycondition_->get_particle_types_subjected_to_dirichlet_bc_set();

  // no particle types are subjected to dirichlet boundary conditions
  if (typessubjectedtodirichletbc.empty()) dirichletboundarycondition_.reset();
}

void PARTICLEALGORITHM::TimInt::init_temperature_boundary_condition()
{
  // create temperature boundary condition handler
  temperatureboundarycondition_ =
      std::unique_ptr<PARTICLEALGORITHM::TemperatureBoundaryConditionHandler>(
          new PARTICLEALGORITHM::TemperatureBoundaryConditionHandler(params_));

  // init temperature boundary condition handler
  temperatureboundarycondition_->Init();

  // get reference to set of particle types subjected to temperature boundary conditions
  const std::set<PARTICLEENGINE::TypeEnum>& typessubjectedtotempbc =
      temperatureboundarycondition_->get_particle_types_subjected_to_temperature_bc_set();

  // no particle types are subjected to temperature boundary conditions
  if (typessubjectedtotempbc.empty()) temperatureboundarycondition_.reset();
}

void PARTICLEALGORITHM::TimInt::add_initial_random_noise_to_position()
{
  // init vector of initial position amplitude for each spatial direction
  std::vector<double> amplitude;
  double value;
  std::istringstream amplitudestream(
      Teuchos::getNumericStringParameter(params_, "INITIAL_POSITION_AMPLITUDE"));

  while (amplitudestream >> value) amplitude.push_back(value);

  // safety check
  if (static_cast<int>(amplitude.size()) != 3)
    FOUR_C_THROW("dimension (dim = %d) of initial position amplitude vector is wrong!",
        static_cast<int>(amplitude.size()));

  // safety check
  for (double a : amplitude)
    if (a < 0.0)
      FOUR_C_THROW(
          "no negative initial position amplitude allowed (set a positive or zero value)!");

  // get magnitude of initial position amplitude
  double temp = 0.0;
  for (double a : amplitude) temp += a * a;
  const double amplitude_norm = std::sqrt(temp);

  // no initial position amplitude defined
  if (not(amplitude_norm > 0.0)) return;

  // safety check
  const double max_amplitude = *std::max_element(amplitude.begin(), amplitude.end());
  if (max_amplitude > particleengineinterface_->MinBinSize())
    FOUR_C_THROW(
        "amplitude of noise added to initial position larger than minimum relevant bin size!");

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  // iterate over particle types
  for (auto& particleType : typestointegrate_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->get_specific_container(particleType, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get pointer to particle state
    double* pos = container->GetPtrToState(PARTICLEENGINE::Position, 0);

    // get particle state dimension
    int statedim = container->GetStateDim(PARTICLEENGINE::Position);

    // iterate over owned particles of current type
    for (int i = 0; i < particlestored; ++i)
    {
      // iterate over spatial dimension
      for (int dim = 0; dim < statedim; ++dim)
      {
        // generate random value
        const double randomvalue = GLOBAL::Problem::Instance()->Random()->Uni();

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
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLERIGIDBODY::RigidBodyHandlerInterface> particlerigidbodyinterface)
{
  // call base class setup
  PARTICLEALGORITHM::TimInt::Setup(particleengineinterface, particlerigidbodyinterface);

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  // iterate over particle types
  for (auto& particleType : typestointegrate_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->get_specific_container(particleType, PARTICLEENGINE::Owned);

    // safety check
    if (container->HaveStoredState(PARTICLEENGINE::ModifiedVelocity) or
        container->HaveStoredState(PARTICLEENGINE::ModifiedAcceleration))
      FOUR_C_THROW(
          "modified velocity and acceleration states not implemented yet for semi-implicit Euler "
          "time integration scheme!");
  }
}

void PARTICLEALGORITHM::TimIntSemiImplicitEuler::pre_interaction_routine()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::TimIntSemiImplicitEuler::pre_interaction_routine");

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  // iterate over particle types
  for (auto& particleType : typestointegrate_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->get_specific_container(particleType, PARTICLEENGINE::Owned);

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

  if (particlerigidbodyinterface_)
  {
    // update velocities of all rigid bodies and corresponding rigid particles
    particlerigidbodyinterface_->UpdateVelocities(dt_);

    // clear accelerations of all rigid bodies
    particlerigidbodyinterface_->ClearAccelerations();

    // update positions of all rigid bodies and corresponding rigid particles
    particlerigidbodyinterface_->UpdatePositions(dt_);
  }

  // evaluate dirichlet boundary condition
  if (dirichletboundarycondition_)
    dirichletboundarycondition_->evaluate_dirichlet_boundary_condition(time_, true, true, true);

  // evaluate temperature boundary condition
  if (temperatureboundarycondition_)
    temperatureboundarycondition_->evaluate_temperature_boundary_condition(time_);
}

void PARTICLEALGORITHM::TimIntSemiImplicitEuler::post_interaction_routine()
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
      particleengineinterface_->get_particle_container_bundle();

  // iterate over particle types
  for (auto& particleType : typestointegrate_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->get_specific_container(particleType, PARTICLEENGINE::Owned);

    // modified velocity states
    if (container->HaveStoredState(PARTICLEENGINE::ModifiedVelocity))
    {
      // update modified velocity of all particles
      container->UpdateState(0.0, PARTICLEENGINE::ModifiedVelocity, 1.0, PARTICLEENGINE::Velocity);
    }
  }
}

void PARTICLEALGORITHM::TimIntVelocityVerlet::pre_interaction_routine()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::TimIntVelocityVerlet::pre_interaction_routine");

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  // iterate over particle types
  for (auto& particleType : typestointegrate_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->get_specific_container(particleType, PARTICLEENGINE::Owned);

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

  if (particlerigidbodyinterface_)
  {
    // update velocities of all rigid bodies and corresponding rigid particles
    particlerigidbodyinterface_->UpdateVelocities(dthalf_);

    // clear accelerations of all rigid bodies
    particlerigidbodyinterface_->ClearAccelerations();

    // update positions of all rigid bodies and corresponding rigid particles
    particlerigidbodyinterface_->UpdatePositions(dt_);
  }

  // evaluate dirichlet boundary condition
  if (dirichletboundarycondition_)
  {
    dirichletboundarycondition_->evaluate_dirichlet_boundary_condition(time_, true, false, true);
    dirichletboundarycondition_->evaluate_dirichlet_boundary_condition(
        time_ - dthalf_, false, true, false);
  }

  // evaluate temperature boundary condition
  if (temperatureboundarycondition_)
    temperatureboundarycondition_->evaluate_temperature_boundary_condition(time_);
}

void PARTICLEALGORITHM::TimIntVelocityVerlet::post_interaction_routine()
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEALGORITHM::TimIntVelocityVerlet::post_interaction_routine");

  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  // iterate over particle types
  for (auto& particleType : typestointegrate_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->get_specific_container(particleType, PARTICLEENGINE::Owned);

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

  // update velocities of all rigid bodies and corresponding rigid particles
  if (particlerigidbodyinterface_) particlerigidbodyinterface_->UpdateVelocities(dthalf_);

  // evaluate dirichlet boundary condition
  if (dirichletboundarycondition_)
    dirichletboundarycondition_->evaluate_dirichlet_boundary_condition(time_, false, true, false);
}

FOUR_C_NAMESPACE_CLOSE
