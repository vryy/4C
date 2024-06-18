/*---------------------------------------------------------------------------*/
/*! \file
\brief initial field handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_algorithm_initial_field.hpp"

#include "4C_global_data.hpp"
#include "4C_particle_algorithm_utils.hpp"
#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_container_bundle.hpp"
#include "4C_particle_engine_enums.hpp"
#include "4C_particle_engine_interface.hpp"
#include "4C_utils_function.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::InitialFieldHandler::InitialFieldHandler(const Teuchos::ParameterList& params)
    : params_(params)
{
  // empty constructor
}

void PARTICLEALGORITHM::InitialFieldHandler::init()
{
  // get control parameters for initial/boundary conditions
  const Teuchos::ParameterList& params_conditions =
      params_.sublist("INITIAL AND BOUNDARY CONDITIONS");

  // relate particle state to input name
  std::map<std::string, PARTICLEENGINE::StateEnum> initialfieldtostateenum = {
      std::make_pair("INITIAL_TEMP_FIELD", PARTICLEENGINE::Temperature),
      std::make_pair("INITIAL_VELOCITY_FIELD", PARTICLEENGINE::Velocity),
      std::make_pair("INITIAL_ANGULAR_VELOCITY_FIELD", PARTICLEENGINE::AngularVelocity),
      std::make_pair("INITIAL_ACCELERATION_FIELD", PARTICLEENGINE::Acceleration),
      std::make_pair("INITIAL_ANGULAR_ACCELERATION_FIELD", PARTICLEENGINE::AngularAcceleration)};

  // iterate over particle states
  for (auto& stateIt : initialfieldtostateenum)
  {
    // get reference to sub-map
    std::map<PARTICLEENGINE::TypeEnum, int>& currentstatetypetofunctidmap =
        statetotypetofunctidmap_[stateIt.second];

    // read parameters relating particle types to values
    PARTICLEALGORITHM::UTILS::ReadParamsTypesRelatedToValues(
        params_conditions, stateIt.first, currentstatetypetofunctidmap);
  }
}

void PARTICLEALGORITHM::InitialFieldHandler::setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

void PARTICLEALGORITHM::InitialFieldHandler::SetInitialFields()
{
  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  for (auto& stateIt : statetotypetofunctidmap_)
  {
    // get state of particles
    PARTICLEENGINE::StateEnum particleState = stateIt.first;

    // iterate over particle types
    for (auto& initialFieldIt : stateIt.second)
    {
      // get type of particles
      PARTICLEENGINE::TypeEnum particleType = initialFieldIt.first;

      // get container of owned particles of current particle type
      PARTICLEENGINE::ParticleContainer* container =
          particlecontainerbundle->get_specific_container(particleType, PARTICLEENGINE::Owned);

      // get number of particles stored in container
      const int particlestored = container->ParticlesStored();

      // no owned particles of current particle type
      if (particlestored <= 0) continue;

      if (not container->HaveStoredState(particleState)) continue;

      // get id of function
      const int functid = initialFieldIt.second;

      // get reference to function
      const auto& function =
          Global::Problem::Instance()->FunctionById<Core::UTILS::FunctionOfSpaceTime>(functid - 1);

      // get pointer to particle states
      const double* pos = container->GetPtrToState(PARTICLEENGINE::Position, 0);
      double* state = container->GetPtrToState(particleState, 0);

      // get particle state dimensions
      int posstatedim = container->GetStateDim(PARTICLEENGINE::Position);
      int statedim = container->GetStateDim(particleState);

      // safety check
      if (static_cast<std::size_t>(statedim) != function.NumberComponents())
        FOUR_C_THROW(
            "dimensions of function defining initial field and of state '%s' not matching!",
            PARTICLEENGINE::EnumToStateName(particleState).c_str());

      // iterate over owned particles of current type
      for (int i = 0; i < particlestored; ++i)
      {
        // evaluate function to set initial field
        for (int dim = 0; dim < statedim; ++dim)
          state[statedim * i + dim] = function.evaluate(&(pos[posstatedim * i]), 0.0, dim);
      }
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
