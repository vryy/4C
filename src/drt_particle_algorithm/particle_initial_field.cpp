/*---------------------------------------------------------------------------*/
/*!
\brief initial field handler for particle simulations

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_initial_field.H"

#include "particle_algorithm_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_container_bundle.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_globalproblem.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::InitialFieldHandler::InitialFieldHandler(const Teuchos::ParameterList& params)
    : params_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init initial field handler                                 sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::InitialFieldHandler::Init()
{
  // get control parameters for initial/boundary conditions
  const Teuchos::ParameterList& params_conditions =
      params_.sublist("INITIAL AND BOUNDARY CONDITIONS");

  // relate particle state to input name
  std::map<std::string, PARTICLEENGINE::StateEnum> initialfieldtostateenum = {
      std::make_pair("INITIAL_TEMP_FIELD", PARTICLEENGINE::Temperature),
      std::make_pair("INITIAL_VELOCITY_FIELD", PARTICLEENGINE::Velocity),
      std::make_pair("INITIAL_ACCELERATION_FIELD", PARTICLEENGINE::Acceleration)};

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

/*---------------------------------------------------------------------------*
 | setup initial field handler                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::InitialFieldHandler::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

/*---------------------------------------------------------------------------*
 | set initial fields                                         sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::InitialFieldHandler::SetInitialFields()
{
  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

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
          particlecontainerbundle->GetSpecificContainer(particleType, PARTICLEENGINE::Owned);

      // get number of particles stored in container
      const int particlestored = container->ParticlesStored();

      // no owned particles of current particle type
      if (particlestored <= 0) continue;

      // get id of function
      const int functid = initialFieldIt.second;

      // get reference to function
      DRT::UTILS::Function& function = DRT::Problem::Instance()->Funct(functid - 1);

      // get pointer to particle position
      const double* pos = container->GetPtrToParticleState(PARTICLEENGINE::Position, 0);

      // get particle state dimension
      int posstatedim = container->GetParticleStateDim(PARTICLEENGINE::Position);

      // get pointer to particle state
      double* state = container->GetPtrToParticleState(particleState, 0);

#ifdef DEBUG
      if (not container->GetStoredStates().count(particleState))
        dserror("particle state '%s' not stored in container!",
            PARTICLEENGINE::EnumToStateName(particleState).c_str());
#endif

      // get particle state dimension
      int statedim = container->GetParticleStateDim(particleState);

      // safety check
      if (statedim != function.NumberComponents())
        dserror("dimensions of function defining initial field and of state '%s' not matching!",
            PARTICLEENGINE::EnumToStateName(particleState).c_str());

      // iterate over owned particles of current type
      for (int i = 0; i < particlestored; ++i)
      {
        // evaluate function to set initial field
        for (int dim = 0; dim < statedim; ++dim)
          state[statedim * i + dim] = function.Evaluate(dim, &(pos[posstatedim * i]), 0.0);
      }
    }
  }
}
