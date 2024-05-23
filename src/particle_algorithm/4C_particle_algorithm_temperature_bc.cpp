/*---------------------------------------------------------------------------*/
/*! \file
\brief temperature boundary condition handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_algorithm_temperature_bc.hpp"

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
PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::TemperatureBoundaryConditionHandler(
    const Teuchos::ParameterList& params)
    : params_(params)
{
  // empty constructor
}

void PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::Init()
{
  // get control parameters for conditions
  const Teuchos::ParameterList& params_conditions =
      params_.sublist("INITIAL AND BOUNDARY CONDITIONS");

  // read parameters relating particle types to values
  PARTICLEALGORITHM::UTILS::ReadParamsTypesRelatedToValues(
      params_conditions, "TEMPERATURE_BOUNDARY_CONDITION", temperaturebctypetofunctid_);

  // iterate over particle types and insert into set
  for (auto& typeIt : temperaturebctypetofunctid_)
    typessubjectedtotemperaturebc_.insert(typeIt.first);
}

void PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

void PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::
    insert_particle_states_of_particle_types(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const
{
  // iterate over particle types subjected to temperature boundary conditions
  for (auto& particleType : typessubjectedtotemperaturebc_)
  {
    // insert states for types subjected to temperature boundary conditions
    particlestatestotypes[particleType].insert(PARTICLEENGINE::ReferencePosition);
  }
}

void PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::set_particle_reference_position() const
{
  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  // iterate over particle types subjected to temperature boundary conditions
  for (auto& particleType : typessubjectedtotemperaturebc_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->get_specific_container(particleType, PARTICLEENGINE::Owned);

    // set particle reference position
    container->UpdateState(0.0, PARTICLEENGINE::ReferencePosition, 1.0, PARTICLEENGINE::Position);
  }
}

void PARTICLEALGORITHM::TemperatureBoundaryConditionHandler::
    evaluate_temperature_boundary_condition(const double& evaltime) const
{
  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  // iterate over particle types subjected to temperature boundary conditions
  for (auto& typeIt : temperaturebctypetofunctid_)
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum particleType = typeIt.first;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->get_specific_container(particleType, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get id of function
    const int functid = typeIt.second;

    // get reference to function
    const auto& function =
        GLOBAL::Problem::Instance()->FunctionById<CORE::UTILS::FunctionOfSpaceTime>(functid - 1);

    // get pointer to particle states
    const double* refpos = container->GetPtrToState(PARTICLEENGINE::ReferencePosition, 0);
    double* temp = container->GetPtrToState(PARTICLEENGINE::Temperature, 0);

    // get particle state dimension
    int statedim = container->GetStateDim(PARTICLEENGINE::Position);

    // safety check
    if (function.NumberComponents() != 1)
      FOUR_C_THROW("dimension of function defining temperature boundary condition is not one!");

    // iterate over owned particles of current type
    for (int i = 0; i < particlestored; ++i)
    {
      // evaluate function
      temp[i] = function.Evaluate(&(refpos[statedim * i]), evaltime, 0);
    }
  }
}

FOUR_C_NAMESPACE_CLOSE
