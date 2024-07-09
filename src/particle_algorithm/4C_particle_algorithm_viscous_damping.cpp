/*---------------------------------------------------------------------------*/
/*! \file
\brief viscous damping handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_algorithm_viscous_damping.hpp"

#include "4C_particle_engine_container.hpp"
#include "4C_particle_engine_container_bundle.hpp"
#include "4C_particle_engine_enums.hpp"
#include "4C_particle_engine_interface.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::ViscousDampingHandler::ViscousDampingHandler(const double viscdampfac)
    : viscdampfac_(viscdampfac)
{
  // empty constructor
}

void PARTICLEALGORITHM::ViscousDampingHandler::init()
{
  // nothing to do
}

void PARTICLEALGORITHM::ViscousDampingHandler::setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

void PARTICLEALGORITHM::ViscousDampingHandler::apply_viscous_damping()
{
  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->get_particle_container_bundle();

  // iterate over particle types
  for (const auto& type : particlecontainerbundle->get_particle_types())
  {
    // no viscous damping contribution for boundary or rigid particles
    if (type == PARTICLEENGINE::BoundaryPhase or type == PARTICLEENGINE::RigidPhase) continue;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->get_specific_container(type, PARTICLEENGINE::Owned);

    // apply viscous damping contribution
    container->update_state(
        1.0, PARTICLEENGINE::Acceleration, -viscdampfac_, PARTICLEENGINE::Velocity);
  }
}

FOUR_C_NAMESPACE_CLOSE
