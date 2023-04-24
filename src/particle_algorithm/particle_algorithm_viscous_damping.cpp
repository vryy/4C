/*---------------------------------------------------------------------------*/
/*! \file
\brief viscous damping handler for particle simulations
\level 2
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_algorithm_viscous_damping.H"

#include "particle_engine_interface.H"
#include "particle_engine_enums.H"
#include "particle_engine_container_bundle.H"
#include "particle_engine_container.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::ViscousDampingHandler::ViscousDampingHandler(const double viscdampfac)
    : viscdampfac_(viscdampfac)
{
  // empty constructor
}

void PARTICLEALGORITHM::ViscousDampingHandler::Init()
{
  // nothing to do
}

void PARTICLEALGORITHM::ViscousDampingHandler::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

void PARTICLEALGORITHM::ViscousDampingHandler::ApplyViscousDamping()
{
  // get particle container bundle
  PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle =
      particleengineinterface_->GetParticleContainerBundle();

  // iterate over particle types
  for (const auto& type : particlecontainerbundle->GetParticleTypes())
  {
    // no viscous damping contribution for boundary or rigid particles
    if (type == PARTICLEENGINE::BoundaryPhase or type == PARTICLEENGINE::RigidPhase) continue;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle->GetSpecificContainer(type, PARTICLEENGINE::Owned);

    // apply viscous damping contribution
    container->UpdateState(
        1.0, PARTICLEENGINE::Acceleration, -viscdampfac_, PARTICLEENGINE::Velocity);
  }
}
