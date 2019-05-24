/*---------------------------------------------------------------------------*/
/*!
\brief viscous damping handler for particle simulations

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
#include "particle_viscous_damping.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_container_bundle.H"
#include "../drt_particle_engine/particle_container.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEALGORITHM::ViscousDampingHandler::ViscousDampingHandler(const double viscdampfac)
    : viscdampfac_(viscdampfac)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init viscous damping handler                               sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ViscousDampingHandler::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup viscous damping handler                              sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ViscousDampingHandler::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;
}

/*---------------------------------------------------------------------------*
 | write restart of viscous damping handler                   sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ViscousDampingHandler::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of viscous damping handler                    sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEALGORITHM::ViscousDampingHandler::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | apply viscous damping contribution                         sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
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
