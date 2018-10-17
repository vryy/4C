/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_pressure.cpp

\brief pressure handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_pressure.H"

#include "particle_interaction_material_handler.H"
#include "particle_interaction_sph_equationofstate.H"
#include "particle_interaction_sph_equationofstate_bundle.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHPressure::SPHPressure()
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init pressure handler                                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHPressure::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup pressure handler                                     sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHPressure::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
    const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set particle material handler
  particlematerial_ = particlematerial;

  // set equation of state handler
  equationofstatebundle_ = equationofstatebundle;
}

/*---------------------------------------------------------------------------*
 | write restart of pressure handler                          sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHPressure::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of pressure handler                           sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHPressure::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | compute pressure using equation of state and density       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHPressure::ComputePressure() const
{
  // iterate over particle types
  for (auto& typeIt : particlecontainerbundle_->GetRefToAllContainersMap())
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // no pressure computation for boundary particles
    if (type == PARTICLEENGINE::BoundaryPhase) continue;

    // get container of owned particles of current particle type
    auto statusIt = (typeIt.second).find(PARTICLEENGINE::Owned);
    if (statusIt == (typeIt.second).end())
      dserror("particle status '%s' not found!",
          PARTICLEENGINE::EnumToStatusName(PARTICLEENGINE::Owned).c_str());
    PARTICLEENGINE::ParticleContainerShrdPtr container = statusIt->second;

    // get number of particles stored in container
    int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get pointer to particle state
    const double* dens = container->GetPtrToParticleState(PARTICLEENGINE::Density, 0);
    double* press = container->GetPtrToParticleState(PARTICLEENGINE::Pressure, 0);

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material =
        particlematerial_->GetPtrToParticleMatParameter(type);

    // get equation of state for current particle type
    std::shared_ptr<SPHEquationOfStateBase> equationofstate =
        equationofstatebundle_->GetSpecificEquationOfState(type);

    // iterate over owned particles of current type
    for (int i = 0; i < particlestored; ++i)
      press[i] = equationofstate->DensityToPressure(dens[i], material->initDensity_);
  }

  // refresh pressure of ghosted particles
  RefreshPressure();
}

/*---------------------------------------------------------------------------*
 | refresh pressure of ghosted particles                      sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHPressure::RefreshPressure() const
{
  // init map
  std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>> particlestatestotypes;

  // iterate over particle types
  for (auto& typeIt : particlecontainerbundle_->GetRefToAllContainersMap())
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum type = typeIt.first;

    // no refreshing of pressure states for boundary particles
    if (type == PARTICLEENGINE::BoundaryPhase) continue;

    // set state enums to map
    particlestatestotypes[type].insert(PARTICLEENGINE::Pressure);
  }

  // refresh specific states of particles of specific types
  particleengineinterface_->RefreshSpecificStatesOfParticlesOfSpecificTypes(particlestatestotypes);
}
