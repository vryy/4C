/*---------------------------------------------------------------------------*/
/*! \file
\brief pressure handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_pressure.H"

#include "particle_interaction_material_handler.H"
#include "particle_interaction_sph_equationofstate.H"
#include "particle_interaction_sph_equationofstate_bundle.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHPressure::SPHPressure()
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHPressure::Init()
{
  // init with potential fluid particle types
  fluidtypes_ = {PARTICLEENGINE::Phase1, PARTICLEENGINE::Phase2};
}

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

  // update with actual fluid particle types
  const auto fluidtypes = fluidtypes_;
  for (const auto& type_i : fluidtypes)
    if (not particlecontainerbundle_->GetParticleTypes().count(type_i)) fluidtypes_.erase(type_i);

  // setup pressure of ghosted particles to refresh
  {
    std::vector<PARTICLEENGINE::StateEnum> states{PARTICLEENGINE::Pressure};

    for (const auto& type_i : fluidtypes_)
      pressuretorefresh_.push_back(std::make_pair(type_i, states));
  }
}

void PARTICLEINTERACTION::SPHPressure::ComputePressure() const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHPressure::ComputePressure");

  // iterate over fluid particle types
  for (const auto& type_i : fluidtypes_)
  {
    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get pointer to particle state
    const double* dens = container->GetPtrToParticleState(PARTICLEENGINE::Density, 0);
    double* press = container->GetPtrToParticleState(PARTICLEENGINE::Pressure, 0);

    // get material for current particle type
    const MAT::PAR::ParticleMaterialBase* material =
        particlematerial_->GetPtrToParticleMatParameter(type_i);

    // get equation of state for current particle type
    const PARTICLEINTERACTION::SPHEquationOfStateBase* equationofstate =
        equationofstatebundle_->GetPtrToSpecificEquationOfState(type_i);

    // iterate over owned particles of current type
    for (int i = 0; i < particlestored; ++i)
      press[i] = equationofstate->DensityToPressure(dens[i], material->initDensity_);
  }

  // refresh pressure of ghosted particles
  particleengineinterface_->RefreshParticlesOfSpecificStatesAndTypes(pressuretorefresh_);
}
