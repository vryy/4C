/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_heatsource.cpp

\brief heat source handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_heatsource.H"

#include "particle_interaction_material_handler.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include "../drt_lib/drt_globalproblem.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHHeatSourceBase::SPHHeatSourceBase(const Teuchos::ParameterList& params)
    : heatsourcefctnumber_(params.get<int>("HEATSOURCE_FUNCT"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init heat source handler                                   sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHHeatSourceBase::Init()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | setup heat source handler                                  sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHHeatSourceBase::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();

  // set particle material handler
  particlematerial_ = particlematerial;

  // determine size of vectors indexed by particle types
  const int typevectorsize = *(--particlecontainerbundle_->GetParticleTypes().end()) + 1;

  // allocate memory to hold particle types
  thermomaterial_.resize(typevectorsize);

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    thermomaterial_[type_i] = dynamic_cast<const MAT::PAR::ParticleMaterialThermo*>(
        particlematerial_->GetPtrToParticleMatParameter(type_i));

    // safety check
    if (not(thermomaterial_[type_i]->thermalAbsorptivity_ > 0.0))
      dserror("thermal absorptivity for particles of type '%s' not positive!",
          PARTICLEENGINE::EnumToTypeName(type_i).c_str());
  }
}

/*---------------------------------------------------------------------------*
 | write restart of heat source handler                       sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHHeatSourceBase::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of heat source handler                        sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHHeatSourceBase::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHHeatSourceVolume::SPHHeatSourceVolume(const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHHeatSourceBase(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | evaluate heat source                                       sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHHeatSourceVolume::EvaluateHeatSource(const double& evaltime) const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHHeatSourceVolume::EvaluateHeatSource");

  // init vector containing evaluated function
  std::vector<double> funct(1);

  // iterate over particle types
  for (const auto& type_i : particlecontainerbundle_->GetParticleTypes())
  {
    // no heat source evaluation for boundary particles
    if (type_i == PARTICLEENGINE::BoundaryPhase) continue;

    // get container of owned particles of current particle type
    PARTICLEENGINE::ParticleContainer* container_i =
        particlecontainerbundle_->GetSpecificContainer(type_i, PARTICLEENGINE::Owned);

    // get number of particles stored in container
    const int particlestored = container_i->ParticlesStored();

    // no owned particles of current particle type
    if (particlestored <= 0) continue;

    // get material for current particle type
    const MAT::PAR::ParticleMaterialThermo* thermomaterial_i = thermomaterial_[type_i];

    // get reference to function
    DRT::UTILS::Function& function = DRT::Problem::Instance()->Funct(heatsourcefctnumber_ - 1);

    // safety check
    if (function.NumberComponents() != 1)
      dserror("dimension of function defining heat source is not one!");

    // iterate over particles in container
    for (int particle_i = 0; particle_i < particlestored; ++particle_i)
    {
      // declare pointer variables for particle i
      const double* pos_i;
      double* tempdot_i;

      // get pointer to particle states
      pos_i = container_i->GetPtrToParticleState(PARTICLEENGINE::Position, particle_i);
      tempdot_i = container_i->GetPtrToParticleState(PARTICLEENGINE::TemperatureDot, particle_i);

      // evaluate function defining heat source
      funct = function.EvaluateTimeDerivative(0, &(pos_i[0]), evaltime, 0);

      // add contribution of heat source
      tempdot_i[0] +=
          thermomaterial_i->thermalAbsorptivity_ * funct[0] * thermomaterial_i->invThermalCapacity_;
    }
  }
}

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHHeatSourceSurface::SPHHeatSourceSurface(
    const Teuchos::ParameterList& params)
    : PARTICLEINTERACTION::SPHHeatSourceBase(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | evaluate heat source                                       sfuchs 02/2019 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHHeatSourceSurface::EvaluateHeatSource(const double& evaltime) const
{
  TEUCHOS_FUNC_TIME_MONITOR("PARTICLEINTERACTION::SPHHeatSourceSurface::EvaluateHeatSource");

  dserror("method 'EvaluateHeatSource(...)' not implemented yet!");
}
