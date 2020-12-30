/*---------------------------------------------------------------------------*/
/*! \file
\brief evaporation induced heat loss handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_heatloss_evaporation.H"

#include "particle_interaction_material_handler.H"

#include "particle_interaction_utils.H"

#include "../drt_particle_engine/particle_engine_interface.H"
#include "../drt_particle_engine/particle_container.H"

#include "../drt_lib/drt_dserror.H"

#include <Teuchos_TimeMonitor.hpp>

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHHeatLossEvaporation::SPHHeatLossEvaporation(
    const Teuchos::ParameterList& params)
    : params_sph_(params),
      evaporatingphase_(PARTICLEENGINE::Phase1),
      recoilboilingtemp_(params_sph_.get<double>("VAPOR_RECOIL_BOILINGTEMPERATURE")),
      recoil_pfac_(params_sph_.get<double>("VAPOR_RECOIL_PFAC")),
      recoil_tfac_(params_sph_.get<double>("VAPOR_RECOIL_TFAC")),
      latentheat_(params_sph_.get<double>("VAPOR_HEATLOSS_LATENTHEAT")),
      enthalpyreftemp_(params_sph_.get<double>("VAPOR_HEATLOSS_ENTHALPY_REFTEMP")),
      heatloss_pfac_(params_sph_.get<double>("VAPOR_HEATLOSS_PFAC")),
      heatloss_tfac_(params_sph_.get<double>("VAPOR_HEATLOSS_TFAC"))
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHHeatLossEvaporation::Init()
{
  // safety check
  if (DRT::INPUT::IntegralValue<INPAR::PARTICLE::SurfaceTensionFormulation>(
          params_sph_, "SURFACETENSIONFORMULATION") == INPAR::PARTICLE::NoSurfaceTension)
    dserror("surface tension evaluation needed for evaporation induced heat loss!");
}

void PARTICLEINTERACTION::SPHHeatLossEvaporation::Setup(
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
    thermomaterial_[type_i] = dynamic_cast<const MAT::PAR::ParticleMaterialThermo*>(
        particlematerial_->GetPtrToParticleMatParameter(type_i));
}

void PARTICLEINTERACTION::SPHHeatLossEvaporation::EvaluateEvaporationInducedHeatLoss() const
{
  TEUCHOS_FUNC_TIME_MONITOR(
      "PARTICLEINTERACTION::SPHHeatLossEvaporation::EvaluateEvaporationInducedHeatLoss");

  // get container of owned particles of evaporating phase
  PARTICLEENGINE::ParticleContainer* container_i =
      particlecontainerbundle_->GetSpecificContainer(evaporatingphase_, PARTICLEENGINE::Owned);

  const MAT::PAR::ParticleMaterialThermo* thermomaterial_i = thermomaterial_[evaporatingphase_];

  // iterate over particles in container
  for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
  {
    const double* dens_i = container_i->GetPtrToState(PARTICLEENGINE::Density, particle_i);
    const double* temp_i = container_i->GetPtrToState(PARTICLEENGINE::Temperature, particle_i);
    const double* cfg_i =
        container_i->GetPtrToState(PARTICLEENGINE::ColorfieldGradient, particle_i);
    const double* ifn_i = container_i->GetPtrToState(PARTICLEENGINE::InterfaceNormal, particle_i);
    double* tempdot_i = container_i->GetPtrToState(PARTICLEENGINE::TemperatureDot, particle_i);

    // evaluation only for non-zero interface normal
    if (not(UTILS::vec_norm2(ifn_i) > 0.0)) continue;

    // heat loss contribution only for temperature above boiling temperature
    if (not(temp_i[0] > recoilboilingtemp_)) continue;

    // compute evaporation induced recoil pressure
    const double recoil_press_i =
        recoil_pfac_ * std::exp(-recoil_tfac_ * (1.0 / temp_i[0] - 1.0 / recoilboilingtemp_));

    // compute vapor mass flow
    const double m_dot_i = heatloss_pfac_ * recoil_press_i * std::sqrt(heatloss_tfac_ / temp_i[0]);

    // evaluate specific enthalpy
    const double specificenthalpy_i =
        thermomaterial_i->thermalCapacity_ * (temp_i[0] - enthalpyreftemp_);

    // add contribution of heat loss
    tempdot_i[0] -= UTILS::vec_norm2(cfg_i) * m_dot_i * (latentheat_ + specificenthalpy_i) *
                    thermomaterial_i->invThermalCapacity_ / dens_i[0];
  }
}
