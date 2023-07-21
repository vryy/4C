/*---------------------------------------------------------------------------*/
/*! \file
\brief evaporation induced recoil pressure handler for smoothed particle hydrodynamics (SPH)
       interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_particle_interaction_sph_surface_tension_recoilpressure_evaporation.H"

#include "baci_particle_interaction_utils.H"

#include "baci_particle_engine_interface.H"
#include "baci_particle_engine_container.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHRecoilPressureEvaporation::SPHRecoilPressureEvaporation(
    const Teuchos::ParameterList& params)
    : params_sph_(params),
      evaporatingphase_(PARTICLEENGINE::Phase1),
      recoilboilingtemp_(params_sph_.get<double>("VAPOR_RECOIL_BOILINGTEMPERATURE")),
      recoil_pfac_(params_sph_.get<double>("VAPOR_RECOIL_PFAC")),
      recoil_tfac_(params_sph_.get<double>("VAPOR_RECOIL_TFAC"))
{
  // empty constructor
}

void PARTICLEINTERACTION::SPHRecoilPressureEvaporation::Init()
{
  // nothing to do
}

void PARTICLEINTERACTION::SPHRecoilPressureEvaporation::Setup(
    const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface)
{
  // set interface to particle engine
  particleengineinterface_ = particleengineinterface;

  // set particle container bundle
  particlecontainerbundle_ = particleengineinterface_->GetParticleContainerBundle();
}

void PARTICLEINTERACTION::SPHRecoilPressureEvaporation::ComputeRecoilPressureContribution() const
{
  // get container of owned particles of evaporating phase
  PARTICLEENGINE::ParticleContainer* container_i =
      particlecontainerbundle_->GetSpecificContainer(evaporatingphase_, PARTICLEENGINE::Owned);

  // iterate over particles in container
  for (int particle_i = 0; particle_i < container_i->ParticlesStored(); ++particle_i)
  {
    const double* dens_i = container_i->GetPtrToState(PARTICLEENGINE::Density, particle_i);
    const double* temp_i = container_i->GetPtrToState(PARTICLEENGINE::Temperature, particle_i);
    const double* cfg_i =
        container_i->GetPtrToState(PARTICLEENGINE::ColorfieldGradient, particle_i);
    const double* ifn_i = container_i->GetPtrToState(PARTICLEENGINE::InterfaceNormal, particle_i);
    double* acc_i = container_i->GetPtrToState(PARTICLEENGINE::Acceleration, particle_i);

    // evaluation only for non-zero interface normal
    if (not(UTILS::VecNormTwo(ifn_i) > 0.0)) continue;

    // recoil pressure contribution only for temperature above boiling temperature
    if (not(temp_i[0] > recoilboilingtemp_)) continue;

    // compute evaporation induced recoil pressure
    const double recoil_press_i =
        recoil_pfac_ * std::exp(-recoil_tfac_ * (1.0 / temp_i[0] - 1.0 / recoilboilingtemp_));

    // add contribution to acceleration
    UTILS::VecAddScale(acc_i, -recoil_press_i / dens_i[0], cfg_i);
  }
}
