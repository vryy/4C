/*---------------------------------------------------------------------------*/
/*! \file
\brief class holding all equation of state handlers
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "4C_particle_interaction_sph_equationofstate_bundle.hpp"

#include "4C_inpar_particle.hpp"
#include "4C_particle_interaction_material_handler.hpp"
#include "4C_particle_interaction_sph_equationofstate.hpp"
#include "4C_utils_exceptions.hpp"

FOUR_C_NAMESPACE_OPEN

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
ParticleInteraction::SPHEquationOfStateBundle::SPHEquationOfStateBundle(
    const Teuchos::ParameterList& params)
    : params_sph_(params)
{
  // empty constructor
}

void ParticleInteraction::SPHEquationOfStateBundle::Init(
    const std::shared_ptr<ParticleInteraction::MaterialHandler> particlematerial)
{
  // get type of smoothed particle hydrodynamics equation of state
  Inpar::PARTICLE::EquationOfStateType equationofstatetype =
      Core::UTILS::IntegralValue<Inpar::PARTICLE::EquationOfStateType>(
          params_sph_, "EQUATIONOFSTATE");

  // determine size of vector indexed by particle types
  const int typevectorsize = *(--particlematerial->GetParticleTypes().end()) + 1;

  // allocate memory to hold particle types
  phasetypetoequationofstate_.resize(typevectorsize);

  // iterate over particle types
  for (const auto& type_i : particlematerial->GetParticleTypes())
  {
    // no equation of state for boundary or rigid particles
    if (type_i == PARTICLEENGINE::BoundaryPhase or type_i == PARTICLEENGINE::RigidPhase) continue;

    // add to set of particle types of stored equation of state handlers
    storedtypes_.insert(type_i);

    // get material for current particle type
    const Mat::PAR::ParticleMaterialSPHFluid* material =
        dynamic_cast<const Mat::PAR::ParticleMaterialSPHFluid*>(
            particlematerial->get_ptr_to_particle_mat_parameter(type_i));

    // create equation of state handler
    switch (equationofstatetype)
    {
      case Inpar::PARTICLE::GenTait:
      {
        const double speedofsound = material->SpeedOfSound();
        const double refdensfac = material->refDensFac_;
        const double exponent = material->exponent_;

        phasetypetoequationofstate_[type_i] =
            std::unique_ptr<ParticleInteraction::SPHEquationOfStateGenTait>(
                new ParticleInteraction::SPHEquationOfStateGenTait(
                    speedofsound, refdensfac, exponent));
        break;
      }
      case Inpar::PARTICLE::IdealGas:
      {
        const double speedofsound = material->SpeedOfSound();

        phasetypetoequationofstate_[type_i] =
            std::unique_ptr<ParticleInteraction::SPHEquationOfStateIdealGas>(
                new ParticleInteraction::SPHEquationOfStateIdealGas(speedofsound));
        break;
      }
      default:
      {
        FOUR_C_THROW("unknown equation of state type!");
        break;
      }
    }

    // init equation of state handler
    phasetypetoequationofstate_[type_i]->Init();
  }
}

void ParticleInteraction::SPHEquationOfStateBundle::Setup()
{
  for (PARTICLEENGINE::TypeEnum type_i : storedtypes_) phasetypetoequationofstate_[type_i]->Setup();
}

FOUR_C_NAMESPACE_CLOSE
