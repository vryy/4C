/*---------------------------------------------------------------------------*/
/*!
\brief class holding all equation of state handlers

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_sph_equationofstate_bundle.H"

#include "particle_interaction_sph_equationofstate.H"
#include "particle_interaction_material_handler.H"

#include "../drt_inpar/inpar_particle.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::SPHEquationOfStateBundle::SPHEquationOfStateBundle(
    const Teuchos::ParameterList& params)
    : params_sph_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init equation of state bundle                              sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHEquationOfStateBundle::Init(
    const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial)
{
  // get type of smoothed particle hydrodynamics equation of state
  INPAR::PARTICLE::EquationOfStateType equationofstatetype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::EquationOfStateType>(
          params_sph_, "EQUATIONOFSTATE");

  // determine size of vector indexed by particle types
  const int typevectorsize = *(--particlematerial->GetParticleTypes().end()) + 1;

  // allocate memory to hold particle types
  phasetypetoequationofstate_.resize(typevectorsize);

  // iterate over particle types
  for (const auto& typeEnum : particlematerial->GetParticleTypes())
  {
    // no equation of state for boundary or rigid particles
    if (typeEnum == PARTICLEENGINE::BoundaryPhase or typeEnum == PARTICLEENGINE::RigidPhase)
      continue;

    // add to set of particle types of stored equation of state handlers
    storedtypes_.insert(typeEnum);

    // get material for current particle type
    const MAT::PAR::ParticleMaterialSPHFluid* material =
        dynamic_cast<const MAT::PAR::ParticleMaterialSPHFluid*>(
            particlematerial->GetPtrToParticleMatParameter(typeEnum));

    // create equation of state handler
    switch (equationofstatetype)
    {
      case INPAR::PARTICLE::GenTait:
      {
        const double speedofsound = material->SpeedOfSound();
        const double refdensfac = material->refDensFac_;
        const double exponent = material->exponent_;

        phasetypetoequationofstate_[typeEnum] =
            std::unique_ptr<PARTICLEINTERACTION::SPHEquationOfStateGenTait>(
                new PARTICLEINTERACTION::SPHEquationOfStateGenTait(
                    speedofsound, refdensfac, exponent));
        break;
      }
      case INPAR::PARTICLE::IdealGas:
      {
        const double speedofsound = material->SpeedOfSound();

        phasetypetoequationofstate_[typeEnum] =
            std::unique_ptr<PARTICLEINTERACTION::SPHEquationOfStateIdealGas>(
                new PARTICLEINTERACTION::SPHEquationOfStateIdealGas(speedofsound));
        break;
      }
      default:
      {
        dserror("unknown equation of state type!");
        break;
      }
    }

    // init equation of state handler
    phasetypetoequationofstate_[typeEnum]->Init();
  }
}

/*---------------------------------------------------------------------------*
 | setup equation of state bundle                             sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHEquationOfStateBundle::Setup()
{
  for (PARTICLEENGINE::TypeEnum typeEnum : storedtypes_)
    phasetypetoequationofstate_[typeEnum]->Setup();
}

/*---------------------------------------------------------------------------*
 | write restart of equation of state bundle                  sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHEquationOfStateBundle::WriteRestart(
    const int step, const double time) const
{
  for (PARTICLEENGINE::TypeEnum typeEnum : storedtypes_)
    phasetypetoequationofstate_[typeEnum]->WriteRestart(step, time);
}

/*---------------------------------------------------------------------------*
 | read restart of equation of state bundle                   sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHEquationOfStateBundle::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  for (PARTICLEENGINE::TypeEnum typeEnum : storedtypes_)
    phasetypetoequationofstate_[typeEnum]->ReadRestart(reader);
}
