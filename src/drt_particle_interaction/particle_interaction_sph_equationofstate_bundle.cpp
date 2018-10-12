/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_equationofstate_bundle.cpp

\brief class holding all equation of state handlers for smoothed particle hydrodynamics (SPH)
interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

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
    : params_sph_(params.sublist("SPH"))
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init equation of state bundle                              sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHEquationOfStateBundle::Init(
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial)
{
  // get type of smoothed particle hydrodynamics equation of state
  INPAR::PARTICLE::EquationOfStateType equationofstatetype =
      DRT::INPUT::IntegralValue<INPAR::PARTICLE::EquationOfStateType>(
          params_sph_, "EQUATIONOFSTATE");

  // iterate over particle types
  for (auto& typeIt : particlematerial->GetRefToParticleMatParMap())
  {
    // get type of particles
    PARTICLEENGINE::TypeEnum particleType = typeIt.first;

    // no equation of state for boundary particles
    if (particleType == PARTICLEENGINE::BoundaryPhase) continue;

    // get material for current particle type
    const MAT::PAR::ParticleMaterialSPHFluid* material =
        dynamic_cast<const MAT::PAR::ParticleMaterialSPHFluid*>(typeIt.second);

    // create equation of state handler
    switch (equationofstatetype)
    {
      case INPAR::PARTICLE::GenTait:
      {
        const double speedofsound = material->SpeedOfSound();
        const double refdensfac = material->refDensFac_;
        const double exponent = material->exponent_;

        phasetypetoequationofstate_[particleType] =
            std::make_shared<PARTICLEINTERACTION::SPHEquationOfStateGenTait>(
                speedofsound, refdensfac, exponent);
        break;
      }
      case INPAR::PARTICLE::IdealGas:
      {
        const double speedofsound = material->SpeedOfSound();

        phasetypetoequationofstate_[particleType] =
            std::make_shared<PARTICLEINTERACTION::SPHEquationOfStateIdealGas>(speedofsound);
        break;
      }
      default:
      {
        dserror("unknown equation of state type!");
        break;
      }
    }

    // init equation of state handler
    phasetypetoequationofstate_[particleType]->Init();
  }
}

/*---------------------------------------------------------------------------*
 | setup equation of state bundle                             sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHEquationOfStateBundle::Setup()
{
  for (auto& typeIt : phasetypetoequationofstate_) (typeIt.second)->Setup();
}

/*---------------------------------------------------------------------------*
 | write restart of equation of state bundle                  sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHEquationOfStateBundle::WriteRestart(
    const int step, const double time) const
{
  for (auto& typeIt : phasetypetoequationofstate_) (typeIt.second)->WriteRestart(step, time);
}

/*---------------------------------------------------------------------------*
 | read restart of equation of state bundle                   sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::SPHEquationOfStateBundle::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  for (auto& typeIt : phasetypetoequationofstate_) (typeIt.second)->ReadRestart(reader);
}

/*---------------------------------------------------------------------------*
 | get specific equation of state                             sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBase>
PARTICLEINTERACTION::SPHEquationOfStateBundle::GetSpecificEquationOfState(
    PARTICLEENGINE::TypeEnum particleType) const
{
  auto typeIt = phasetypetoequationofstate_.find(particleType);
  if (typeIt == phasetypetoequationofstate_.end())
    dserror("equation of state for particle type '%s' not found!",
        PARTICLEENGINE::EnumToTypeName(particleType).c_str());
  return typeIt->second;
}
