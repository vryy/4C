/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_material_handler.cpp

\brief particle material handler for particle simulations

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
#include "particle_interaction_material_handler.H"

#include "../drt_particle_algorithm/particle_algorithm_utils.H"

#include "../drt_mat/matpar_bundle.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | constructor                                                sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::MaterialHandler::MaterialHandler(const Teuchos::ParameterList& params)
    : params_(params)
{
  // empty constructor
}

/*---------------------------------------------------------------------------*
 | init particle material handler                             sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::MaterialHandler::Init()
{
  // init map relating particle types to material ids
  std::map<PARTICLEENGINE::TypeEnum, int> typetomatidmap;

  // read parameters relating particle types to IDs
  PARTICLEALGORITHM::UTILS::ReadParamsTypesRelatedToIDs(
      params_, "PHASE_TO_MATERIAL_ID", typetomatidmap);

  // map particle types to particle material parameters
  for (auto& typeIt : typetomatidmap) MapParticleTypeToMatParameter(typeIt.first, typeIt.second);
}

/*---------------------------------------------------------------------------*
 | setup particle material handler                            sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::MaterialHandler::Setup()
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | write restart of particle material handler                 sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::MaterialHandler::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | read restart of particle material handler                  sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::MaterialHandler::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}

/*---------------------------------------------------------------------------*
 | return pointer to particle material parameter              sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
const MAT::PAR::ParticleMaterialBase*
PARTICLEINTERACTION::MaterialHandler::GetPtrToParticleMatParameter(
    PARTICLEENGINE::TypeEnum particleType) const
{
  auto typeIt = phasetypetoparticlematpar_.find(particleType);
  if (typeIt == phasetypetoparticlematpar_.end())
    dserror("particle material parameters of phase '%s' not found!",
        PARTICLEENGINE::EnumToTypeName(particleType).c_str());
  return typeIt->second;
}

/*---------------------------------------------------------------------------*
 | map particle types to particle material parameters         sfuchs 07/2018 |
 *---------------------------------------------------------------------------*/
void PARTICLEINTERACTION::MaterialHandler::MapParticleTypeToMatParameter(
    PARTICLEENGINE::TypeEnum particleType, int matID)
{
  // get material parameters and cast to particle material parameter
  const MAT::PAR::Parameter* matparameter =
      DRT::Problem::Instance()->Materials()->ParameterById(matID);
  const MAT::PAR::ParticleMaterialBase* particlematparameter =
      dynamic_cast<const MAT::PAR::ParticleMaterialBase*>(matparameter);

  // safety check
  if (particlematparameter == NULL) dserror("cast to specific particle material failed!");

  phasetypetoparticlematpar_.insert(std::make_pair(particleType, particlematparameter));
}
