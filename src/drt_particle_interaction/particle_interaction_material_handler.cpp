/*---------------------------------------------------------------------------*/
/*! \file
\brief particle material handler for particle simulations
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_material_handler.H"

#include "../drt_particle_algorithm/particle_algorithm_utils.H"

#include "../drt_mat/matpar_bundle.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
PARTICLEINTERACTION::MaterialHandler::MaterialHandler(const Teuchos::ParameterList& params)
    : params_(params)
{
  // empty constructor
}

void PARTICLEINTERACTION::MaterialHandler::Init()
{
  // init map relating particle types to material ids
  std::map<PARTICLEENGINE::TypeEnum, int> typetomatidmap;

  // read parameters relating particle types to values
  PARTICLEALGORITHM::UTILS::ReadParamsTypesRelatedToValues(
      params_, "PHASE_TO_MATERIAL_ID", typetomatidmap);

  // determine size of vector indexed by particle types
  const int typevectorsize = ((--typetomatidmap.end())->first) + 1;

  // allocate memory to hold particle types
  phasetypetoparticlematpar_.resize(typevectorsize);

  // relate particle types to particle material parameters
  for (auto& typeIt : typetomatidmap)
  {
    // get type of particle
    PARTICLEENGINE::TypeEnum type_i = typeIt.first;

    // add to set of particle types of stored particle material parameters
    storedtypes_.insert(type_i);

    // get material parameters and cast to particle material parameter
    const MAT::PAR::Parameter* matparameter =
        DRT::Problem::Instance()->Materials()->ParameterById(typeIt.second);
    const MAT::PAR::ParticleMaterialBase* particlematparameter =
        dynamic_cast<const MAT::PAR::ParticleMaterialBase*>(matparameter);

    // safety check
    if (particlematparameter == NULL) dserror("cast to specific particle material failed!");

    // relate particle types to particle material parameters
    phasetypetoparticlematpar_[type_i] = particlematparameter;
  }
}

void PARTICLEINTERACTION::MaterialHandler::Setup()
{
  // nothing to do
}

void PARTICLEINTERACTION::MaterialHandler::WriteRestart(const int step, const double time) const
{
  // nothing to do
}

void PARTICLEINTERACTION::MaterialHandler::ReadRestart(
    const std::shared_ptr<IO::DiscretizationReader> reader)
{
  // nothing to do
}
