/*---------------------------------------------------------------------------*/
/*! \file
\brief particle type definitions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_TYPEDEFS_H
#define PARTICLE_TYPEDEFS_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_enums.H"

#include <map>
#include <set>
#include <tuple>
#include <vector>
#include <memory>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleContainer;
  class ParticleContainerBundle;
  class ParticleObject;
}  // namespace PARTICLEENGINE

namespace DRT
{
  class Element;
}

/*---------------------------------------------------------------------------*
 | type definitions                                           sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  typedef enum PARTICLEENGINE::ParticleType TypeEnum;

  typedef enum PARTICLEENGINE::ParticleStatus StatusEnum;

  typedef enum PARTICLEENGINE::ParticleState StateEnum;

  typedef std::vector<std::vector<double>> ParticleStates;

  typedef std::shared_ptr<PARTICLEENGINE::ParticleObject> ParticleObjShrdPtr;

  typedef std::vector<std::vector<std::shared_ptr<PARTICLEENGINE::ParticleContainer>>>
      TypeStatusContainers;

  typedef std::shared_ptr<PARTICLEENGINE::ParticleContainerBundle> ParticleContainerBundleShrdPtr;

  typedef std::tuple<TypeEnum, StatusEnum, int> LocalIndexTuple;

  typedef std::shared_ptr<LocalIndexTuple> LocalIndexTupleShrdPtr;

  typedef std::vector<std::vector<std::pair<PARTICLEENGINE::ParticleType, int>>> ParticlesToBins;

  typedef std::vector<std::pair<LocalIndexTuple, LocalIndexTuple>> PotentialParticleNeighbors;

  typedef std::vector<std::pair<TypeEnum, std::vector<StateEnum>>> StatesOfTypesToRefresh;

  typedef std::vector<std::vector<int>> BinsToColWallEles;

  typedef std::vector<std::pair<LocalIndexTuple, DRT::Element*>> PotentialWallNeighbors;

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
