/*---------------------------------------------------------------------------*/
/*!
\file particle_typedefs.H

\brief particle type definitions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

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

/*---------------------------------------------------------------------------*
 | type definitions                                           sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  typedef enum PARTICLEENGINE::ParticleType TypeEnum;

  typedef enum PARTICLEENGINE::ParticleStatus StatusEnum;

  typedef enum PARTICLEENGINE::ParticleState StateEnum;

  typedef std::shared_ptr<PARTICLEENGINE::ParticleObject> ParticleObjShrdPtr;

  typedef std::shared_ptr<PARTICLEENGINE::ParticleContainer> ParticleContainerShrdPtr;

  typedef std::shared_ptr<PARTICLEENGINE::ParticleContainerBundle> ParticleContainerBundleShrdPtr;

  typedef std::map<StateEnum, std::vector<double>> ParticleStates;

  typedef std::vector<std::vector<ParticleContainerShrdPtr>> TypeStatusContainers;

  typedef std::tuple<PARTICLEENGINE::TypeEnum, PARTICLEENGINE::StatusEnum, int> LocalIndexTuple;

  typedef std::shared_ptr<LocalIndexTuple> LocalIndexTupleShrdPtr;

  typedef std::vector<std::vector<std::vector<LocalIndexTuple>>> ParticleNeighbors;
}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
