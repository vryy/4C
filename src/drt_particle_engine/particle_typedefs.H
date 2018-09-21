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

  typedef std::map<TypeEnum, std::map<StatusEnum, ParticleContainerShrdPtr>> TypeStatusContainerMap;

  typedef std::map<TypeEnum, std::map<StatusEnum, std::set<int>>> TypeStatusIndexMap;

  typedef std::map<TypeEnum, std::map<int, TypeStatusIndexMap>> ParticleNeighborsMap;

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
