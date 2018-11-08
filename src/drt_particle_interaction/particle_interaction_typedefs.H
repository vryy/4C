/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_typedefs.H

\brief particle interaction type definitions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_TYPEDEFS_H
#define PARTICLE_INTERACTION_TYPEDEFS_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include <map>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  struct ParticlePairSPH;
  struct ParticlePairDEM;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | type definitions                                           sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  typedef std::map<PARTICLEENGINE::TypeEnum,
      std::map<int, std::map<PARTICLEENGINE::TypeEnum,
                        std::map<PARTICLEENGINE::StatusEnum,
                            std::map<int, PARTICLEINTERACTION::ParticlePairSPH>>>>>
      SPHNeighborPairsMap;

  typedef std::map<PARTICLEENGINE::TypeEnum,
      std::map<int, std::map<PARTICLEENGINE::TypeEnum,
                        std::map<PARTICLEENGINE::StatusEnum,
                            std::map<int, PARTICLEINTERACTION::ParticlePairDEM>>>>>
      DEMNeighborPairsMap;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
