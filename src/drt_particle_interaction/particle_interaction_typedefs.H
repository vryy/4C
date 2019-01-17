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

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  struct SPHParticlePair;
  struct DEMParticlePair;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | type definitions                                           sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  //! smoothed particle hydrodynamics (SPH) specific typedefs

  typedef std::pair<PARTICLEENGINE::LocalIndexTuple, PARTICLEINTERACTION::SPHParticlePair>
      SPHNeighborPair;

  typedef std::vector<std::vector<std::vector<SPHNeighborPair>>> SPHNeighborPairData;


  //! discrete element method (DEM) specific typedefs

  typedef std::pair<PARTICLEENGINE::LocalIndexTuple, PARTICLEINTERACTION::DEMParticlePair>
      DEMNeighborPair;

  typedef std::vector<std::vector<std::vector<DEMNeighborPair>>> DEMNeighborPairData;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
