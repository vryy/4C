/*---------------------------------------------------------------------------*/
/*!
\file particle_wall_interface.H

\brief pure virtual interface to provide access to particle wall handler

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_WALL_INTERFACE_H
#define PARTICLE_WALL_INTERFACE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_typedefs.H"

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEALGORITHM
{
  class WallHandlerInterface
  {
   public:
    //! virtual destructor
    virtual ~WallHandlerInterface(){};

    //! get reference to potential wall neighbors
    virtual const PARTICLEENGINE::PotentialWallNeighbors& GetPotentialWallNeighbors() const = 0;
  };

}  // namespace PARTICLEALGORITHM

/*---------------------------------------------------------------------------*/
#endif
