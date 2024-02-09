/*---------------------------------------------------------------------------*/
/*! \file
\brief neighbor pair struct for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_INTERACTION_SPH_NEIGHBOR_PAIR_STRUCT_HPP
#define BACI_PARTICLE_INTERACTION_SPH_NEIGHBOR_PAIR_STRUCT_HPP

#include "baci_config.hpp"

#include "baci_particle_engine_typedefs.hpp"

BACI_NAMESPACE_OPEN

namespace DRT
{
  class Element;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  //! struct to store quantities of interacting particles
  struct SPHParticlePair final
  {
    //! local index tuple of particles i and j
    PARTICLEENGINE::LocalIndexTuple tuple_i_;
    PARTICLEENGINE::LocalIndexTuple tuple_j_;

    //! absolute distance between particles
    double absdist_;

    //! versor from particle j to i
    double e_ij_[3];

    //! kernel
    double Wij_;
    double Wji_;

    //! first derivative of kernel
    double dWdrij_;
    double dWdrji_;
  };

  //! struct to store quantities of interacting particles and wall elements
  struct SPHParticleWallPair final
  {
    //! local index tuple of particle i
    PARTICLEENGINE::LocalIndexTuple tuple_i_;

    //! pointer to column wall element
    DRT::Element* ele_;

    //! absolute distance between particle and wall contact point
    double absdist_;

    //! versor from wall contact point j to particle i
    double e_ij_[3];

    //! parameter space coordinates of wall contact point
    double elecoords_[2];
  };
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif
