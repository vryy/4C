/*---------------------------------------------------------------------------*/
/*! \file
\brief neighbor pair struct for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef FOUR_C_PARTICLE_INTERACTION_SPH_NEIGHBOR_PAIR_STRUCT_HPP
#define FOUR_C_PARTICLE_INTERACTION_SPH_NEIGHBOR_PAIR_STRUCT_HPP

#include "4C_config.hpp"

#include "4C_particle_engine_typedefs.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::Elements
{
  class Element;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace ParticleInteraction
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
    Core::Elements::Element* ele_;

    //! absolute distance between particle and wall contact point
    double absdist_;

    //! versor from wall contact point j to particle i
    double e_ij_[3];

    //! parameter space coordinates of wall contact point
    double elecoords_[2];
  };
}  // namespace ParticleInteraction

/*---------------------------------------------------------------------------*/
FOUR_C_NAMESPACE_CLOSE

#endif
