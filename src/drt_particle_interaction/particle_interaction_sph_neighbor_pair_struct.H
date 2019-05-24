/*---------------------------------------------------------------------------*/
/*!
\brief neighbor pair struct for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_NEIGHBOR_PAIR_STRUCT_H
#define PARTICLE_INTERACTION_SPH_NEIGHBOR_PAIR_STRUCT_H

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  //! struct to store quantities of interacting particles
  struct SPHNeighborPair
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
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
