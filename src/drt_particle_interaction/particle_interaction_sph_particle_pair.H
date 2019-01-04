/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_particle_pair.H

\brief particle pair structs for smoothed particle hydrodynamics (SPH) interactions

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
#ifndef PARTICLE_INTERACTION_SPH_PARTICLE_PAIR_H
#define PARTICLE_INTERACTION_SPH_PARTICLE_PAIR_H

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  //! struct to store quantities of interacting particles
  struct SPHParticlePair
  {
    //! absolute distance between particles
    double absdist_;

    //! versor from particle j to i
    double e_ij_[3];

    //! kernel
    double Wij_;

    //! first derivative of kernel
    double dWdrij_;
  };
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
