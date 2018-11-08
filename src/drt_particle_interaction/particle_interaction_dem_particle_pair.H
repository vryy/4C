/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_dem_particle_pair.H

\brief particle pair structs for discrete element method (DEM) interactions

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
#ifndef PARTICLE_INTERACTION_DEM_PARTICLE_PAIR_H
#define PARTICLE_INTERACTION_DEM_PARTICLE_PAIR_H

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  //! struct to store quantities of interacting particles
  struct ParticlePairDEM
  {
    //! penetration between particles
    double gap_;

    //! versor from particle i to j
    double e_ji_[3];

    //! effective mass of particles i and j
    double m_eff_;
  };
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
