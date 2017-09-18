/*----------------------------------------------------------------------*/
/*!
\file particle_enums.H

\brief enumerators for particle dynamics

\level 3

\maintainer  Christoph Meier
             meier@lnm.mw.tum.de
             http://www.lnm.mw.tum.de

*-----------------------------------------------------------------------*/

#ifndef SRC_DRT_PARTICLE_PARTICLE_ENUMS_H_
#define SRC_DRT_PARTICLE_PARTICLE_ENUMS_H_


namespace PARTICLE
{
  enum StateVectorType
  {
    Dis,
    Vel,
    VelConv,
    ColorFieldGradient,
    Radius,
    Density,
    DensityDot,
    SpecEnthalpy,
    Mass,
    Pressure,
    Temperature,
    mGradW,
    mHessW,
    ColorField,
    DensitySum
  };

  enum NeighbourType
  {
    NT_Particle,
    NT_Wall,
    NT_HeatSource
  };

  enum FreeSurfaceParticleType
  {
    FS_NONE,        //Particle that is not near a free surface
    FS_DIRECT,      //Particle that is directly near a free surface (low color function value)
    FS_INDIRECT     //Particle that has a particle of type FS_DIRECT as neighbor
  };
}


#endif /* SRC_DRT_PARTICLE_PARTICLE_ENUMS_H_ */
