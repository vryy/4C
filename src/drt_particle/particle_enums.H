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
    mHessW
  };

  enum NeighbourType
  {
    NT_Particle,
    NT_Wall,
    NT_HeatSource
  };
}


#endif /* SRC_DRT_PARTICLE_PARTICLE_ENUMS_H_ */
