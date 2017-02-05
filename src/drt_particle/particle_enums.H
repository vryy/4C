/*---------------------------------------------------------------------*/
/*!
\file particle_enums.H

\brief enumerators for particle dynamics

\level 3

\maintainer Alessandro Cattabiani
*----------------------------------------------------------------------*/

#ifndef SRC_DRT_PARTICLE_PARTICLE_ENUMS_H_
#define SRC_DRT_PARTICLE_PARTICLE_ENUMS_H_


namespace PARTICLE
{
  enum StateVectorType
  {
    Dis,
    Vel,
    ColorFieldGradient,
    Radius,
    Density,
    DensityDot,
    SpecEnthalpy,
    Mass,
    Pressure,
    Temperature,
    Alpha,
//    densityApprox,
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
