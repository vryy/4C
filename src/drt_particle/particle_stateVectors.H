/*---------------------------------------------------------------------*/
/*!
\file particle_stateVectors.H

\brief state vector enumerators

\level 3

\maintainer Alessandro Cattabiani
*----------------------------------------------------------------------*/

#ifndef SRC_DRT_PARTICLE_PARTICLE_STATEVECTORS_H_
#define SRC_DRT_PARTICLE_PARTICLE_STATEVECTORS_H_


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
}


#endif /* SRC_DRT_PARTICLE_PARTICLE_STATEVECTORS_H_ */
