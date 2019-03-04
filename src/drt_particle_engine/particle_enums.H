/*---------------------------------------------------------------------------*/
/*!
\file particle_enums.H

\brief defining enums of particle states

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_ENUMS_H_
#define PARTICLE_ENUMS_H_

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_lib/drt_dserror.H"

namespace PARTICLEENGINE
{
  //! enums of particle states
  enum ParticleState
  {
    Radius,
    Mass,
    Density,
    DensitySum,
    DensityDot,
    Pressure,
    Temperature,
    TemperatureDot,
    Position,
    Velocity,
    Acceleration,
    Force,
    LastTransferPosition,
    ReferencePosition,
    ModifiedVelocity,
    ModifiedAcceleration,
    BoundaryPressure,
    BoundaryVelocity,
    Colorfield,
    ColorfieldGradient,
    InterfaceNormal,
    UnitWallNormal,
    WallDistance,
    Curvature,
    TemperatureGradient
  };

  //! get dimension of particle state
  int EnumToStateDim(const enum PARTICLEENGINE::ParticleState& stateEnum);

  //! get name of particle state
  std::string EnumToStateName(const enum PARTICLEENGINE::ParticleState& stateEnum);

  //! get enum of particle states
  enum PARTICLEENGINE::ParticleState EnumFromStateName(const std::string& stateName);


  //! enums of particle types
  enum ParticleType
  {
    Phase1,
    Phase2,
    BoundaryPhase,
    RigidPhase
  };

  //! get name of particle types
  std::string EnumToTypeName(const enum PARTICLEENGINE::ParticleType& typeEnum);

  //! get enum of particle types
  enum PARTICLEENGINE::ParticleType EnumFromTypeName(const std::string& typeName);


  //! enums of particle status
  enum ParticleStatus
  {
    Owned,
    Ghosted
  };

  //! get name of particle status
  std::string EnumToStatusName(const enum PARTICLEENGINE::ParticleStatus& statusEnum);

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
