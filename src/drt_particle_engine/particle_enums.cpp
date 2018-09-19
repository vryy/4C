/*---------------------------------------------------------------------------*/
/*!
\file particle_enums.cpp

\brief defining enums of particle states

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_enums.H"

/*---------------------------------------------------------------------------*
 | get dimension of particle state                            sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
int PARTICLEENGINE::EnumToStateDim(const enum PARTICLEENGINE::ParticleState& stateEnum)
{
  int dim = 0;

  switch (stateEnum)
  {
    // scalar states
    case PARTICLEENGINE::Radius:
    case PARTICLEENGINE::Mass:
    case PARTICLEENGINE::Density:
    case PARTICLEENGINE::Pressure:
      dim = 1;
      break;

    // vectorial states
    case PARTICLEENGINE::Position:
    case PARTICLEENGINE::Velocity:
    case PARTICLEENGINE::Acceleration:
    case PARTICLEENGINE::ReferencePosition:
    case PARTICLEENGINE::ModifiedVelocity:
    case PARTICLEENGINE::ModifiedAcceleration:
      dim = 3;
      break;

    default:
      dserror("particle state enum unknown!");
  }

  return dim;
}

/*---------------------------------------------------------------------------*
 | get name of particle state                                 sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
std::string PARTICLEENGINE::EnumToStateName(const enum PARTICLEENGINE::ParticleState& stateEnum)
{
  std::string name;

  switch (stateEnum)
  {
    case PARTICLEENGINE::Radius:
      name = "radius";
      break;
    case PARTICLEENGINE::Mass:
      name = "mass";
      break;
    case PARTICLEENGINE::Density:
      name = "density";
      break;
    case PARTICLEENGINE::Pressure:
      name = "pressure";
      break;
    case PARTICLEENGINE::Position:
      name = "position";
      break;
    case PARTICLEENGINE::Velocity:
      name = "velocity";
      break;
    case PARTICLEENGINE::Acceleration:
      name = "acceleration";
      break;
    case PARTICLEENGINE::ReferencePosition:
      name = "reference position";
      break;
    case PARTICLEENGINE::ModifiedVelocity:
      name = "modified velocity";
      break;
    case PARTICLEENGINE::ModifiedAcceleration:
      name = "modified acceleration";
      break;
    default:
      dserror("particle state enum unknown!");
  }

  return name;
}

/*---------------------------------------------------------------------------*
 | get name of particle type                                  sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
std::string PARTICLEENGINE::EnumToTypeName(const enum PARTICLEENGINE::ParticleType& typeEnum)
{
  std::string name;

  switch (typeEnum)
  {
    case PARTICLEENGINE::Phase1:
      name = "phase1";
      break;
    case PARTICLEENGINE::Phase2:
      name = "phase2";
      break;
    default:
      dserror("particle type enum unknown!");
  }

  return name;
}

/*---------------------------------------------------------------------------*
 | get enum of particle types                                 sfuchs 06/2018 |
 *---------------------------------------------------------------------------*/
enum PARTICLEENGINE::ParticleType PARTICLEENGINE::EnumFromTypeName(const std::string& typeName)
{
  // attention: this method is expensive (comparison of strings)
  //            and should be used only for initialization or result testing

  enum PARTICLEENGINE::ParticleType type;

  if (typeName == "phase1")
    type = PARTICLEENGINE::Phase1;
  else if (typeName == "phase2")
    type = PARTICLEENGINE::Phase2;
  else
    dserror("particle type '%s' unknown!", typeName.c_str());

  return type;
}

/*---------------------------------------------------------------------------*
 | get name of particle status                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
std::string PARTICLEENGINE::EnumToStatusName(const enum PARTICLEENGINE::ParticleStatus& statusEnum)
{
  std::string name;

  switch (statusEnum)
  {
    case PARTICLEENGINE::Owned:
      name = "owned";
      break;
    case PARTICLEENGINE::Ghosted:
      name = "ghosted";
      break;
    default:
      dserror("particle status enum unknown!");
  }

  return name;
}
