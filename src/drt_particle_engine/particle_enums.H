/*---------------------------------------------------------------------------*/
/*! \file
\brief defining enums for particle problem
\level 1
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_ENUMS_H
#define PARTICLE_ENUMS_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_lib/drt_dserror.H"

/*---------------------------------------------------------------------------*
 | enum definitions                                                          |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  //! \name definition of particle states
  //! @{

  /*!
   * \brief enums of particle states
   *
   * \author Sebastian Fuchs \date 03/2018
   */
  enum ParticleState
  {
    /*---------------------------------------------------------------------------*/
    // default particle states
    Position,              //!< position
    Velocity,              //!< velocity
    Acceleration,          //!< acceleration
    LastTransferPosition,  //!< position after last particle transfer
    /*---------------------------------------------------------------------------*/
    // particle states for advanced time integration schemes
    ModifiedVelocity,      //!< modified velocity
    ModifiedAcceleration,  //!< modified acceleration
    /*---------------------------------------------------------------------------*/
    // particle states for boundary conditions
    ReferencePosition,  //!< reference position
    /*---------------------------------------------------------------------------*/
    // common particle states
    Radius,       //!< radius
    Mass,         //!< mass
    Density,      //!< density
    Pressure,     //!< pressure
    Temperature,  //!< temperature
    /*---------------------------------------------------------------------------*/
    // particle states for rigid body evaluation
    RigidBodyColor,             //!< rigid body color
    RelativePosition,           //!< position relative to rigid body position in rotating frame
    ReferenceRelativePosition,  //!< position relative to rigid body position in reference frame
    /*---------------------------------------------------------------------------*/
    // particle states for SPH interaction
    DensitySum,           //!< density from summation of weighted masses
    DensityDot,           //!< first time derivative of density
    TemperatureDot,       //!< first time derivative of temperature
    BoundaryPressure,     //!< boundary pressure
    BoundaryVelocity,     //!< boundary velocity
    Colorfield,           //!< colorfield
    ColorfieldGradient,   //!< colorfield gradient
    InterfaceNormal,      //!< interface normal
    UnitWallNormal,       //!< unit wall normal
    WallDistance,         //!< wall distance
    TemperatureGradient,  //!< temperature gradient
    /*---------------------------------------------------------------------------*/
    // particle states for DEM interaction
    AngularVelocity,      //!< angular velocity
    AngularAcceleration,  //!< angular acceleration
    Force,                //!< force
    Moment,               //!< moment
    /*---------------------------------------------------------------------------*/
    // particle states for two way coupled partitioned PASI
    LastIterPosition,              //!< position after last converged iteration
    LastIterVelocity,              //!< velocity after last converged iteration
    LastIterAcceleration,          //!< acceleration after last converged iteration
    LastIterAngularVelocity,       //!< angular velocity after last converged iteration
    LastIterAngularAcceleration,   //!< angular acceleration after last converged iteration
    LastIterModifiedAcceleration,  //!< modified acceleration after last converged iteration
    LastIterDensity,               //!< density after last converged iteration
    LastIterTemperature            //!< temperature after last converged iteration
    /*---------------------------------------------------------------------------*/
  };

  /*!
   * \brief convert particle state enum to dimension
   *
   * \note This method should be used only for initialization.
   *
   * \author Sebastian Fuchs \date 03/2018
   *
   * \param[in] stateEnum particle state enum
   *
   * \return particle state dimension
   */
  int EnumToStateDim(const enum PARTICLEENGINE::ParticleState& stateEnum);

  /*!
   * \brief convert particle state enum to name
   *
   * \author Sebastian Fuchs \date 03/2018
   *
   * \param[in] stateEnum particle state enum
   *
   * \return particle state name
   */
  std::string EnumToStateName(const enum PARTICLEENGINE::ParticleState& stateEnum);

  /*!
   * \brief convert particle state name to enum
   *
   * \note Only relevant particle states, i.e., currently needed for initialization are processed.
   *       Extend any additional states if needed!
   *
   * \warning This method is computational expensive due to comparison of strings and should only be
   *          used for initialization.
   *
   * \author Sebastian Fuchs \date 03/2018
   *
   * \param[in] stateName particle state name
   *
   * \return particle state enum
   */
  enum PARTICLEENGINE::ParticleState EnumFromStateName(const std::string& stateName);

  //! @}

  //! \name definition of particle types
  //! @{

  /*!
   * \brief enums of particle type
   *
   * Enum for respective particle type used to relate and distinguish particles to different phases.
   *
   * \author Sebastian Fuchs \date 03/2018
   */
  enum ParticleType
  {
    /*---------------------------------------------------------------------------*/
    // default particle types
    Phase1,  //!< particle type for particles of first phase
    Phase2,  //!< particle type for particles of second phase
    /*---------------------------------------------------------------------------*/
    // particle types for SPH interaction
    BoundaryPhase,   //!< particle type for boundary phase particles
    RigidPhase,      //!< particle type for rigid phase particles
    DirichletPhase,  //!< particle type for dirichlet phase particles (open boundary)
    NeumannPhase     //!< particle type for neumann phase particles (open boundary)
    /*---------------------------------------------------------------------------*/
  };

  /*!
   * \brief convert particle type enum to name
   *
   * \author Sebastian Fuchs \date 03/2018
   *
   * \param[in] typeEnum particle type enum
   *
   * \return particle type name
   */
  std::string EnumToTypeName(const enum PARTICLEENGINE::ParticleType& typeEnum);

  /*!
   * \brief convert particle type name to enum
   *
   * \warning This method is computational expensive due to comparison of strings and should only be
   *          used for initialization.
   *
   * \author Sebastian Fuchs \date 03/2018
   *
   * \param[in] typeName particle type name
   *
   * \return particle type enum
   */
  enum PARTICLEENGINE::ParticleType EnumFromTypeName(const std::string& typeName);

  //! @}

  //! \name definition of particle status
  //! @{

  /*!
   * \brief enums of particle status
   *
   * Enum for respective particle status used to distinguish owned and ghosted particles on each
   * processor.
   *
   * \author Sebastian Fuchs \date 03/2018
   */
  enum ParticleStatus
  {
    Owned,   //!< particle status for particles being owned on processors
    Ghosted  //!< particle status for particles being ghosted on processors
  };

  /*!
   * \brief convert particle status enum to name
   *
   * \author Sebastian Fuchs \date 03/2018
   *
   * \param[in] statusEnum particle status enum
   *
   * \return particle status name
   */
  std::string EnumToStatusName(const enum PARTICLEENGINE::ParticleStatus& statusEnum);

  //! @}

}  // namespace PARTICLEENGINE

/*---------------------------------------------------------------------------*/
#endif
