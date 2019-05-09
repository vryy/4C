/*---------------------------------------------------------------------------*/
/*!
\brief input parameters for particle problems

\level 3

\maintainer  Sebastian Fuchs

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#ifndef INPAR_PARTICLE_H
#define INPAR_PARTICLE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
#include "inpar_parameterlist_utils.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

/*---------------------------------------------------------------------------*
 | input parameters for particle problems                     sfuchs 03/2018 |
 *---------------------------------------------------------------------------*/
namespace INPAR
{
  namespace PARTICLE
  {
    /*---------------------------------------------------------------------------*
     | general control parameters for particle simulations                       |
     *---------------------------------------------------------------------------*/

    //! type of particle time integration
    enum DynamicType
    {
      dyna_semiimpliciteuler,  //! semi implicit euler scheme (explicit)
      dyna_velocityverlet      //! velocity verlet scheme (explicit)
    };

    //! type of particle interaction
    enum InteractionType
    {
      interaction_none,  //! no particle interaction
      interaction_sph,   //! smoothed particle hydrodynamics
      interaction_dem    //! discrete element method
    };

    //! data format for written numeric data via vtp
    enum OutputDataFormat
    {
      binary,
      ascii
    };

    //! type of particle wall source
    enum ParticleWallSource
    {
      NoParticleWall,    //! no particle wall
      DiscretCondition,  //! particle wall from discretization condition
      BoundingBox        //! particle wall from bounding box
    };

    /*---------------------------------------------------------------------------*
     | smoothed particle hydrodynamics (SPH) specific control parameters         |
     *---------------------------------------------------------------------------*/

    //! type of smoothed particle hydrodynamics kernel
    enum KernelType
    {
      CubicSpline,
      QuinticSpline
    };

    //! kernel space dimension number
    enum KernelSpaceDimension
    {
      Kernel1D,
      Kernel2D,
      Kernel3D
    };

    //! type of smoothed particle hydrodynamics equation of state
    enum EquationOfStateType
    {
      GenTait,
      IdealGas
    };

    //! type of smoothed particle hydrodynamics momentum formulation
    enum MomentumFormulationType
    {
      AdamiMomentumFormulation,
      MonaghanMomentumFormulation
    };

    //! type of density evaluation scheme
    enum DensityEvaluationScheme
    {
      DensitySummation,
      DensityIntegration,
      DensityPredictCorrect
    };

    //! type of density correction scheme
    enum DensityCorrectionScheme
    {
      NoCorrection,
      InteriorCorrection,
      NormalizedCorrection,
      RandlesCorrection
    };

    //! type of boundary particle formulation
    enum BoundaryParticleFormulationType
    {
      NoBoundaryFormulation,
      AdamiBoundaryFormulation
    };

    //! type of boundary particle interaction
    enum BoundaryParticleInteraction
    {
      NoSlipBoundaryParticle,
      FreeSlipBoundaryParticle
    };

    //! type of transport velocity formulation
    enum TransportVelocityFormulation
    {
      NoTransportVelocity,
      StandardTransportVelocity,
      GeneralizedTransportVelocity
    };

    //! type of temperature evaluation scheme
    enum TemperatureEvaluationScheme
    {
      NoTemperatureEvaluation,
      TemperatureIntegration
    };

    //! type of heat source
    enum HeatSourceType
    {
      NoHeatSource,
      VolumeHeatSource,
      SurfaceHeatSource
    };

    //! type of surface tension formulation
    enum SurfaceTensionFormulation
    {
      NoSurfaceTension,
      ContinuumSurfaceForce
    };

    //! type of phase change
    enum PhaseChangeType
    {
      NoPhaseChange,
      TwoWayScalarPhaseChange
    };

    /*---------------------------------------------------------------------------*
     | discrete element method (DEM) specific control parameters                 |
     *---------------------------------------------------------------------------*/

    //! type of normal contact law
    enum NormalContact
    {
      NormalLinSpring,
      NormalLinSpringDamp,
      NormalHertz,
      NormalLeeHerrmann,
      NormalKuwabaraKono,
      NormalTsuji
    };

    //! type of tangential contact law
    enum TangentialContact
    {
      NoTangentialContact,
      TangentialLinSpringDamp
    };


    //! set the particle parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace PARTICLE

}  // namespace INPAR

/*---------------------------------------------------------------------------*/
#endif
