/*----------------------------------------------------------------------*/
/*!
\file inpar_particle.H

\brief Input parameters for particle problems

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*-----------------------------------------------------------------------*/
#ifndef INPAR_PARTICLE_H_
#define INPAR_PARTICLE_H_

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

/*----------------------------------------------------------------------*/

namespace INPAR
{
  /// Input parameters for particle problems
  namespace PARTICLE
  {

    /// Type of time integrator
    enum DynamicType
    {
      dyna_expleuler,                  ///< forward Euler (explicit)
      dyna_centrdiff,                  ///< central differences (explicit)
      dyna_kickdrift,                  ///< kick-dirft scheme (explicit)
      dyna_rk2,                        ///< Runge-Kutta 2 (explicit)
      dyna_rk4                         ///< Runge-Kutta 4 (explicit)
    };

    /// Contact type for particles
    enum ParticleInteractions
    {
      None,
      Normal_DEM,
      Normal_MD,
      NormalAndTang_DEM,
      SPH
    };

    /// Normal contact law
    enum NormalContact
    {
      LinSpring,
      Hertz,
      LinSpringDamp,
      LeeHerrmann,
      KuwabaraKono,
      Tsuji
    };

    /// adhesion law governing normal contact of particles
    enum NormalAdhesion
    {
      adhesion_none,
      adhesion_linspring,
      adhesion_linspringdamp,
      adhesion_vdWDMT,
      adhesion_regDMT
    };

    //! type of collision
    enum coltype
    {
      particle_particle,
      particle_wall,
      undefined
    };

    /// problem dimension in case of quasi 2D
    enum ParticleDim
    {
      particle_3D,
      particle_2Dx,
      particle_2Dy,
      particle_2Dz
    };

    /// random distribution of particle radii
    enum RadiusDistribution
    {
      radiusdistribution_none,
      radiusdistribution_lognormal,
      radiusdistribution_normal
    };

    /// Weight function for SPH interactions
    enum WeightFunction
    {
      CubicBspline,
      QuinticBspline,
      SqrtHyperbola,
      HyperbolaNoRsz
    };

    /// Number of weight function space dimensions for SPH interactions
    enum WeightFunctionDim
    {
      WF_1D,
      WF_2D,
      WF_3D
    };

    /// Equation of state for SPH interactions
    enum EquationOfState
    {
      GenTait,
      IdealGas
    };

    enum ExtendedGhosting
    {
      StandardGhosting,
      BdryParticleGhosting,
      WallElementGhosting,
      AddLayerGhosting
    };

    /// wall interaction type for SPH interactions
    enum WallInteractionType
    {
      BoundarParticle_NoSlip,
      BoundarParticle_FreeSlip,
      NoWallInteraction
    };

    /// type of free-surface treatment for SPH interactions
    enum FreeSurfaceType
    {
      FreeSurface_None,               //no free surfaces
      DensityIntegration,             //density determined via time integration exclusively
      InteriorReinitialization,       //density via time integration and re-initialization in the interior via density summation (in every time step)
      NormalizedReinitialization,     //like 'InteriorReinitialization' but with additional re-initialization of free-surface particles via normalized density summation (Shepard Filter)
      RandlesReinitialization,        //like 'NormalizedReinitialization' but with re-initialization of free-surface particles inspired by Randles 1996 instead of Shepard Filter
      TwoPhase                        //void phase explicitly modeled / discretized
    };

    /// type of surface tension forces for SPH interactions
    enum SurfaceTensionType
    {
      ST_NONE,             //no surface tension
      ST_VDW_DIRECT,       //surface tension via van der Waals forces according to Tartakovsky et al. 2016, direct definition of surface tension and contact angle
      ST_VDW_INDIRECT,     //surface tension via van der Waals forces according to Tartakovsky et al. 2016, definition of surface tension and contact angle via s_ff, s_sf
      ST_CONTI_ADAMI      //surface tension via continuum approach (pressure jump ~ interface curvature) according to Adami et al. 2010
    };

    /// rendering for smoothed particle hydrodynamics
    enum RenderingType
    {
      NoRendering,
      StandardRendering,
      NormalizedRendering
    };

    enum RenderingOutput
    {
      DiscretAndMatlab,
      Discret,
      Matlab
    };

    enum RenderingBdryPart
    {
      WithBdryParticle,
      NoBdryParticle
    };

    /// set the particle parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set the particle conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

  }
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif /* INPAR_PARTICLE_H_ */
