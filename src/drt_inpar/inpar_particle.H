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
      dyna_rk2,                        ///< Runge-Kutta 2 (explicit)
      dyna_rk4,                        ///< Runge-Kutta 4 (explicit)
      dyna_hybridMeshFreeDivFree,      ///< divergence-free corrections (hybrid implicit-explicit)
      dyna_genAlpha                    ///< general alpha scheme (implicit)
    };

    /// Contact type for particles
    enum ParticleInteractions
    {
      None,
      Normal_DEM,
      Normal_MD,
      NormalAndTang_DEM,
      Normal_DEM_thermo,
      MeshFree
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
      adhesion_lennardjones,
      adhesion_linspring,
      adhesion_linspringdamp
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

    /// Weight function for MeshFree interactions
    enum WeightFunction
    {
      CubicBspline,
      SqrtHyperbola,
      HyperbolaNoRsz
    };

    /// Number of weight function space dimensions for MeshFree interactions
    enum WeightFunctionDim
    {
      WF_3D,
      WF_2D,
      WF_1D
    };

    /// wall interaction type for MeshFree interactions
    enum WallInteractionType
    {
      InitParticle,
      Mirror,
      Custom,
      BoundarParticle_NoSlip,
      BoundarParticle_FreeSlip,
      NoWallInteraction
    };

    /// time step choice
    enum TimeStepType
    {
      Manual,      /// dt is the prescribed time step
      Auto_CFL     /// Replace the prescribed dt with the  CFL condition if it is more strict. Courant-Friedrich-Levy condition: dt = 0.49 * radius_min / str::abs(vel_max)
    };

    enum MidAverageEnum
    {
      midavg_vague,
      midavg_imrlike,  ///< alphaf-mid-averaging is done IMR-like, i.e.
                       ///< \f$F_{int,m}\f$
                       ///< \f$= F_{int}(D_m)\f$
                       ///< \f$= F_{int}( (1-\alpha_f)*D_{n+1} + \alpha_f*D_n )\f$
                       ///< (IMR means implicit mid-point rule.)
      midavg_trlike    ///< alphaf-mid-averaging is done TR-like, i.e.
                       ///< \f$F_{int,m}\f$
                       ///< \f$= (1-\alpha_f)*F_{int,n+1} + \alpha_f*F_{int,n}\f$
                       ///<  \f$= (1-\alpha_f)*F_{int}(D_{n+1}) + \alpha_f*F_{int}(D_n)\f$
                       ///<  (TR means trapezoidal rule.)
    };

    /// Map mid-averaging to std::string
    static inline const std::string MidAverageString
    (
      const enum MidAverageEnum name  ///< enum to convert
    )
    {
      switch (name)
      {
        break;
      case midavg_imrlike :
        return "IMR-like";
      case midavg_trlike :
        return "TR-like";
        break;
      default :
        dserror("Cannot make std::string for time integrator %d", name);
        return "";
      }
    }

    /// approximations of the generalized-alpha scheme
    enum GAapproxType
    {
      gaapprox_full,
      gaapprox_onlyhess
    };

    /// set the particle parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set the particle conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

  }
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif /* INPAR_PARTICLE_H_ */
