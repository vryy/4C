/*----------------------------------------------------------------------*/
/*! \file
\brief Input parameters for particle problems (old implementation)

\level 3

\maintainer  Sebastian Fuchs

*-----------------------------------------------------------------------*/
#ifndef INPAR_PARTICLE_OLD_H_
#define INPAR_PARTICLE_OLD_H_

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT
{
  namespace INPUT
  {
    class ConditionDefinition;
  }
}  // namespace DRT

/*----------------------------------------------------------------------*/

namespace INPAR
{
  namespace PARTICLEOLD
  {
    /// Type of time integrator
    enum DynamicType
    {
      dyna_expleuler,  ///< forward Euler (explicit)
      dyna_centrdiff,  ///< central differences (explicit)
      dyna_rk2,        ///< Runge-Kutta 2 (explicit)
      dyna_rk4         ///< Runge-Kutta 4 (explicit)
    };

    /// Type of repartitioning strategy
    enum RepartitionStrategy
    {
      repstr_adaptive,  ///< only do repartitioning in case physically necessary
      repstr_everydt    ///< do repartitioning every time step
    };

    /// Contact type for particles
    enum ParticleInteractions
    {
      None,
      Normal_DEM,
      Normal_MD,
      NormalAndTang_DEM
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

    /// Rolling contact law
    enum RollingContact
    {
      rolling_none,
      rolling_viscous,
      rolling_constant
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

    /// random distribution of surface energy
    enum AdhesionSurfaceEnergyDistribution
    {
      adhesionsurfaceenergydistribution_none,
      adhesionsurfaceenergydistribution_lognormal,
      adhesionsurfaceenergydistribution_normal
    };

    enum ExtendedGhosting
    {
      StandardGhosting,
      WallElementGhosting,
      AddLayerGhosting
    };

    /// set the particle parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set the particle conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition>>& condlist);

  }  // namespace PARTICLEOLD
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif /* INPAR_PARTICLE_OLD_H_ */
