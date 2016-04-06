/*----------------------------------------------------------------------*/
/*!
\file inpar_cavitation.H

\brief Input parameters for particle and cavitation problems

\maintainer Georg Hammerl
*-----------------------------------------------------------------------*/
#ifndef INPAR_CAVITATION_H
#define INPAR_CAVITATION_H

#include "inpar_parameterlist_utils.H"

// forward declaration
namespace DRT { namespace INPUT { class ConditionDefinition; } }

/*----------------------------------------------------------------------*/

namespace INPAR
{
  /// Input parameters for particle/cavitation problems
  namespace PARTICLE
  {

    /// Type of time integrator
    enum DynamicType
    {
      dyna_expleuler,                  ///< forward Euler (explicit)
      dyna_centrdiff,                  ///< central differences (explicit)
      dyna_rk2,                        ///< Runge-Kutta 2 (explicit)
      dyna_rk4                         ///< Runge-Kutta 4 (explicit)
    };

    /// Contact type for particles
    enum ContactStrategy
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

    /// set the particle parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set the particle conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

  }

  /// Input parameters for cavitation
  namespace CAVITATION
  {

    /// Type of coupling strategy for cavitation problems
    enum CouplingStrategyOverFields
    {
      OneWay,
      TwoWayMomentum,
      VoidFracOnly,
      TwoWayFull_weak,
      TwoWayFull_strong
    };

    /// Strategy for calculation of void fraction of fluid elements
    enum VoidFractionCalculation
    {
      analytical_constpoly,
      analytical_quadraticpoly,
      gaussian_integration
    };

    //! reconstruction type of fluid fraction
    enum FluidFracReconstructionMethod
    {
      fluidfracreco_spr,
      fluidfracreco_l2
    };

    /// set the cavitation parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

  }  // namespace CAVITATION

}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif  // #ifndef INPAR_CAVITATION_H
