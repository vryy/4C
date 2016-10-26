/*----------------------------------------------------------------------*/
/*!
\file inpar_particle.H

\brief Input parameters for particle problems

\level 3

\maintainer Georg Hammerl
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
      dyna_rk4                         ///< Runge-Kutta 4 (explicit)
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

    /// Weight function for MeshFree interactions
    enum WeightFunction
    {
      CubicBspline
    };

    /// wall interaction type for MeshFree interactions
    enum WallInteractionType
    {
      InitParticle,
      Mirror,
      Custom
    };

    /// set the particle parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set the particle conditions
    void SetValidConditions(std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist);

  }
}  // namespace INPAR

/*----------------------------------------------------------------------*/
#endif /* INPAR_PARTICLE_H_ */
