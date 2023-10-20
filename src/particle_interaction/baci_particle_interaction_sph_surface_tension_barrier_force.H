/*---------------------------------------------------------------------------*/
/*! \file
\brief barrier force handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef BACI_PARTICLE_INTERACTION_SPH_SURFACE_TENSION_BARRIER_FORCE_H
#define BACI_PARTICLE_INTERACTION_SPH_SURFACE_TENSION_BARRIER_FORCE_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "baci_inpar_particle.H"
#include "baci_particle_engine_enums.H"
#include "baci_particle_engine_typedefs.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
  class ParticleContainerBundle;
}  // namespace PARTICLEENGINE

namespace PARTICLEINTERACTION
{
  class SPHNeighborPairs;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHBarrierForce
  {
   public:
    //! constructor
    explicit SPHBarrierForce(const Teuchos::ParameterList& params);

    //! init barrier force handler
    void Init();

    //! setup barrier force handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs);

    //! compute barrier force contribution
    void ComputeBarrierForceContribution() const;

   protected:
    //! compute barrier force contribution (particle contribution)
    void ComputeBarrierForceParticleContribution() const;

    //! compute barrier force contribution (particle-boundary contribution)
    void ComputeBarrierForceParticleBoundaryContribution() const;

    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs_;

    //! liquid particle type
    PARTICLEENGINE::TypeEnum liquidtype_;

    //! gas particle type
    PARTICLEENGINE::TypeEnum gastype_;

    //! set of fluid particle types
    std::set<PARTICLEENGINE::TypeEnum> fluidtypes_;

    //! set of boundary particle types
    std::set<PARTICLEENGINE::TypeEnum> boundarytypes_;

    //! barrier force distance
    const double dist_;

    //! barrier force temperature scaling
    const double cr_;

    //! transition reference temperature
    const double trans_ref_temp_;

    //! transition temperature difference for barrier force evaluation
    const double trans_dT_barrier_;

    //! barrier force stiffness of heavy phase
    const double stiff_h_;

    //! barrier force damping parameter of heavy phase
    const double damp_h_;

    //! barrier force stiffness of gas phase
    const double stiff_g_;

    //! barrier force damping parameter of gas phase
    const double damp_g_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
