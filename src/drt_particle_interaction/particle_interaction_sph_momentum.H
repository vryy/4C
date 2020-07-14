/*---------------------------------------------------------------------------*/
/*! \file
\brief momentum handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_MOMENTUM_H
#define PARTICLE_INTERACTION_SPH_MOMENTUM_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_inpar/inpar_particle.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
  class ParticleContainerBundle;
}  // namespace PARTICLEENGINE

namespace PARTICLEWALL
{
  class WallHandlerInterface;
}

namespace PARTICLEINTERACTION
{
  class SPHKernelBase;
  class MaterialHandler;
  class InteractionWriter;
  class SPHEquationOfStateBundle;
  class SPHNeighborPairs;
  class SPHVirtualWallParticle;
  class SPHMomentumFormulationBase;
  class SPHArtificialViscosity;
}  // namespace PARTICLEINTERACTION

namespace MAT
{
  namespace PAR
  {
    class ParticleMaterialSPHFluid;
  }
}  // namespace MAT

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHMomentum final
  {
   public:
    //! constructor
    explicit SPHMomentum(const Teuchos::ParameterList& params);

    //! destructor
    ~SPHMomentum();

    //! init momentum handler
    void Init();

    //! setup momentum handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::InteractionWriter> particleinteractionwriter,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs,
        const std::shared_ptr<PARTICLEINTERACTION::SPHVirtualWallParticle> virtualwallparticle);

    //! write restart of momentum handler
    void WriteRestart(const int step, const double time) const;

    //! read restart of momentum handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! insert momentum evaluation dependent states
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    //! add momentum contribution to acceleration field
    void AddAccelerationContribution() const;

   private:
    //! init momentum formulation handler
    void InitMomentumFormulationHandler();

    //! init artificial viscosity handler
    void InitArtificialViscosityHandler();

    //! setup particle interaction writer
    void SetupParticleInteractionWriter();

    //! momentum equation (particle contribution)
    void MomentumEquationParticleContribution() const;

    //! momentum equation (particle-wall contribution)
    void MomentumEquationParticleWallContribution() const;

    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! interface to particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface_;

    //! kernel handler
    std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! particle interaction writer
    std::shared_ptr<PARTICLEINTERACTION::InteractionWriter> particleinteractionwriter_;

    //! equation of state bundle
    std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs_;

    //! virtual wall particle handler
    std::shared_ptr<PARTICLEINTERACTION::SPHVirtualWallParticle> virtualwallparticle_;

    //! momentum formulation handler
    std::unique_ptr<PARTICLEINTERACTION::SPHMomentumFormulationBase> momentumformulation_;

    //! artificial viscosity handler
    std::unique_ptr<PARTICLEINTERACTION::SPHArtificialViscosity> artificialviscosity_;

    //! type of boundary particle interaction
    INPAR::PARTICLE::BoundaryParticleInteraction boundaryparticleinteraction_;

    //! type of transport velocity formulation
    INPAR::PARTICLE::TransportVelocityFormulation transportvelocityformulation_;

    //! apply transport velocity formulation
    bool applytransportvelocity_;

    //! no convection of momentum with relative velocity in transport velocity formulation
    bool norelativevelocitycontribution_;

    //! pointer to fluid material of particle types
    std::vector<const MAT::PAR::ParticleMaterialSPHFluid*> fluidmaterial_;

    //! write particle-wall interaction output
    const bool writeparticlewallinteraction_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
