/*---------------------------------------------------------------------------*/
/*! \file
\brief momentum handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer Sebastian Fuchs
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

namespace PARTICLEINTERACTION
{
  class SPHKernelBase;
  class MaterialHandler;
  class SPHEquationOfStateBundle;
  class SPHNeighborPairs;
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
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs);

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

    //! momentum equation (particle contribution)
    void MomentumEquationParticleContribution() const;

    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! kernel handler
    std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! equation of state bundle
    std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs_;

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
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
