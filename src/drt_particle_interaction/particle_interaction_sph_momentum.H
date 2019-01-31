/*---------------------------------------------------------------------------*/
/*!
\file particle_interaction_sph_momentum.H

\brief momentum handler for smoothed particle hydrodynamics (SPH) interactions

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_MOMENTUM_H
#define PARTICLE_INTERACTION_SPH_MOMENTUM_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_inpar/inpar_particle.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 08/2018 |
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

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 08/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHMomentum
  {
   public:
    //! constructor
    explicit SPHMomentum(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHMomentum();

    //! init momentum handler
    virtual void Init();

    //! setup momentum handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::SPHEquationOfStateBundle> equationofstatebundle,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs);

    //! write restart of momentum handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of momentum handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! insert momentum evaluation dependent states
    virtual void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    //! add momentum contribution to acceleration field
    virtual void AddAccelerationContribution() const;

   protected:
    //! init momentum formulation handler
    virtual void InitMomentumFormulationHandler();

    //! init artificial viscosity handler
    virtual void InitArtificialViscosityHandler();

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

    //! artificial viscous damping factor
    const double dampingfactor_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
