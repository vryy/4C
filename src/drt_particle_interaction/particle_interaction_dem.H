/*---------------------------------------------------------------------------*/
/*!
\brief discrete element method (DEM) interaction handler

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_DEM_H
#define PARTICLE_INTERACTION_DEM_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include "particle_interaction_base.H"

#include "../drt_inpar/inpar_particle.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMNeighborPairs;
  class DEMHistoryPairs;
  class DEMContact;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class ParticleInteractionDEM : public ParticleInteractionBase
  {
   public:
    //! constructor
    explicit ParticleInteractionDEM(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~ParticleInteractionDEM();

    //! init particle interaction handler
    virtual void Init();

    //! setup particle interaction handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEALGORITHM::WallHandlerInterface> particlewallinterface);

    //! write restart of particle interaction handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of particle interaction handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! insert interaction dependent states of all particle types
    virtual void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) override;

    //! set initial states
    virtual void SetInitialStates() override;

    //! evaluate particle interactions
    virtual void EvaluateInteractions() override;

    //! maximum interaction distance (on this processor)
    virtual double MaxInteractionDistance() const override;

    //! distribute interaction history
    virtual void DistributeInteractionHistory() const override;

    //! communicate interaction history
    virtual void CommunicateInteractionHistory() const override;

    //! set current step size
    virtual void SetCurrentStepSize(const double currentstepsize);

   private:
    //! init neighbor pair handler
    virtual void InitNeighborPairHandler();

    //! init history pair handler
    virtual void InitHistoryPairHandler();

    //! init contact handler
    virtual void InitContactHandler();

    //! set initial radius
    virtual void SetInitialRadius();

    //! set initial mass
    virtual void SetInitialMass();

    //! clear force and moment states of particles
    virtual void ClearForceAndMomentStates() const;

    //! compute acceleration from force and moment
    virtual void ComputeAcceleration() const;

    //! discrete element method specific parameter list
    const Teuchos::ParameterList& params_dem_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::DEMNeighborPairs> neighborpairs_;

    //! history pair handler
    std::shared_ptr<PARTICLEINTERACTION::DEMHistoryPairs> historypairs_;

    //! contact handler
    std::unique_ptr<PARTICLEINTERACTION::DEMContact> contact_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
