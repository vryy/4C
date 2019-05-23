/*---------------------------------------------------------------------------*/
/*!

\brief base particle interaction handler

\level 3

\maintainer  Sebastian Fuchs
             fuchs@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289 -15262

*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_BASE_H
#define PARTICLE_INTERACTION_BASE_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "../linalg/linalg_fixedsizematrix.H"

#include <Epetra_Comm.h>
#include <Teuchos_ParameterList.hpp>

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

namespace PARTICLEINTERACTION
{
  class MaterialHandler;
}

namespace PARTICLEALGORITHM
{
  class WallHandlerInterface;
}

namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
}

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 05/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class ParticleInteractionBase
  {
   public:
    //! constructor
    explicit ParticleInteractionBase(const Epetra_Comm& comm, const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~ParticleInteractionBase(){};

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
            particlestatestotypes) = 0;

    //! set initial states
    virtual void SetInitialStates() = 0;

    //! evaluate particle interactions
    virtual void EvaluateInteractions() = 0;

    //! maximum interaction distance (on this processor)
    virtual double MaxInteractionDistance() const = 0;

    //! distribute interaction history
    virtual void DistributeInteractionHistory() const = 0;

    //! communicate interaction history
    virtual void CommunicateInteractionHistory() const = 0;

    //! set current time
    virtual void SetCurrentTime(const double currenttime);

    //! set current step size
    virtual void SetCurrentStepSize(const double currentstepsize);

    //! set gravity
    virtual void SetGravity(std::vector<double>& gravity) final;

   private:
    //! init particle material handler
    virtual void InitParticleMaterialHandler();

   protected:
    //! maximum particle radius (on this processor)
    virtual double MaxParticleRadius() const;

    //! communication
    const Epetra_Comm& comm_;

    //! processor id
    const int myrank_;

    //! particle simulation parameter list
    const Teuchos::ParameterList& params_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! interface to particle wall handler
    std::shared_ptr<PARTICLEALGORITHM::WallHandlerInterface> particlewallinterface_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! current time
    double time_;

    //! time step size
    double dt_;

    //! current gravity vector
    std::vector<double> gravity_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
