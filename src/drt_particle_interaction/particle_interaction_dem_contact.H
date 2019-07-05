/*---------------------------------------------------------------------------*/
/*!
\brief contact handler for discrete element method (DEM) interactions

\level 3

\maintainer  Sebastian Fuchs
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_DEM_CONTACT_H
#define PARTICLE_INTERACTION_DEM_CONTACT_H

/*---------------------------------------------------------------------------*
 | headers                                                    sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_inpar/inpar_particle.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                       sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
namespace IO
{
  class DiscretizationReader;
}

namespace PARTICLEALGORITHM
{
  class WallHandlerInterface;
}

namespace PARTICLEENGINE
{
  class ParticleEngineInterface;
  class ParticleContainerBundle;
}  // namespace PARTICLEENGINE

namespace PARTICLEINTERACTION
{
  class MaterialHandler;
  class DEMNeighborPairs;
  class DEMHistoryPairs;
  class DEMContactNormalBase;
  class DEMContactTangentialBase;
  class DEMContactRollingBase;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class definitions                                          sfuchs 11/2018 |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMContact
  {
   public:
    //! constructor
    explicit DEMContact(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~DEMContact();

    //! init contact handler
    virtual void Init();

    //! setup contact handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEALGORITHM::WallHandlerInterface> particlewallinterface,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::DEMNeighborPairs> neighborpairs,
        const std::shared_ptr<PARTICLEINTERACTION::DEMHistoryPairs> historypairs);

    //! write restart of contact handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of contact handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! set current step size
    virtual void SetCurrentStepSize(const double currentstepsize) final;

    //! insert contact evaluation dependent states
    virtual void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    //! have tangential contact
    virtual bool HaveTangentialContact() const { return contacttangential_ ? true : false; };

    //! have rolling contact
    virtual bool HaveRollingContact() const { return contactrolling_ ? true : false; };

    //! check critical time step (on this processor)
    virtual void CheckCriticalTimeStep() const;

    //! add contact contribution to force and moment field
    virtual void AddForceAndMomentContribution();

   private:
    //! init normal contact handler
    virtual void InitNormalContactHandler();

    //! init tangential contact handler
    virtual void InitTangentialContactHandler();

    //! init rolling contact handler
    virtual void InitRollingContactHandler();

    //! get maximum density of all materials
    virtual double GetMaxDensityOfAllMaterials() const final;

    //! evaluate particle contact contribution
    virtual void EvaluateParticleContact();

    //! evaluate particle-wall contact contribution
    virtual void EvaluateParticleWallContact();

    //! discrete element method specific parameter list
    const Teuchos::ParameterList& params_dem_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! interface to particle wall handler
    std::shared_ptr<PARTICLEALGORITHM::WallHandlerInterface> particlewallinterface_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::DEMNeighborPairs> neighborpairs_;

    //! history pair handler
    std::shared_ptr<PARTICLEINTERACTION::DEMHistoryPairs> historypairs_;

    //! normal contact handler
    std::unique_ptr<PARTICLEINTERACTION::DEMContactNormalBase> contactnormal_;

    //! tangential contact handler
    std::unique_ptr<PARTICLEINTERACTION::DEMContactTangentialBase> contacttangential_;

    //! rolling contact handler
    std::unique_ptr<PARTICLEINTERACTION::DEMContactRollingBase> contactrolling_;

    //! time step size
    double dt_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
