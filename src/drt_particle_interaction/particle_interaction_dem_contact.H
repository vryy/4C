/*---------------------------------------------------------------------------*/
/*! \file
\brief contact handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_DEM_CONTACT_H
#define PARTICLE_INTERACTION_DEM_CONTACT_H

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
  class MaterialHandler;
  class InteractionWriter;
  class DEMNeighborPairs;
  class DEMHistoryPairs;
  class DEMContactNormalBase;
  class DEMContactTangentialBase;
  class DEMContactRollingBase;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMContact final
  {
   public:
    //! constructor
    explicit DEMContact(const Teuchos::ParameterList& params);

    //! destructor
    ~DEMContact();

    //! init contact handler
    void Init();

    //! setup contact handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
        const std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial,
        const std::shared_ptr<PARTICLEINTERACTION::InteractionWriter> particleinteractionwriter,
        const std::shared_ptr<PARTICLEINTERACTION::DEMNeighborPairs> neighborpairs,
        const std::shared_ptr<PARTICLEINTERACTION::DEMHistoryPairs> historypairs);

    //! write restart of contact handler
    void WriteRestart(const int step, const double time) const;

    //! read restart of contact handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! set current step size
    void SetCurrentStepSize(const double currentstepsize);

    //! insert contact evaluation dependent states
    void InsertParticleStatesOfParticleTypes(
        std::map<PARTICLEENGINE::TypeEnum, std::set<PARTICLEENGINE::StateEnum>>&
            particlestatestotypes) const;

    //! get normal contact stiffness
    double GetNormalContactStiffness() const;

    //! have tangential contact
    bool HaveTangentialContact() const { return contacttangential_ ? true : false; };

    //! have rolling contact
    bool HaveRollingContact() const { return contactrolling_ ? true : false; };

    //! check critical time step (on this processor)
    void CheckCriticalTimeStep() const;

    //! add contact contribution to force and moment field
    void AddForceAndMomentContribution();

   private:
    //! init normal contact handler
    void InitNormalContactHandler();

    //! init tangential contact handler
    void InitTangentialContactHandler();

    //! init rolling contact handler
    void InitRollingContactHandler();

    //! setup particle interaction writer
    void SetupParticleInteractionWriter();

    //! get maximum density of all materials
    double GetMaxDensityOfAllMaterials() const;

    //! evaluate particle contact contribution
    void EvaluateParticleContact();

    //! evaluate particle-wall contact contribution
    void EvaluateParticleWallContact();

    //! discrete element method specific parameter list
    const Teuchos::ParameterList& params_dem_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! interface to particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface_;

    //! particle material handler
    std::shared_ptr<PARTICLEINTERACTION::MaterialHandler> particlematerial_;

    //! particle interaction writer
    std::shared_ptr<PARTICLEINTERACTION::InteractionWriter> particleinteractionwriter_;

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

    //! write particle-wall interaction output
    const bool writeparticlewallinteraction_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
