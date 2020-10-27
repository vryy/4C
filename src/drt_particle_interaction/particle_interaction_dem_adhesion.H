/*---------------------------------------------------------------------------*/
/*! \file
\brief adhesion handler for discrete element method (DEM) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_DEM_ADHESION_H
#define PARTICLE_INTERACTION_DEM_ADHESION_H

/*---------------------------------------------------------------------------*
 | headers                                                                   |
 *---------------------------------------------------------------------------*/
#include "../drt_particle_engine/particle_enums.H"
#include "../drt_particle_engine/particle_typedefs.H"

#include "../drt_inpar/inpar_particle.H"

/*---------------------------------------------------------------------------*
 | forward declarations                                                      |
 *---------------------------------------------------------------------------*/
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
  class InteractionWriter;
  class DEMNeighborPairs;
  class DEMHistoryPairs;
  class DEMAdhesionLawBase;
  class DEMAdhesionSurfaceEnergyBase;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class DEMAdhesion final
  {
   public:
    //! constructor
    explicit DEMAdhesion(const Teuchos::ParameterList& params);

    /*!
     * \brief destructor
     *
     * \author Sebastian Fuchs \date 07/2019
     *
     * \note At compile-time a complete type of class T as used in class member
     *       std::unique_ptr<T> ptr_T_ is required
     */
    ~DEMAdhesion();

    //! init contact handler
    void Init();

    //! setup contact handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
        const std::shared_ptr<PARTICLEINTERACTION::InteractionWriter> particleinteractionwriter,
        const std::shared_ptr<PARTICLEINTERACTION::DEMNeighborPairs> neighborpairs,
        const std::shared_ptr<PARTICLEINTERACTION::DEMHistoryPairs> historypairs,
        const double& k_normal);

    //! get adhesion distance
    inline double GetAdhesionDistance() const { return adhesion_distance_; };

    //! add adhesion contribution to force field
    void AddForceContribution();

   private:
    //! init adhesion law handler
    void InitAdhesionLawHandler();

    //! init adhesion surface energy handler
    void InitAdhesionSurfaceEnergyHandler();

    //! setup particle interaction writer
    void SetupParticleInteractionWriter();

    //! evaluate particle adhesion contribution
    void EvaluateParticleAdhesion();

    //! evaluate particle-wall adhesion contribution
    void EvaluateParticleWallAdhesion();

    //! discrete element method specific parameter list
    const Teuchos::ParameterList& params_dem_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! interface to particle wall handler
    std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface_;

    //! particle interaction writer
    std::shared_ptr<PARTICLEINTERACTION::InteractionWriter> particleinteractionwriter_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::DEMNeighborPairs> neighborpairs_;

    //! history pair handler
    std::shared_ptr<PARTICLEINTERACTION::DEMHistoryPairs> historypairs_;

    //! adhesion law handler
    std::unique_ptr<PARTICLEINTERACTION::DEMAdhesionLawBase> adhesionlaw_;

    //! adhesion surface energy handler
    std::unique_ptr<PARTICLEINTERACTION::DEMAdhesionSurfaceEnergyBase> adhesionsurfaceenergy_;

    //! adhesion distance
    const double adhesion_distance_;

    //! write particle-wall interaction output
    const bool writeparticlewallinteraction_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
