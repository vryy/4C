/*---------------------------------------------------------------------------*/
/*! \file
\brief rigid particle contact handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_RIGID_PARTICLE_CONTACT_H
#define PARTICLE_INTERACTION_SPH_RIGID_PARTICLE_CONTACT_H

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

namespace PARTICLEINTERACTION
{
  class SPHNeighborPairs;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHRigidParticleContactBase
  {
   public:
    //! constructor
    explicit SPHRigidParticleContactBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHRigidParticleContactBase() = default;

    //! init rigid particle contact handler
    virtual void Init();

    //! setup rigid particle contact handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs);

    //! add rigid particle contact contribution to force field
    virtual void AddForceContribution() = 0;

   protected:
    //! smoothed particle hydrodynamics specific parameter list
    const Teuchos::ParameterList& params_sph_;

    //! interface to particle engine
    std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface_;

    //! particle container bundle
    PARTICLEENGINE::ParticleContainerBundleShrdPtr particlecontainerbundle_;

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs_;

    //! set of boundary particle types
    std::set<PARTICLEENGINE::TypeEnum> boundarytypes_;
  };

  class SPHRigidParticleContactElastic : public SPHRigidParticleContactBase
  {
   public:
    //! constructor
    explicit SPHRigidParticleContactElastic(const Teuchos::ParameterList& params);

    //! setup rigid particle contact handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs) override;

    //! add rigid particle contact contribution to force field
    void AddForceContribution() override;

   private:
    //! contact stiffness
    const double stiff_;

    //! contact damping parameter
    const double damp_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
