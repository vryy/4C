/*---------------------------------------------------------------------------*/
/*! \file
\brief boundary particle handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_BOUNDARY_PARTICLE_H
#define PARTICLE_INTERACTION_SPH_BOUNDARY_PARTICLE_H

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
  class SPHNeighborPairs;
}

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHBoundaryParticleBase
  {
   public:
    //! constructor
    explicit SPHBoundaryParticleBase(const Teuchos::ParameterList& params);

    //! virtual destructor
    virtual ~SPHBoundaryParticleBase() = default;

    //! init boundary particle handler
    virtual void Init();

    //! setup boundary particle handler
    virtual void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs);

    //! write restart of boundary particle handler
    virtual void WriteRestart(const int step, const double time) const;

    //! read restart of boundary particle handler
    virtual void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! init boundary particle states
    virtual void InitBoundaryParticleStates(std::vector<double>& gravity) = 0;

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

  class SPHBoundaryParticleAdami : public SPHBoundaryParticleBase
  {
   public:
    //! constructor
    explicit SPHBoundaryParticleAdami(const Teuchos::ParameterList& params);

    //! setup boundary particle handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs) override;

    //! init boundary particle states
    void InitBoundaryParticleStates(std::vector<double>& gravity) override;

   private:
    //! modified states of ghosted boundary particles to refresh
    PARTICLEENGINE::StatesOfTypesToRefresh boundarystatestorefresh_;

    //! contributions of neighboring particles
    std::vector<std::vector<double>> sumj_Wij_;
    std::vector<std::vector<double>> sumj_press_j_Wij_;
    std::vector<std::vector<std::vector<double>>> sumj_dens_j_r_ij_Wij_;
    std::vector<std::vector<std::vector<double>>> sumj_vel_j_Wij_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
