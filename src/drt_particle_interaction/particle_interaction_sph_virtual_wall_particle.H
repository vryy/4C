/*---------------------------------------------------------------------------*/
/*! \file
\brief virtual wall particle handler for smoothed particle hydrodynamics (SPH) interactions
\level 3
*/
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 | definitions                                                               |
 *---------------------------------------------------------------------------*/
#ifndef PARTICLE_INTERACTION_SPH_VIRTUAL_WALL_PARTICLE_H
#define PARTICLE_INTERACTION_SPH_VIRTUAL_WALL_PARTICLE_H

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
  class SPHNeighborPairs;
}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*
 | class declarations                                                        |
 *---------------------------------------------------------------------------*/
namespace PARTICLEINTERACTION
{
  class SPHVirtualWallParticle
  {
   public:
    //! constructor
    explicit SPHVirtualWallParticle(const Teuchos::ParameterList& params);

    //! destructor
    ~SPHVirtualWallParticle() = default;

    //! init virtual wall particle handler
    void Init();

    //! setup virtual wall particle handler
    void Setup(
        const std::shared_ptr<PARTICLEENGINE::ParticleEngineInterface> particleengineinterface,
        const std::shared_ptr<PARTICLEWALL::WallHandlerInterface> particlewallinterface,
        const std::shared_ptr<PARTICLEINTERACTION::SPHKernelBase> kernel,
        const std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs);

    //! write restart of virtual wall particle handler
    void WriteRestart(const int step, const double time) const;

    //! read restart of virtual wall particle handler
    void ReadRestart(const std::shared_ptr<IO::DiscretizationReader> reader);

    //! get reference to relative positions of virtual particles
    inline const std::vector<std::vector<double>>& GetRelativePositionsOfVirtualParticles() const
    {
      return virtualparticles_;
    };

    //! get reference to weighted fluid particle pressure
    inline const std::vector<double>& GetWeightedPressure() const { return weightedpressure_; };

    //! get reference to weighted fluid particle pressure gradient
    inline const std::vector<std::vector<double>>& GetWeightedPressureGradient() const
    {
      return weightedpressuregradient_;
    };

    //! get reference to weighted fluid particle distance vector
    inline const std::vector<std::vector<double>>& GetWeightedDistanceVector() const
    {
      return weighteddistancevector_;
    };

    //! get reference to weighted fluid particle velocity
    inline const std::vector<std::vector<double>>& GetWeightedVelocity() const
    {
      return weightedvelocity_;
    };

    //! init relative positions of virtual particles
    void InitRelativePositionsOfVirtualParticles(const double maxinteractiondistance);

    //! init states at wall contact points
    void InitStatesAtWallContactPoints(std::vector<double>& gravity);

   private:
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

    //! neighbor pair handler
    std::shared_ptr<PARTICLEINTERACTION::SPHNeighborPairs> neighborpairs_;

    //! relative positions of virtual particles
    std::vector<std::vector<double>> virtualparticles_;

    //! weighted fluid particle pressure
    std::vector<double> weightedpressure_;

    //! weighted fluid particle pressure gradient
    std::vector<std::vector<double>> weightedpressuregradient_;

    //! weighted fluid particle distance vector
    std::vector<std::vector<double>> weighteddistancevector_;

    //! weighted fluid particle velocity
    std::vector<std::vector<double>> weightedvelocity_;

    //! set of all fluid particle types
    std::set<PARTICLEENGINE::TypeEnum> allfluidtypes_;

    //! set of integrated fluid particle types
    std::set<PARTICLEENGINE::TypeEnum> intfluidtypes_;
  };

}  // namespace PARTICLEINTERACTION

/*---------------------------------------------------------------------------*/
#endif
